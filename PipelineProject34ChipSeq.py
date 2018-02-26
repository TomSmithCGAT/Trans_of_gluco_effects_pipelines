'''
PipelineProject34ChipSeq.py - Utility functions for Project 34 ChIP-Seq analysis
==================================================================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Code
----

'''

import os
import pandas as pd
import sqlite3
import copy
import re
import itertools
import glob
import numpy as np
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
import scipy.stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects import pandas2ri

from CGATPipelines.Pipeline import cluster_runnable
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineTracks as PipelineTracks

def database2DataFrame(database, table, index_col):
    # Note: this is currently very specific to generating a df
    # from the context_stats table
    dbh = sqlite3.connect(database)
    df = pd.read_sql("Select * from %(table)s" % locals(), dbh,
                     index_col=index_col)

    df2 = copy.deepcopy(df)

    columns = df2.columns.tolist()
    for ix in df2.index:
        total = df2['total'].ix[ix]
        for column in columns:
            value = df2.ix[ix, column]
            df2.ix[ix, column] = float(value)/total
            df_mean = pd.DataFrame(df2.apply(np.mean, axis=0))
            df_mean.columns = ["mean"]
            df_mean["element"] = df2.columns.tolist()

    threshold = df2.mean() > 0.0001
    df_covered = df2.ix[:, threshold]
    df_covered.drop("total", axis=1, inplace=True)

    return df_covered

def intersectBedFiles(infiles, outfile):
    '''merge :term:`bed` formatted *infiles* by intersection
    and write to *outfile*.

    Only intervals that overlap in all files are retained.
    Interval coordinates are given by the first file in *infiles*.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered starting from 1.
    '''

    if len(infiles) == 1:

        shutil.copyfile(infiles[0], outfile)

    elif len(infiles) == 2:

        if IOTools.isEmpty(infiles[0]) or IOTools.isEmpty(infiles[1]):
            P.touch(outfile)
        else:
            statement = '''
        intersectBed -u -a %s -b %s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        | bgzip > %%(outfile)s 
        ''' % (infiles[0], infiles[1])
            P.run()

    else:

        tmpfile = P.getTempFilename(".")

        # need to merge incrementally
        fn = infiles[0]
        if IOTools.isEmpty(infiles[0]):
            P.touch(outfile)
            return

        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run()

        for fn in infiles[1:]:
            if IOTools.isEmpty(infiles[0]):
                P.touch(outfile)
                os.unlink(tmpfile)
                return

            statement = '''mergeBed -i %(fn)s | intersectBed -u -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run()

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        | bgzip
        > %(outfile)s '''
        P.run()

        os.unlink(tmpfile)

@cluster_runnable
def enrichmentTTest(mapping_database, CONTROLTRACKS, TRACKS, MODIFICATIONS,
                    mapper, outfile):

    df = database2DataFrame(mapping_database, "context_stats", "track")

    features = df.columns

    final_df_t_test = pd.DataFrame()

    for mark in MODIFICATIONS:
        samples = MODIFICATIONS[mark]

        group1 = [str(x) for x in samples if "Saline" in str(x)]
        group2 = [str(x) for x in samples if "Dex" in str(x)]

        control_g1 = [str(x) + ".bwa" for x in CONTROLTRACKS
                      if "Saline" in str(x)][0]
        control_g2 = [str(x) + ".bwa" for x in CONTROLTRACKS
                      if "Dex" in str(x)][0]

        tmp_df_g1 = df.loc[[x + ".bwa" for x in group1]]
        tmp_df_g1_enrichment = (tmp_df_g1/df.loc[control_g1]).applymap(np.log2)

        tmp_df_g2 = df.loc[[x + ".bwa" for x in group2]]
        tmp_df_g2_enrichment = (tmp_df_g2/df.loc[control_g2]).applymap(np.log2)

        p_values = []

        for feature in features:
            p_values.append(scipy.stats.ttest_ind(
                tmp_df_g1_enrichment[feature],
                tmp_df_g2_enrichment[feature],
                axis=0, equal_var=True)[1])

        tmp_df = pd.DataFrame({
            "p-value": p_values,
            "feature": features,
            "mark": str(mark).split("-")[0],
            "Saline_enrichment_mean": tmp_df_g1_enrichment.mean(axis=0),
            "Saline_enrichment_se": tmp_df_g1_enrichment.std(axis=0, ddof=0),
            "Dex_enrichment_mean": tmp_df_g2_enrichment.mean(axis=0),
            "Dex_enrichment_se": tmp_df_g2_enrichment.std(axis=0, ddof=0)
        })

        final_df_t_test = final_df_t_test.append(tmp_df)

    stats = importr('stats')
    final_df_t_test['FDR'] = stats.p_adjust(
        FloatVector(final_df_t_test['p-value']), method='BH')
    final_df_t_test.sort('p-value', inplace=True)

    final_df_t_test.to_csv(outfile, index=False, sep="\t")


@cluster_runnable
def enrichmentBlockedANOVA(mapping_database, CONTROLTRACKS, TRACKS,
                           EXPERIMENTS, mapper, outfile, filename_base):

    df = database2DataFrame(mapping_database, "context_stats", "track")

    features = df.columns

    final_df = pd.DataFrame()

    for sample in EXPERIMENTS:
        control = [str(x) + "." + mapper for x in CONTROLTRACKS
                   if sample.condition in str(x)][0]

        tmp_df = df.loc[[str(x) + "." + mapper for x in EXPERIMENTS[sample]]]
        tmp_df_enrichment = tmp_df/df.loc[control]

        tmp_df_enrichment_log = tmp_df_enrichment.applymap(np.log2)
        final_df = final_df.append(tmp_df_enrichment_log)

    r_df = pandas2ri.py2ri(final_df)

    blockedAnova = R('''
    function(df){

        suppressMessages(library(ggplot2))

        rownames(df) <- gsub(".bwa", "", rownames(df))

        features <- colnames(df)
        features <- features[which(features!="track")]

        df$id_1 = sapply(strsplit(rownames(df), "-"), "[", 1)
        df$id_2 = sapply(strsplit(rownames(df), "-"), "[", 2)
        df$id_3 = sapply(strsplit(rownames(df), "-"), "[", 3)

        blocked_anova_df = data.frame()

        for (modification in c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3")){
            tmp_df = df[grep(modification, rownames(df)), ]

            p = NULL

                for(feature in features){
                    av = aov(tmp_df[[feature]] ~ tmp_df$id_2 + tmp_df$id_3)
                    p = c(p, summary(av)[[1]][["Pr(>F)"]][[1]])
                    }

            p_adj = p.adjust(p, method = "BH")

            tmp_df = data.frame("feature" = features,
                                "mark" = modification,
                                "p_value" = p)

            blocked_anova_df = rbind(blocked_anova_df, tmp_df)
        }

        blocked_anova_df = blocked_anova_df[order(blocked_anova_df$p_value),]
        blocked_anova_df$FDR <- p.adjust(blocked_anova_df$p_value, method = "BH")
        write.table(blocked_anova_df, "%(outfile)s",
                    quote=FALSE, sep="\t", row.names=FALSE)

        candidates = unique(blocked_anova_df[blocked_anova_df$FDR < 0.2,
                            "feature"])

        write.table(df, "%(outfile)s_test.tsv", quote=FALSE, sep="\t")

        l_txt = element_text(size=20)

        for(feature in candidates){
            tmp_df <- data.frame("y"=df[[feature]], "id_3"=df$id_3,
                                 "id_2"=df$id_2, "id_1"=df$id_1)
            p = ggplot(tmp_df, aes(x=id_3, y=y, colour=id_2)) +
                       geom_point(size=5) +
                       facet_grid(.~id_1) +
                       theme(axis.text.x = l_txt, axis.text.y = l_txt,
                             axis.title.x = l_txt, axis.title.y = l_txt,
                             strip.text.x = l_txt, legend.text = l_txt,
                             title = l_txt) +
                       scale_colour_discrete(name="") +
                       ggtitle(feature) +
                       xlab("Litter number") +
                       ylab("log2 enrichment vs. H3 Input")
            ggsave(file=paste0("%(filename_base)s_", feature,".png"),
                   width=10, height=10)
        }

    }''' % locals())

    blockedAnova(r_df)


@cluster_runnable
def enrichmentPlots(mapping_database, CONTROLTRACKS, TRACKS,
                    EXPERIMENTS, mapper, filename_base):

    df = database2DataFrame(mapping_database, "context_stats", "track")

    features = df.columns

    final_df = pd.DataFrame()
    final_df_enrichment = pd.DataFrame()

    featureDict = {
        "cpgisland": "CpG Island",
        "Alu": "SINE - Alu",
        "ID": "SINE - ID",
        "pseudogene": "Pseudogene",
        "protein_coding": "Protein coding"}

    for sample in EXPERIMENTS:
        control = [str(x) + "." + mapper for x in CONTROLTRACKS
                   if sample.condition in str(x)][0]

        tmp_df = df.loc[[str(x) + "." + mapper
                         for x in EXPERIMENTS[sample]]]
        tmp_df_enrichment = tmp_df/df.loc[control]

        tmp_df_enrichment_log = tmp_df_enrichment.applymap(np.log2)
        final_df = final_df.append(tmp_df_enrichment_log)

        array_enrichment_log_se = tmp_df_enrichment_log.std(axis=0, ddof=0)
        array_se = tmp_df.std(axis=0, ddof=0)

        tmp_df = pd.DataFrame({
            "log_enrichment_se": array_enrichment_log_se,
            "mark": sample.mark,
            "condition": sample.condition,
            "log_enrichment_mean": tmp_df_enrichment_log.mean(axis=0),
            "log_enrichment_max": tmp_df_enrichment_log.max(axis=0),
            "log_enrichment_min": tmp_df_enrichment_log.min(axis=0),
            "feature": array_se.index,
            "mean": tmp_df.mean(axis=0),
            "se": tmp_df.std(axis=0, ddof=0)})

        tmp_df['feature'] = [featureDict[x] if x in featureDict else x
                             for x in tmp_df['feature']]

        final_df_enrichment = final_df_enrichment.append(tmp_df)

    final_df_enrichment['condition'] = [x.replace("Saline", "Veh") for x in
                                        final_df_enrichment['condition']]

    final_df_enrichment.reset_index(inplace=True)
    #final_df_enrichment.set_index(
    #    [range(0, len(final_df_enrichment.index))], inplace=True)
    final_df_enrichment.to_csv(filename_base + ".tsv",
                               index=False, sep="\t")

    r_df = pandas2ri.py2ri(final_df)
    r_df_enrichment = pandas2ri.py2ri(final_df_enrichment)

    plotEnrichmentsPCA = R('''
    function(df){

        suppressMessages(library(ggplot2))
        suppressMessages(library(grid))
        suppressMessages(library(scales))

        l_txt = element_text(size=20)
        m_txt = element_text(size=15)
        s_txt = element_text(size=10)

        pca <- prcomp(df, center = TRUE, scale=TRUE)

        variance = pca$sdev^2
        variance_explained = round(variance/sum(variance), 5)

        variance_df = data.frame("Variance_explained" = variance_explained,
                                 "PC" = seq(1, length(variance)))

        make_plot_variance <- function(variance_df){
            p_variance = ggplot(variance_df, aes(x=PC, y=Variance_explained))+
            geom_point(size=5)+
            geom_line()+
            theme_bw()+
            ylab("Variance explained (%%)")+
            theme(axis.text.x = l_txt,
                axis.title.y = l_txt,
                axis.title.x = l_txt,
                axis.text.y = l_txt)
           return(p_variance)
           }

        p_variance = make_plot_variance(variance_df)
        p_variance = p_variance + ggtitle("PC variance explained for all modifications")

        ggsave("%(filename_base)s_pca_all_variance.png", width=10, height=10)

        make_PC_df <- function(pca){

            PCs_df = data.frame(pca$x)
            rownames(PCs_df) = gsub(".bwa", "", rownames(PCs_df))

            PCs_df$id_1 = sapply(strsplit(rownames(PCs_df), "-"), "[", 1)
            PCs_df$id_2 = sapply(strsplit(rownames(PCs_df), "-"), "[", 2)
            PCs_df$id_3 = sapply(strsplit(rownames(PCs_df), "-"), "[", 3)

            return(PCs_df)}

        PCs_df = make_PC_df(pca)

        PCs_df$id_3[PCs_df$id_3 == 3] <- 1
        PCs_df$id_3[PCs_df$id_3 == 7] <- 2
        PCs_df$id_3[PCs_df$id_3 == 9] <- 3

        write.table(PCs_df,
                    "/ifs/projects/proj034/ChIP_Seq/full/project_pipeline/test.tsv",
                    sep="\t")

        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        three_tones = c(cbPalette[c(3,4,6,8)])

        point <- geom_point(colour="black",
                            aes(shape=id_3, size=id_2,fill=id_1))
        fill <- scale_fill_manual(values=three_tones,
                                  name=guide_legend(title='Modification'))
        shape <- scale_shape_manual(values=c(21,22,24,25),
                                    name=guide_legend(title='Replicate'))
        size <- scale_size_manual(values=c(15,5),
                                  name=guide_legend(title='Treatment'))
        my_theme <- theme(
             axis.text.x = l_txt, axis.text.y = l_txt, title = l_txt,
             legend.text = l_txt, legend.title = l_txt,
             legend.key.size = unit(7, "mm"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

        fill_guide <- guides(fill = guide_legend(override.aes=list(
                             colour = three_tones, size=10)))
        fill_shape <- guides(shape = guide_legend(override.aes=list(size=10)))

        p_pca1 = ggplot(PCs_df, aes(x=PC1, y=PC2)) +
        point + fill + shape + size + theme_bw() + my_theme + fill_guide + fill_shape +
        xlab(paste0('PC1 (Variance explained = ' ,
                   round(100 * variance_explained[1], 1), '%%)')) +
        ylab(paste0('PC2 (Variance explained = ' ,
                   round(100 * variance_explained[2], 1), '%%)'))

        ggsave("%(filename_base)s_pca_all_PC1_PC2.png", width=10, height=10)

        p_pca2 = ggplot(PCs_df, aes(x=PC3, y=PC4)) +
        point + fill + shape + size + theme_bw() + my_theme + fill_guide + fill_shape +
        xlab(paste0('PC3 (Variance explained = ' ,
                   round(100 * variance_explained[3], 1), '%%)')) +
        ylab(paste0('PC4 (Variance explained = ' ,
                   round(100 * variance_explained[4], 1), '%%)'))

        #ggsave("%(filename_base)s_pca_all_PC3_PC4.png", width=10, height=10)


        point_mod <- geom_point(size=15,  aes(shape=id_3, colour=id_2))
        shape_mod <- scale_shape_manual(name=guide_legend(title='Replicate'),
                                        values=c(15,16,17,18))
        colour_mod <-scale_fill_discrete(name=guide_legend(title='Treatment'))

        for(modification in c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3")){

          pca_mod <- prcomp(df[grep(modification, rownames(df)),],
                            center = TRUE, scale=TRUE)

          variance_mod = pca_mod$sdev^2
          variance_mod_explained = round(variance_mod/sum(variance_mod), 5)

          variance_df = data.frame(
           "Variance_explained" = variance_mod_explained,
           "PC" = seq(1, length(variance_mod)))

          p_variance = make_plot_variance(variance_df)

          p_variance = p_variance + ggtitle(paste0(
            "variance explained in PCs for ", modification))

          ggsave(paste0("%(filename_base)s_pca_",modification,"_variance.png"),
                 width=10, height=10)

          PCs_mod_df = make_PC_df(pca_mod)

          p_pca1 = ggplot(PCs_mod_df, aes(x=PC1, y=PC2)) +
          point_mod + shape_mod + colour_mod + theme_bw() + my_theme +
          xlab(paste0('PC1 (Variance explained = ',
                round(100 * variance_mod_explained[1], 1), '%%')) +
          ylab(paste0('PC2 (Variance explained = ' ,
                      round(100 * variance_mod_explained[2], 1), '%%'))+
          ggtitle(modification)

          ggsave(paste0("%(filename_base)s_pca_", modification , "PC1_PC2.png"),
                 width=10, height=10)

          p_pca1 = ggplot(PCs_mod_df, aes(x=PC3, y=PC4)) +
          point_mod + shape_mod + colour_mod + theme_bw() + my_theme +
          xlab(paste0('PC3 (Variance explained = ',
                round(100 * variance_mod_explained[3], 1), '%%')) +
          ylab(paste0('PC4 (Variance explained = ' ,
                      round(100 * variance_mod_explained[4], 1), '%%'))+
          ggtitle(modification)

          ggsave(paste0("%(filename_base)s_pca_", modification , "PC3_PC4.png"),
                 width=10, height=10)

        }

    }''' % locals())

    plotEnrichmentsPCA(r_df)

    plotEnrichmentsBar = R('''
    function(df){
        suppressMessages(library(grid))
        suppressMessages(library(ggplot2))

        top_enrichments = unique(df[2**df$log_enrichment_mean > 1.2, "feature"])
        top_depletions = unique(df[2**df$log_enrichment_mean <0.8, "feature"])

        top_both = c(as.character(top_enrichments),
                     as.character(top_depletions))

        limits <- aes(ymax = log_enrichment_mean + log_enrichment_se,
                      ymin = log_enrichment_mean - log_enrichment_se)
        #limits <- aes(ymax = log_enrichment_max, ymin = log_enrichment_min)
        dodge <- position_dodge(width=0.9)

        df = df[df$feature %%in%% top_both,]
        write.table(df, "./test.tsv", sep="\t")
        df$feature = factor(df$feature, levels=c(
            "SINE - Alu", "Protein coding", "rRNA", "Pseudogene", "CpG Island"))
        df$mark = factor(df$mark, levels=c(
            "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3"))

        df = df[!is.na(df$feature),]
        df = df[!is.na(df$mark),]

        l_txt = element_text(size=25)

        p = ggplot(df, aes(x=as.factor(feature),
                           y=as.numeric(log_enrichment_mean),
                           group=as.factor(condition))) +
        geom_bar(aes(fill=as.factor(condition)),
                     position=dodge, stat="identity") +

        geom_errorbar(limits, position=dodge, width=0.25) +
        xlab("") + ylab("Enrichment vs. H3-INPUT (log2)") +
        theme_bw() +
        theme(axis.text.x=element_text(size=20, angle=90, vjust=0.5, hjust=1),
              axis.text.y=l_txt, axis.title.x=l_txt, axis.title.y=l_txt,
              legend.text=l_txt, legend.title=l_txt, strip.text.x = l_txt,
              legend.key.size=unit(1, "cm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_fill_discrete(name="Treatment") +
        facet_grid(. ~ mark) +
        scale_fill_manual(name="", values=c("#9e66ab", "#f9a65a")) +
        scale_y_continuous(breaks=seq(-1,20,0.5))

        ggsave("%(filename_base)s.png", width=10, height=10)

    }''' % locals())

    plotEnrichmentsBar(r_df_enrichment)


@cluster_runnable
def runMMDiff(infiles, sample_ids, tsv_outfile, plot_outfile):

    header = ",".join(("sampleID", "Tissue", "Factor", "Condition",
                       "Replicate", "bamReads", "bamControl", "Peaks",
                       "PeakCaller"))

    tmp_info_file = P.getTempFilename(".")

    with IOTools.openFile(tmp_info_file, "w") as outf:
        outf.write(header + "\n")
        for infile in infiles:
            sample = sample_ids[infiles.index(infile)]
            tissue = "WT"

            treatment, marker, replicate = sample.split("-")

            if treatment == "Saline":
                condition = 0
            elif treatment == "Dex":
                condition = 1

            bam = "bams/%s.call.bam" % sample
            bamControl = "bams/%s-H3-7.call.bam" % treatment
            Peaks = infile
            caller = "bed"

            outf.write(",".join(map(str, (
                sample, tissue, marker, condition, replicate, bam,
                bamControl, Peaks, caller))) + "\n")

    outf.close()

    gr = importr('grDevices')

    MMDiff = R('''function(){
    suppressMessages(library(MMDiff))
    suppressMessages(library(DiffBind))
    Cfp1 <- dba(sampleSheet="%(tmp_info_file)s", minOverlap=2, bCorPlot=FALSE)

    Cfp1Profiles <- getPeakProfiles(
        Cfp1, draw.on=FALSE, Peaks=dba.peakset(Cfp1,bRetrieve=TRUE),
        bin.length=50, run.parallel=FALSE, keep.extra=TRUE)

    Cfp1Norm <- getNormFactors(Cfp1Profiles, method = "DESeq")

    Cfp1Dists <- compHistDists(Cfp1Norm, method='MMD', overWrite=FALSE,
                               NormMethod='DESeq', run.parallel=FALSE )

    sample_ids = Cfp1Dists$sample$sampleID
    group.1 <- seq(1,length(sample_ids))[grepl("Saline", sample_ids)]
    group.2 <- seq(1,length(sample_ids))[grepl("Dex", sample_ids)]

    Cfp1Pvals <- detPeakPvals(Cfp1Dists, group1=group.1, group2=group.2,
    name1='Saline', name2='Dex', method='MMD')

    idxMMD <- which(as.data.frame(Cfp1Pvals$MD$Pvals$MMD)$combined<0.05)

    png(filename="%(plot_outfile)s")
    plotHistDists(Cfp1Pvals, group1=group.1, group2=group.2,
                  thresh = 0.1, method='MMD')
    dev.off()

    adj_p_value = as.data.frame(Cfp1Pvals$MD$Pvals$MMD)$combined
    peaks = rownames(as.data.frame(Cfp1Pvals$MD$Pvals$MMD))

    df = data.frame("peak" = peaks, adjusted_p_value = adj_p_value)
    write.table(df, file="%(tsv_outfile)s", quote=FALSE,
                row.names=FALSE, sep="\t")
    write.table(as.data.frame(Cfp1Pvals$MD$DISTS$MMD),
                file="%(tsv_outfile)s_all.tsv", quote=FALSE, sep="\t")

    }''' % locals())
    
    MMDiff()

    #os.unlink(tmp_info_file)


def plotGeneProfiles(infiles, outfile):
    ''' take all geneprofile flatfiles and generate plots
    for tss profile and geneprofile:
    * All modifications per sample overlaid for each profile
    * All samples per modification overlaid for each profile

    Note: the infile name is the basename
    For example Saline-H3INPUT-7_geneset_coding_exons.tsv.gz
    --> Saline-H3INPUT-7_geneset_coding_exons.tsv.gz.geneprofile.matrix.tsv.gz
    '''
    pandas2ri.activate()

    def getEnrichment(inf_ip, df_input1, df_input2):
        '''
        Get enrichment from area in ip and input
        '''
        df = pd.read_table(inf_ip)

        df['input'] = np.divide((df_input1['area'] + df_input2['area']), 2)

        # use the average of the two inputs
        df['enrichment'] = np.divide(df['area'].astype(float), df['input'])

        return df

    profiles = ["tssprofile", "geneprofile"]

    samples = set()
    for infile in infiles:
        sample = os.path.basename(infile).split("_")[0]
        condition = sample.split("-")[0]
        replicate = sample.split("-")[2]
        samples.update(("%s-%s" % (condition, replicate),))

    with IOTools.openFile(outfile, "w") as log:
        for profile in profiles:
            profile_infiles = ["%s.%s.matrix.tsv.gz" % (x, profile) for x in infiles]
            input_infiles = [x for x in profile_infiles if "H3-" in x]

            df_input1 = pd.read_table(input_infiles[0])
            df_input2 = pd.read_table(input_infiles[1])

            profile_infiles = [x for x in profile_infiles if "H3-" not in x]

            df = pd.DataFrame()

            for inf in profile_infiles:
                sample = os.path.basename(inf).split("_")[0]
                condition = sample.split("-")[0]
                ip = sample.split("-")[1]
                replicate = sample.split("-")[2]

                tmp_df = getEnrichment(inf, df_input1, df_input2)
                tmp_df['ip'] = ip
                tmp_df['condition'] = condition
                tmp_df['replicate'] = replicate
                tmp_df['region'] = [x.upper() for x in tmp_df['region']]

                df = pd.concat([df, tmp_df])

            df['condition'] = [x.replace("Saline", "Veh") for x in
                               df['condition']]

            df.index = range(0, len(df))

            plot_outfile = P.snip(outfile, ".log") + "_%s_combined.png" % profile
            plot_outfile_per_mark = (P.snip(outfile, ".log") +
                                     "_%s_combined_by_mark.png" % profile)

            plotEnrichments = R('''function(df, profile){

            suppressMessages(library(ggplot2))
            suppressMessages(library(RColorBrewer))
            
            cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

            df$ip <- factor(df$ip, levels=c("H3K4me1", "H3K4me3",
                                            "H3K9me3", "H3K27me3"))
            df$enrichment = as.numeric(df$enrichment)
            df$bin = as.numeric(df$bin)

            m_txt = element_text(size=15)
            s_txt = element_text(size=10)

            p = ggplot(df, aes(colour=ip, group=interaction(replicate, ip))) +
            scale_colour_manual(name="", values=cbPalette) +
            xlab("") + ylab("IP Enrichment") +
            theme_bw() +
            theme(
            axis.text.y=m_txt,
            axis.text.x=s_txt,
            axis.title=m_txt,
            strip.text=m_txt,
            legend.text=m_txt,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
            theme(aspect.ratio=1)

            if (profile=="geneprofile"){
            p = p +
            geom_vline(xintercept=1000, linetype=2, size=0.25) +
            geom_vline(xintercept=2000, linetype=2, size=0.25) +
            geom_line(aes(x=bin, y=enrichment), size=0.25) +
            #scale_y_continuous(breaks=seq(0.8, 1.8, 0.2), limits=c(0.8, 1.8)) +
            scale_x_continuous(breaks=c(500,1500,2550),
                               labels=c("Upstream", "Exons", "Downstream")) +
            theme(
            axis.ticks.x=element_blank()) +
            facet_grid(.~condition) +
            guides(colour = guide_legend(override.aes = list(size=2)))

            ggsave("%(plot_outfile)s")#, width=7, height=4)

            p = p + facet_grid(ip~condition, scales="free") + 
            geom_line(aes(x=bin, y=enrichment), size=0.5) +
            theme(legend.position="none")

            ggsave("%(plot_outfile_per_mark)s")#, width=6, height=12)

            }

            else if(profile=="tssprofile"){
            p = p +
            geom_vline(xintercept=0, linetype=2, size=0.25) +
            geom_line(aes(x=region_bin-3000, y=enrichment), size=0.25) +
            #scale_y_continuous(breaks=seq(0.8, 1.8, 0.2), limits=c(0.8, 1.8)) +
            scale_x_continuous(breaks=c(-3000, -1500, 0, 1500, 3000)) +
            theme(
            axis.text.x=element_text(angle=90, size=20, vjust=0.5, hjust=1)) +
            facet_grid(region~condition, scales="free_x") +
            guides(colour = guide_legend(override.aes = list(size=2)))

            ggsave("%(plot_outfile)s", width=7, height=8)
            }

            }''' % locals())

            plotEnrichments(df, profile)
            log.write("plot generated at %s\n" % plot_outfile)

            for sample in samples:
                cond = sample.split("-")[0]
                rep = sample.split("-")[1]
                pattern = ".*/%s-.*-%s_geneset_coding_exons.tsv.gz" % (cond, rep)
                sample_infiles = [x for x in profile_infiles if re.match(pattern, x)]

                df = pd.DataFrame()

                for inf in sample_infiles:
                    ip = os.path.basename(inf).split("_")[0].split("-")[1]

                    tmp_df = getEnrichment(inf, df_input1, df_input2)
                    tmp_df['condition'] = cond
                    tmp_df['replicate'] = rep
                    tmp_df['ip'] = ip
                    tmp_df['region'] = [x.upper() for x in tmp_df['region']]

                    df = pd.concat([df, tmp_df])

                plot_outfile = P.snip(outfile, ".log") + "_%s_%s.png" % (sample, profile)

                plotEnrichmentsSample = R('''function(df, profile){
                suppressMessages(library(ggplot2))
                df$ip <- factor(df$ip, levels=c("H3K4me1", "H3K4me3",
                                                "H3K9me3", "H3K27me3"))
                df$enrichment = as.numeric(df$enrichment)
                df$bin = as.numeric(df$bin)

                m_txt = element_text(size=15)

                p = ggplot(df, aes(colour=ip, group=interaction(replicate, ip))) +
                scale_colour_discrete(name="") +
                xlab("") + ylab("IP Enrichment") +
                theme_bw() +
                theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text=m_txt,
                axis.title=m_txt,
                strip.text=m_txt,
                legend.text=m_txt) +
                theme(aspect.ratio=1)

                if (profile=="geneprofile"){
                p = p +
                geom_vline(xintercept=1000, linetype=2, size=0.25) +
                geom_vline(xintercept=2000, linetype=2, size=0.25) +
                scale_x_continuous(breaks=c(500,1500,2500),
                                   labels=c("Upstream", "Exons", "Downstream")) +
                scale_y_continuous(breaks=seq(0.8, 1.8, 0.2), limits=c(0.8, 1.8)) +
                theme(
                axis.ticks.x=element_blank()) +
                guides(colour = guide_legend(override.aes = list(size=2))) +
                geom_line(aes(x=bin, y=enrichment), size=0.25)


                ggsave("%(plot_outfile)s", width=8, height=5)
                }

                else if(profile=="tssprofile"){
                p = p +
                geom_vline(xintercept=0, linetype=2, size=0.25) +
                scale_x_continuous(breaks=c(-3000, -1500, 0, 1500, 3000)) +
                scale_y_continuous(breaks=seq(0.8, 1.8, 0.2), limits=c(0.8, 1.8)) +
                theme(
                axis.text.x=element_text(angle=90, size=15, vjust=0.5, hjust=1)) +
                facet_grid(.~region, scales="free_x") +
                guides(colour = guide_legend(override.aes = list(size=2))) +
                geom_line(aes(x=region_bin-3000, y=enrichment), size=0.25)

                ggsave("%(plot_outfile)s", width=11, height=5)
                }

                }''' % locals())

                plotEnrichmentsSample(df, profile)
                log.write("plot generated at %s\n" % plot_outfile)


######################################################################
# move to module file
######################################################################


def intersectPeaks(infiles, samples,
                   intersection_outfile, all_peaks_outfile,
                   plot_outfile_counts, plot_outfile_fraction):
    '''use intersectBed within each group (e.g treatment + histone mark)
    and identify intersections

    Note: infiles and samples must be in the same order'''

    intersect_values = {}

    # first get the numbers of peaks in a single file
    with IOTools.openFile(all_peaks_outfile, "w") as outf:
        empty_files = 0
        for infile in infiles:
            n = 0
            with IOTools.openFile(infile, "r") as inf:
                for line in inf:
                    if line[0:7] == "chrNull":
                        empty_files += 1
                        break
                    else:
                        outf.write(line)
                        n += 1

            intersect_values[infile] = n

    outf.close()
    assert empty_files != len(infiles), "no peaks in any of the infiles"

    for pair in itertools.combinations(infiles, 2):
        inf1, inf2 = pair

        # check if either sample has no peaks, if so, no need to intersect
        n = 0
        if intersect_values[inf1] == 0 or intersect_values[inf2] == 0:
            pass
        else:
            tmp_file = P.getTempFilename(shared=True) + ".gz"
            P34ChipSeq.intersectBedFiles(pair, tmp_file)

            with IOTools.openFile(tmp_file, "r") as tmp_inf:
                for line in tmp_inf:
                    n += 1

        os.unlink(tmp_file)
        intersect_values[inf1 + ":" + inf2] = n

    # if any sample has no peaks, or any combination has no peaks,
    # no need to intersect all
    counts = [intersect_values[infile] for infile in infiles]
    counts.extend([intersect_values[x[0] + ":" + x[1]]
                   for x in itertools.combinations(infiles, 2)])

    n = 0
    if min(counts) == 0:
        P.touch(intersection_outfile)
    else:
        P34ChipSeq.intersectBedFiles(infiles, intersection_outfile)
        with IOTools.openFile(intersection_outfile, "r") as inf:
            for line in inf:
                n += 1

    intersect_values["all"] = n

    if len(infiles) == 3:
        # use intersect values to generate values for venn diagram
        subsets = []
        s1, s2, s3 = infiles
        subsets.append(intersect_values[s1] - intersect_values[s1 + ":" + s2] -
                       intersect_values[s1 + ":" + s3] + intersect_values["all"])
        subsets.append(intersect_values[s2] - intersect_values[s1 + ":" + s2] -
                       intersect_values[s2 + ":" + s3] + intersect_values["all"])
        subsets.append(intersect_values[s1 + ":" + s2] - intersect_values["all"])
        subsets.append(intersect_values[s3] - intersect_values[s1 + ":" + s3] -
                       intersect_values[s2 + ":" + s3] + intersect_values["all"])
        subsets.append(intersect_values[s1 + ":" + s3] - intersect_values["all"])
        subsets.append(intersect_values[s2 + ":" + s3] - intersect_values["all"])

        # want to overlap all three even if no peaks in full intersection
        # use a dummy value (1) and replace in plot label
        if intersect_values["all"] == 0:
            subsets.append(1)
            replace = 1
        else:
            subsets.append(intersect_values["all"])
            replace = 0

        mark = s1.split("-")[1]

        labels = [x.split("-")[0] for x in samples]
        for ix in range(0, len(labels)):
            replicate_number = ix + 1
            labels[ix] = "%s-%i" % (labels[ix].split("-")[0], replicate_number)

        v = venn3(subsets=subsets, set_labels=labels)
        if replace:
            try:
                v.get_label_by_id('111').set_text('0')
            except:
                print("not possible to replace text in venn diagram")
        plt.title("Peak overlap for %s" % mark)
        plt.savefig(plot_outfile_counts)
        plt.close()

        if intersect_values["all"] == 0:
            subsets = [round(np.divide(float(x), sum(subsets)-1), 2)
                       for x in subsets]
        else:
            subsets = [round(np.divide(float(x), sum(subsets)), 2)
                       for x in subsets]

        v = venn3(subsets=subsets, set_labels=labels)
        if replace:
            try:
                v.get_label_by_id('111').set_text('0')
            except:
                print("not possible to replace text in venn diagram")
        plt.title("Peak overlap for %s" % mark)
        plt.savefig(plot_outfile_fraction)
        plt.close()

    elif len(infiles) == 2:
        print("\n\n\n\n Only 2 infiles: %s \n\n\n\n" % ",".join(infiles))
        print(infiles)
        print(intersect_values)
        # use intersect values to generate values for venn diagram
        subsets = []
        s1, s2 = infiles
        subsets.append(intersect_values[s1] - intersect_values["all"])
        subsets.append(intersect_values[s2] - intersect_values["all"])
        subsets.append(intersect_values["all"])

        # want to overlap all three even if no peaks in full intersection
        # use a dummy value (1) and replace in plot label
        if intersect_values["all"] == 0:
            subsets.append(1)
            replace = 1
        else:
            subsets.append(intersect_values["all"])
            replace = 0

        mark = s1.split("-")[1]

        labels = [x.split("-")[0] for x in samples]
        for ix in range(0, len(labels)):
            replicate_number = ix + 1
            labels[ix] = "%s-%i" % (labels[ix].split("-")[0], replicate_number)

        v = venn2(subsets=subsets, set_labels=labels)
        if replace:
            try:
                v.get_label_by_id('111').set_text('0')
            except:
                print("not possible to replace text in venn diagram")
        plt.title("Peak overlap for %s" % mark)
        plt.savefig(plot_outfile_counts)
        plt.close()

        if intersect_values["all"] == 0:
            subsets = [round(np.divide(float(x), sum(subsets)-1), 2)
                       for x in subsets]
        else:
            subsets = [round(np.divide(float(x), sum(subsets)), 2)
                       for x in subsets]

        v = venn2(subsets=subsets, set_labels=labels)
        if replace:
            try:
                v.get_label_by_id('111').set_text('0')
            except:
                print("not possible to replace text in venn diagram")
        plt.title("Peak overlap for %s" % mark)
        plt.savefig(plot_outfile_fraction)
        plt.close()


def intersectPeaksHeatmap(infiles, outfile, plotfile):
    ''' get the intersection of peaks in all samples'''

    inf_counts = {}

    with IOTools.openFile(outfile, "w") as outf:

        outf.write("%s\n" % "\t".join(
            ("sample1", "sample2", "count", "fraction")))

        # first get the numbers of peaks in a single file
        empty_files = 0
        for infile in infiles:

            sample = os.path.basename(infile).split("_")[0]

            print(infile)
            n = 0
            with IOTools.openFile(infile, "r") as inf:
                for line in inf:
                    if line[0:7] == "chrNull":
                        empty_files += 1
                        break
                    else:
                        n += 1

            inf_counts[infile] = n
            outf.write("%s\n" % "\t".join(map(str, (sample, sample, n, 1.0))))

        for pair in itertools.combinations(infiles, 2):
            print(pair)
            inf1, inf2 = pair
            sample1 = os.path.basename(inf1).split("_")[0]
            sample2 = os.path.basename(inf2).split("_")[0]

            # check if either sample has no peaks, if so, no need to intersect
            n = 0
            if inf_counts[inf1] == 0 or inf_counts[inf2] == 0:
                pass
            else:
                tmp_file = P.getTempFilename(shared=True) + ".gz"
                intersectBedFiles(pair, tmp_file)

                with IOTools.openFile(tmp_file, "r") as tmp_inf:
                    for line in tmp_inf:
                        n += 1

            os.unlink(tmp_file)

            fraction1 = float(n)/(inf_counts[inf1])
            fraction2 = float(n)/(inf_counts[inf2])

            outf.write("%s\n" % "\t".join(
                map(str, (sample1, sample2, n, fraction1))))
            outf.write("%s\n" % "\t".join(
                map(str, (sample2, sample1, n, fraction2))))

    plotIntersectionHeatmap = R('''
    df = read.table("%(outfile)s", sep="\t", header=T)
    library(ggplot2)
    m_txt = element_text(size=10)
    m_txt_90 = element_text(size=10, angle=90, vjust=0.5, hjust=1)

    p = ggplot(df, aes(sample1, sample2, fill=100*fraction)) +
    geom_tile() +
    geom_text(aes(label=count), size=10) +
    scale_fill_gradient(name="Intersection (%%)",
                        low="yellow", high="dodgerblue4") +
    theme(axis.text.x = m_txt_90, axis.text.y = m_txt,
          legend.text = m_txt, legend.title = m_txt,
          aspect.ratio=1) +
    xlab("") + ylab("")

    ggsave("%(plotfile)s", width=20, height=20)
    ''' % locals())

    plotIntersectionHeatmap


def plotHistoneRNACorrelation(infiles, outfile, nreads_dir, alpha,
                              DESeq2_fold_inf, suffix, nreads_dict,
                              promoter_pad):

    ''' correlate the sRNA-Seq fold changes and fold change in
    piRNA/miRNA promoter histone enrichment fold changes '''

    pandas2ri.activate()

    def histone2df(infile, count_id):
        histone_df = pd.read_table(
            infile, sep="\t", header=None, usecols=[3, 6])
        histone_df.columns = ["Geneid", count_id]
        return histone_df

    def getHistoneFoldChange(infiles, nreads_dict, suffix,
                             feature_length=1000):

        sample_id = os.path.basename(infiles[0]).replace(suffix, "")
        histone_df = histone2df(infiles[0], sample_id)

        for infile in infiles[1:]:
            sample_id = os.path.basename(infile).replace(suffix, "")
            tmp_df = histone2df(infile, sample_id)
            histone_df = pd.merge(histone_df, tmp_df,
                                  left_on="Geneid", right_on="Geneid")

        histone_df.set_index("Geneid", inplace=True)

        # subset to covered promoters
        histone_df = histone_df[
            histone_df.apply(axis=1, func=max) > feature_length/50.0]
        histone_df = histone_df[
            histone_df.apply(axis=1, func=min) > feature_length/1000.0]

        for sample in histone_df.columns:
            nreads = nreads = float(nreads_dict[sample])
            histone_df[sample] = histone_df[sample] / (nreads/1000000)

        # calculate enrichment
        histone_enrichment_df = pd.DataFrame()

        for sample in histone_df.columns:
            if "H3-" in sample:
                pass
            else:
                tmp_df = histone_df[[sample, "Dex-H3-7", "Saline-H3-7"]]
                histone_enrichment_df[sample] = tmp_df.apply(
                        axis=1, func=lambda x: np.divide(
                            float(x[0]+1), np.divide(x[1]+x[2], 2)+1))

        histone_enrichment_df.reset_index(inplace=True)

        histone_enrichment_df = pd.melt(histone_enrichment_df, id_vars="Geneid")

        # extract sample info from name
        sample_info_df = pd.DataFrame.from_records(
            [x.split("-") for x in histone_enrichment_df['variable'].tolist()])
        sample_info_df.columns = ["condition", "mark", "replicate"]
        histone_enrichment_df = pd.concat(
            [histone_enrichment_df, sample_info_df], axis=1)


        # aggregate replicates
        histone_enrichment_agg_sum_df = histone_enrichment_df[
            ["Geneid", "condition", "mark", "value"]].groupby(
                ["Geneid", "condition", "mark"]).sum()
        histone_enrichment_agg_sd_df = histone_enrichment_df[
            ["Geneid", "condition", "mark", "value"]].groupby(
                ["Geneid", "condition", "mark"]).apply(lambda x: np.std(x))

        histone_enrichment_agg_sum_df.reset_index(inplace=True)
        histone_enrichment_agg_sd_df.reset_index(inplace=True)

        # get fold change
        histone_fold = pd.DataFrame()
        for histone_mark in set(histone_enrichment_agg_sum_df['mark']):
            tmp_df = histone_enrichment_agg_sum_df.loc[
                histone_enrichment_agg_sum_df['mark'] == histone_mark]
            tmp_df = tmp_df.pivot(index="Geneid", columns="condition", values="value")
            tmp_df['fold'] = tmp_df.apply(axis=1, func=lambda x: np.log2(
                np.divide(x['Dex'], x['Saline'])))
            tmp_df.columns = ["_".join((x, histone_mark)) for x in tmp_df.columns]
            histone_fold = pd.concat((histone_fold, tmp_df), axis=1)

        return histone_fold

    def mergeHistonesRNASeq(histone_fold, DESeq2_fold):
        ''' merge the histone fold change and DESeq2 fold change '''
        fold_df = pd.read_table(DESeq2_fold, sep="\t", comment="#")
        fold_df.set_index("test_id", inplace=True)

        fold_df = fold_df[["l2fold", "transformed_l2fold",
                           "control_mean", "p_value"]]

        fold_df.index.name = "Geneid"

        final_df = histone_fold.merge(
            fold_df, left_index=True, right_index=True)

        return final_df

    histone_fold = getHistoneFoldChange(
        infiles, nreads_dict, suffix,
        feature_length=promoter_pad)

    final_df = mergeHistonesRNASeq(histone_fold, DESeq2_fold_inf)

    lm_summary_outfile = P.snip(outfile, ".png") + "_lm_summary.txt"

    final_df.to_csv(P.snip(outfile, ".png") + ".tsv", sep="\t")

    plotCorrelations = R('''function(df){

    library(ggplot2)
    library(gridExtra)
    library(cowplot)

    sink(file="%(lm_summary_outfile)s")

    fit = lm(formula = transformed_l2fold ~ fold_H3K27me3 +
             fold_H3K4me1 + fold_H3K4me3 + fold_H3K9me3, data = df)

    print(summary(fit))

    m_txt = element_text(size=20)

    plots = NULL

    for (col in c("fold_H3K4me1", "fold_H3K4me3", "fold_H3K9me3", "fold_H3K27me3")){

        print(col)

        title = gsub("fold_", "", col)
        tmp_df = df[,c(col, "transformed_l2fold")]
        tmp_df = tmp_df[is.finite(tmp_df[["transformed_l2fold"]]),]

        print(cor(tmp_df, method="spearman"))
        print(cor(tmp_df, method="pearson"))

        fit = lm(formula = formula(paste0("transformed_l2fold ~", col)),
                 data = tmp_df)

        print(summary(fit))

        colnames(tmp_df) <- c("histone", "fold")

        print(head(tmp_df))

        p = ggplot(tmp_df, aes(x=histone, y=fold)) +
        geom_point(alpha=%(alpha)s) +
        ggtitle(title) +
        geom_smooth(method="lm") +
        xlab("Fold change histone enrichment") +
        ylab("Expression fold change") +
        theme_cowplot() +
        theme(
        axis.text.x=m_txt,
        axis.text.y=m_txt,
        axis.title.x=m_txt,
        axis.title.y=m_txt,
        plot.title=m_txt,
        aspect.ratio=1) +
        geom_vline(xintercept=0, linetype=2, colour="grey60", size=0.5) +
        geom_hline(yintercept=0, linetype=2, colour="grey60", size=0.5)

        plots[[title]] = p
    }

    sink()

    png("%(outfile)s", width=30, height=30, units="cm", res=400)
    do.call("grid.arrange", c(c(plots), ncol=2))
    dev.off()

    }''' % locals())

    plotCorrelations(final_df)


def plotChipSeqheatmap(infile, outfile, nreads_dict):
    ''' make a heatmap + dendograms for the ChIP-Seq enrichment scores'''

    pandas2ri.activate()

    df = pd.read_table(infile, sep="\t")
    df_pivot = pd.pivot_table(df, index=["start", "contig", "length"],
                              columns="sample", values="count")

    df_pivot = df_pivot[df_pivot['Dex-H3-7'] > 10]
    df_pivot = df_pivot[df_pivot['Saline-H3-7'] > 10]

    for sample in df_pivot.columns:
        nreads = float(nreads_dict[sample])
        df_pivot[sample] = df_pivot[sample]/(nreads/1000000)

    df_final = pd.DataFrame()
    for sample in df_pivot.columns:
        if "H3-" not in sample:
            condition, mark, replicate = sample.split("-")
            control = "%s-H3-7" % condition
            control1 = "Saline-H3-7"
            control2 = "Dex-H3-7"

            # divide by average of H3 IPs
            df_final[sample] = np.divide(
                df_pivot[sample], ((df_pivot[control1] + df_pivot[control2])/2))

    plotHeatmaps = R('''
    function(df_final){
    dist_matrix_samples = as.dist(1-cor(df_final, method="spearman"))
    dist_matrix_peaks = as.dist(1-cor(t(df_final), method="spearman"))

    clust_samples = hclust(dist_matrix_samples, method="average")
    clust_peaks = hclust(dist_matrix_peaks, method="ward.D2")

    library(cowplot)
    library(ggplot2)
    library(ggdendro)
    library(reshape2)
    library(gridExtra)
    library(grid)
    library(lattice)

    k_columns = 3
    k_means = cutree(clust_samples, k=k_columns)

    ddata <- dendro_data(clust_samples, type = "rectangle")

    k_means = as.data.frame(k_means)
    k_means$label = rownames(k_means)

    labels = merge(label(ddata), k_means, by="label")

    # makes the colours match up with the enrichment plots
    #labels$sample_id = c(rep(2, 3), rep(6, 3), rep(4, 3),
    #                     rep(1, 3), rep(5, 3), rep(3, 3))

    write.table(df_final, "./test1.tsv", sep="\t")
    write.table(segment(ddata), "./test2.tsv", sep="\t")

    p_dend_sample <- ggplot(segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_rect(data = labels, aes(
        xmin=x-0.5,  xmax=x+0.5, ymin=y-0.3, ymax=y-0.1,
        fill=factor(sample_id), colour=factor(sample_id))) +
    scale_fill_brewer(palette="Paired") +
    scale_colour_brewer(palette="Paired") +
    theme_classic() +
    scale_x_continuous(trans="reverse") +
    scale_y_continuous(breaks=seq(0, 0.8, 0.2)) +
    xlab("") + ylab("") +
    theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size=15),
    legend.position="none")

    
    # copied from cowplot (function is now deprecated)
    ggplot_to_gtable <- function(plot)
    {
      if (methods::is(plot, "ggplot")){
        # ggplotGrob must open a device and when a multiple page capable device
        # (e.g. PDF) is open this will save a blank page
        # in order to avoid saving this blank page to the final target device
        # a NULL device is opened and closed here to *absorb* the blank plot

        # commenting this out to see if it was the cause of
        grDevices::pdf(NULL)
        plot <- ggplot2::ggplotGrob(plot)
        grDevices::dev.off()
        plot
      }
      else if (methods::is(plot, "gtable")){
        plot
      }
      else{
        stop('Argument needs to be of class "ggplot" or "gtable"' )
      }
    }

    switch_axis_position <- function(
        plot, axis = c('y', 'x', 'xy'), keep = c('none', 'x', 'y', 'xy', 'yx'))
    {
      keep.x <- switch(keep[1],
                       x = TRUE,
                       y = FALSE,
                       xy = TRUE,
                       yx = TRUE,
                       FALSE
      )
      keep.y <- switch(keep[1],
                       x = FALSE,
                       y = TRUE,
                       xy = TRUE,
                       yx = TRUE,
                       FALSE
      )

      # extract gtable
      gt <- ggplot_to_gtable(plot)

      # extract theme
      theme <- plot_theme(plot)

      result <- switch(axis[1],
                       x = switch_xaxis_position(gt, theme, keep.x),
                       y = switch_yaxis_position(gt, theme, keep.y),
                       switch_xaxis_position(
                                switch_yaxis_position(gt, theme, keep.y),
                                theme, keep.x)
      )
      result
    }
    
    p_dend_sample = ggdraw(switch_axis_position(p_dend_sample))

    k_rows = 5

    k_means = cutree(clust_peaks, k=k_rows)

    ddata <- dendro_data(clust_peaks, type = "rectangle")

    k_means = as.data.frame(k_means)
    k_means$label = rownames(k_means)

    labels = merge(label(ddata), k_means, by="label")

    p_dend_peak <- ggplot(segment(ddata)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_rect(data = labels, aes(
        xmin=x+0.01,  xmax=x+0.01, ymin=y-1, ymax=y,
        fill=paste0("cluster " ,k_means),
        colour=paste0("cluster " ,k_means))) +
    scale_fill_hue(c=60, l=70) +
    scale_colour_hue(c=60, l=70) +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    xlab("") + ylab("") +
    theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position="top")

    p_dend_peak = p_dend_peak + theme(legend.position="none")

    rename_sample_ids <-function(sample_ids){
        sample_ids <- gsub("Saline", "Veh", sample_ids)
        sample_ids <- gsub("\\\.1", " (1)", sample_ids)
        sample_ids <- gsub("\\\.3", " (1)", sample_ids)
        sample_ids <- gsub("\\\.7", " (2)", sample_ids)
        sample_ids <- gsub("\\\.9", " (3)", sample_ids)
        sample_ids <- gsub("\\\.", " - ", sample_ids)
        return(sample_ids)}

    get_legend<-function(myggplot){
        tmp <- ggplot_gtable(ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
        }


    df_final$peak = rownames(df_final)

    df_melt = melt(df_final, id="peak")


    df_melt$peak = factor(df_melt$peak, levels=rownames(df_final)[clust_peaks$order])
    df_melt$variable = factor(
        rename_sample_ids(df_melt$variable),
        levels=rename_sample_ids(rev(colnames(df_final)[clust_samples$order])))

    p = ggplot(df_melt, aes(variable, peak, fill=log(value, 2))) +
    geom_tile() +
    xlab("") + ylab("") +
    theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.direction="horizontal",
    legend.text=element_text(size=20),
    legend.title=element_text(size=20),
    legend.key.size=unit(1, "cm"),
    axis.text.x=element_text(angle=90, size=20, hjust=1, vjust=0.5),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    legend.title.align=0.5)  +
    scale_fill_gradient2(
        high="slateblue4", mid="grey90", low="grey90",
        midpoint=0, breaks=seq(0, 4, 2), limits=c(0,4),
        name="Enrichment (log2)",
        guide = guide_colorbar(
            title.position="top", title.hjust = 0.5, label.position="bottom"))

    legend = get_legend(p)

    p_blank =  textGrob("text")

    p_dend_sample_final = p_dend_sample +
    theme(plot.margin = unit(c(3,-0.51,-1.7,-1.33), "cm"))

    p_dend_peak_final = p_dend_peak +
    theme(plot.margin = unit(c(0,0,6.3,0), "cm"))

    p_no_legend = p + theme(
    legend.position="none",
    plot.margin = unit(c(0.95,2.25,0.5,-0.9), "cm"))

    png("%(outfile)s", width=36,height=36, units="cm", res=400)
    grid.arrange(legend,
                 p_dend_sample_final,
                 p_dend_peak_final,
                 p_no_legend,
                 ncol=2, widths = c(1, 3),  heights = c(1, 4))
    dev.off()

    } ''' % locals())

    plotHeatmaps(df_final)
