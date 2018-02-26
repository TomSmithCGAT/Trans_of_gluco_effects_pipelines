'''
PipelineProject34sRNA.py - Utility functions for Project 34 sRNA-Seq analysis
=============================================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Code
----

'''

import collections
import os
import copy
import pandas as pd
import sqlite3
import numpy as np
from rpy2.robjects import r as R
import rpy2.robjects as ro
import pandas.rpy.common as com
from rpy2.robjects.packages import importr
import MySQLdb
import pysam

from CGATPipelines.Pipeline import cluster_runnable
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq
import CGAT.Counts as Counts
import CGAT.Expression as Expression
import CGAT.FastaIterator as FastaIterator
import CGATPipelines.PipelineTracks as PipelineTracks


@cluster_runnable
def collapseFastq(infile, outfile, collapse_log):

    seq_dict_count = collections.defaultdict(int)
    seq_dict_score = collections.defaultdict(int)
    seq_dict_entry = collections.defaultdict()
    seq2identifiers = collections.defaultdict(list)

    for entry in Fastq.iterate(IOTools.openFile(infile, "r")):
        entry.format = "sanger"

        seq_dict_count[entry.seq] += 1
        seq2identifiers[entry.seq].append(entry.identifier)

        sum_qual = sum(entry.toPhred())

        if sum_qual > seq_dict_score[entry.seq]:
            seq_dict_score[entry.seq] = sum_qual
            seq_dict_entry[entry.seq] = entry

    with IOTools.openFile(outfile, "w") as outf:
        for x in seq_dict_count:
            entry = seq_dict_entry[x]

            children = seq2identifiers[entry.seq]

            entry.identifier = (entry.identifier.split(" ")[0] +
                                "_" + str(seq_dict_count[x]))
            outf.write("%s\n" % entry)


def connectToUCSC(database):
    dbhandle = MySQLdb.Connect(host="genome-mysql.cse.ucsc.edu",
                               user="genome")

    cc = dbhandle.cursor()
    cc.execute("USE %s " % database)

    return dbhandle


def getRepeatsFromUCSC(dbhandle, repclasses, outfile,
                       ensembl_remove_contigs=False):
    '''select repeats from UCSC and write to *outfile* in gff format.

    If *repclasses* is None or an empty list, all repeats will be collected.
    '''

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is
    # queried for tables that end in rmsk.
    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    # now collect repeats
    tmpfile = IOTools.openFile(outfile, "w")

    for table in tables:

        cc = dbhandle.cursor()
        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd,
        '.', strand, '.',
        CONCAT('gene_id \\"', repName, '\\"; transcript_id \\"',
        repName, '\\"; class \\"', repClass, '\\"; family \\"',
        repFamily, '\\"' )
        FROM %(table)s"""

        if repclasses:
            repclasses_str = ",".join(
                ["'" + x.strip() + "'" for x in repclasses])

            sql += ''' WHERE repClass in (%(repclasses_str)s) ''' % locals()

        sql = sql % locals()

        E.debug("executing sql statement: %s" % sql)
        cc.execute(sql)

        for data in cc.fetchall():
            tmpfile.write("\t".join(map(str, data)) + "\n")
    tmpfile.close()


def importRNAAnnotationFromUCSC(outfile, repclasses=None, database="rn5"):
    '''import repetetive RNA from a UCSC formatted file.

    The repetetive RNA are taken from the repeat-masker track.

    The results are stored as a :term:`gff` formatted file.
    '''

    dbhandle = connectToUCSC(database)
    getRepeatsFromUCSC(dbhandle, repclasses, outfile)


def iterativeBowtieMapping(infile, fastas, outfile,
                           mismatches, seed_length, mapping_order,
                           multimapping_max):
    ''' perform iterative mapping using bowtie
    alignment is performed in the order specified in mapping_order param
    Initially only uniquely aligning reads are included'''

    def getCounts(read):
        ''' return counts from pysam read '''
        return int(read.query_name.strip().split("_")[1])

    def getID(read):
        ''' return ID from psyam read '''
        return read.query_name.strip().split("_")[0]

    def getIDCounts(read):
        ''' return ID and Counts from psyam read '''
        x, y = read.query_name.strip().split("_")
        return (x, int(y))

    def getFastaID(entry):
        ''' return ID from fasta entry '''
        return entry.title.split(" ")[0]

    def updateCountsAndCoverage(count, possible_alignments,
                                ensembl2coverage, ensembl2identifier):
        '''take a read and the possible alignments and update the counts and
        coverage dictionary'''

        coverage_array = []

        for alignment in possible_alignments:

            if alignment.is_reverse:
                strand = "reverse"
            else:
                strand = "forward"

            coverage = ensembl2coverage[
                reference2ensembl[alignment.reference_id]][strand][
                alignment.aend - alignment.alen:alignment.aend]

            coverage_array.append(np.mean([x + 1 for x in coverage]))

        coverage_array = [np.divide(float(x), sum(coverage_array))
                          for x in coverage_array]

        random_selection = np.random.choice(len(coverage_array), count,
                                            p=coverage_array)

        for rand in random_selection:

            selected = possible_alignments[rand]

            if selected.is_reverse:
                strand = "reverse"
            else:
                strand = "forward"

            gene_id = reference2ensembl[selected.reference_id]

            for position in range(selected.aend - selected.alen,
                                  selected.aend, 1):
                ensembl2coverage[gene_id][strand][position] += 1

            ensembl2coverage[gene_id]['hits'] += 1
            ensembl2identifier[gene_id].update((read.query_name,))

        return ensembl2coverage, ensembl2identifier

    ensembl2annotation = collections.defaultdict(str)

    # keep track of which fastq identifiers align to which ensembl ids
    ensembl2identifier = collections.defaultdict(set)

    # create a dictionary object to hold the counts per gene and coverage
    ensembl2coverage = collections.defaultdict(
        lambda: collections.defaultdict(np.array))

    # dictionary to hold index name for each annotation
    annotation2bowtieIndex = {}

    for annotation in mapping_order.strip().split(","):
        for fasta in fastas:
            if annotation in fasta:
                annotation2bowtieIndex[annotation] = P.snip(fasta, ".fa")

                for entry in FastaIterator.FastaIterator(
                        IOTools.openFile(
                            annotation2bowtieIndex[annotation] + ".fa", "r")):

                    entryID = getFastaID(entry)
                    seq_len = len(entry.sequence)
                    ensembl2annotation[entryID] = annotation

                    ensembl2coverage[entryID]['forward'] = np.zeros(
                        (seq_len,), dtype=np.int)

                    ensembl2coverage[entryID]['reverse'] = np.zeros(
                        (seq_len,), dtype=np.int)

                    ensembl2coverage[entryID]['hits'] = 0

    # dictionary to hold the number of reads assigned to each annotation per
    # multimapping number
    countsPerAnnotation = collections.defaultdict(
        lambda: collections.defaultdict(int))

    unmapped_fastq_outfile = P.snip(outfile, "_counts.tsv") + "_unmapped.fastq"

    samfile = "%s.sam" % P.snip(outfile, "_counts.tsv")

    for multimap in range(1, multimapping_max + 1):

        for annotation in mapping_order.strip().split(","):

            bowtie_index = annotation2bowtieIndex[annotation]

            statement = '''
            bowtie --sam --strata --best
            -n %(mismatches)s
            -l %(seed_length)s
            -m %(multimap)s
            -a %(bowtie_index)s %(infile)s > %(samfile)s
            '''

            P.run()

            insam = pysam.Samfile(samfile, "r")
            reference2ensembl = {x: y for x, y in enumerate(insam.references)}

            inreads = insam.fetch()

            current_id = ""
            possible_alignments = 0
            with IOTools.openFile(unmapped_fastq_outfile, "w") as outf:

                for read in inreads:

                    if read.is_unmapped:
                        continue

                    else:
                        read_id = getID(read)

                        if read_id == current_id:
                            possible_alignments.append(read)

                        else:
                            # check whether current_id has been set yet,
                            # if so, update coverage data
                            if current_id:
                                count = getCounts(possible_alignments[0])

                                ensembl2coverage, ensembl2identifier = (
                                    updateCountsAndCoverage(
                                        count, possible_alignments,
                                        ensembl2coverage, ensembl2identifier))

                                countsPerAnnotation[annotation][multimap] += count

                            # reset the current_id and possible_alignments
                            current_id = read_id
                            possible_alignments = [read]

                # pick up the last set of alignments
                # need to check any alignments have actually been found
                if possible_alignments:

                    count = getCounts(possible_alignments[0])
                    ensembl2coverage, ensembl2identifier = updateCountsAndCoverage(
                        count, possible_alignments, ensembl2coverage,
                        ensembl2identifier)

            statement = """
            samtools view -f4 %(samfile)s|
            awk '{printf("@%%s\\n%%s\\n+\\n%%s\\n",$1,$10,$11)}' >
            %(unmapped_fastq_outfile)s"""

            P.run()

            # reset infile to unmapped fastq file so next mapping used unmapped
            infile = unmapped_fastq_outfile

    # write out files
    with IOTools.openFile(outfile, "w") as outf:
        for x in ensembl2coverage:
            outf.write("\t".join((ensembl2annotation[x], x,
                                  str(ensembl2coverage[x]['hits']))) + "\n")

    # stop numpy adding line breaks
    np.set_printoptions(linewidth=100000)

    forward_outfile = "%s_cov_forward.tsv" % P.snip(outfile, "_counts.tsv")
    with IOTools.openFile(forward_outfile, "w") as outf:
        for x in ensembl2coverage:
            outf.write("%s\t%s\t%s\n" % (
                (ensembl2annotation[x], x,
                 ",".join(map(str, ensembl2coverage[x]['forward'])))))

    reverse_outfile = "%s_cov_reverse.tsv" % P.snip(outfile, "_counts.tsv")
    with IOTools.openFile(reverse_outfile, "w") as outf:
        for x in ensembl2coverage:
            outf.write("%s\t%s\t%s\n" % (
                (ensembl2annotation[x], x,
                 ",".join(map(str, ensembl2coverage[x]['reverse'])))))

    ensembl2fastq_outfile = "%s_ensembl2fastq.tsv" % P.snip(
        outfile, "_counts.tsv")
    with IOTools.openFile(ensembl2fastq_outfile, "w") as outf:
        for x in ensembl2identifier:
            outf.write("%s\t%s\n" % (
                x, ",".join(list(ensembl2identifier[x]))))

    summary_outfile = "%s_summary.tsv" % P.snip(outfile, "_counts.tsv")
    with IOTools.openFile(summary_outfile, "w") as outf:
        for annotation in countsPerAnnotation:
            for multimap in countsPerAnnotation[annotation]:
                outf.write("%s\t%i\t%i\n" % (
                    annotation, multimap,
                    countsPerAnnotation[annotation][multimap]))


def collapsetRNAs(infile, outfile, length):

    def getID(read):
        ''' return ID from read '''
        return read.title.strip().split(" ")[0]

    def get5PrimeSequence(read, length=18):
        ''' return ID from read '''
        return read.sequence.strip()[0:length]

    fasta = FastaIterator.iterate(IOTools.openFile(infile, "r"))

    FivePrimeSequence2Loci = collections.defaultdict(list)

    for entry in fasta:
        FivePrimeSequence2Loci[get5PrimeSequence(entry, length)].append(
            getID(entry))

    with IOTools.openFile(outfile, "w") as outf:
        for sequence, tRNAs in FivePrimeSequence2Loci.iteritems():
            if len(sequence) == length:
                outf.write("%s\t%s\n" % (sequence, ",".join((tRNAs))))


def countReadsPertRF(infile, outfile, collapsed_tRNAs, length):

    # read the collapsed_tRNAs in as a dictionary
    FivePrimeSequence2Loci = {}
    with IOTools.openFile(collapsed_tRNAs, "r") as inf:
        for line in inf:
            FivePrimeSequence2Loci[line.split("\t")[0]] = [
                x for x in line.strip().split("\t")[1].split(",")]


    # get the counts per tRNA
    # use maximum coverage over 5' region (default first 30 nt)
    tRF_counts = {}

    with IOTools.openFile(infile, "r") as inf:
        for line in inf:
            annot_class, name, values = line.strip().split("\t")

            if annot_class == "tRNA":
                coverage_values = map(int, [x for x in values.split(",")])

                if len(coverage_values) >= length:
                    tRF_counts[name] = max(coverage_values[0:30])

    # create a dataframe and write out
    rownames = []
    total_counts = []
    for sequence, tRNAs in FivePrimeSequence2Loci.iteritems():
        rownames.append(sequence)
        total_counts.append(sum([tRF_counts[x] for x in tRNAs]))

    collapsed_counts = pd.DataFrame(
        {"counts": total_counts, "gene": rownames})
    collapsed_counts.sort(inplace=True)

    collapsed_counts.to_csv(outfile, sep="\t", index=False)


def plottRNAProfile(infile, outfiles, length):

    sum_outfile, norm_outfile = outfiles

    # tRNAs must be at least 40 bp long to be included
    # only the first 40bp will be considered
    seq_max = 40
    seq_range = range(0, seq_max)

    rownames = []
    tRNA_coverage = {x: [] for x in seq_range}

    with IOTools.openFile(infile, "r") as inf:
        for line in inf:
            annot_class, name, values = line.strip().split("\t")

            if annot_class == "tRNA":
                coverage_values = map(int, [x for x in values.split(",")])

                if len(coverage_values) >= seq_max:
                    rownames.append(line.split("\t")[1])
                    for x in seq_range:
                        tRNA_coverage[x].append(coverage_values[x])

    df = pd.DataFrame(tRNA_coverage, index=rownames)

    # only include tRFs with at least 10 counts
    df = df.ix[df.apply(axis=1, func=np.mean) > 10, :]

    df_col_avr = pd.DataFrame(
        {"nt": seq_range, "mean": df.apply(axis=0, func=sum)})

    # Need to normalise to total number of counts
    df_col_avr['mean'] = [float(x)/sum(df_col_avr['mean']) for
                          x in df_col_avr['mean']]

    df_col_avr.to_csv(sum_outfile, sep="\t", index=False)

    def normalise(row):
        return [float(x)/sum(row) for x in row]

    df_norm = df.apply(axis=1, func=normalise)

    df_norm_col_avr = pd.DataFrame(
        {"nt": seq_range,
         "mean": df_norm.ix[
             df.apply(axis=1, func=np.mean) > 100, :].apply(axis=0, func=sum)})

    num_genes = df.shape[0]

    # Need to normalise to total number of genes
    df_norm_col_avr['mean'] = [float(x)/num_genes for
                               x in df_norm_col_avr['mean']]

    df_norm_col_avr.to_csv(norm_outfile, sep="\t", index=False)

    plotCoverage = R('''
    function(df, plot_outfile){
    suppressMessages(library(ggplot2))

    l_txt = element_text(size=20)

    p = ggplot(df, aes(x=nt, y=mean)) +
    geom_line(size=2) +
    geom_point(size=4) +
    theme(axis.text.x = l_txt,
          axis.text.y = l_txt,
          axis.title.x = l_txt,
          axis.title.y = l_txt) +
    xlab("Base pair (from 5')") +
    ylim(0, max(df$mean)) +
    ylab("Coverage")
    ggsave(plot_outfile, width=10, height=10)
    }
    ''')

    sum_plot_outfile = sum_outfile.replace(".tsv", ".png")
    norm_plot_outfile = norm_outfile.replace(".tsv", ".png")
    r_sum_plot_outfile = ro.StrVector(sum_plot_outfile)
    r_norm_plot_outfile = ro.StrVector(norm_plot_outfile)

    r_df_col_avr = com.convert_to_r_dataframe(df_col_avr)
    r_df_norm_col_avr = com.convert_to_r_dataframe(df_norm_col_avr)

    plotCoverage(r_df_col_avr, sum_plot_outfile)
    plotCoverage(r_df_norm_col_avr, norm_plot_outfile)


def runDESeq2(counts_infile, design_infile, outfile):

    plot_base = P.snip(outfile, ".deseq.tsv")
    logfile = plot_base + ".log"
    outfile_summary = plot_base + ".summary.tsv"
    dispersion_plot_filename = plot_base + ".dispersion.png"
    transformations_plot_filename = plot_base + ".transformations.png"
    heatmap_plot_filename = plot_base + ".heatmap.png"
    MA_plot_filename = plot_base + "_MA_plot.png"
    unshrunken_MA_plot_filename = plot_base + "_unshrunken_MA_plot.png"

    # create Counts object
    counts = Counts.Counts(pd.io.parsers.read_csv(
        counts_infile, sep="\t", index_col=0, comment="#"))

    # create Design object
    design = Expression.ExperimentalDesign(pd.read_csv(
        design_infile, sep="\t", index_col=0, comment="#"))

    # validate design against counts
    design.validate(counts)

    # restrict counts to samples in design table
    counts.restrict(design)

    # plot transformations
    # counts.plotTransformations(plot_filename=transformations_plot_filename)

    colData = pd.DataFrame(design.table.loc[:, "treatment"])

    # convert counts object to an r object
    r_counts = com.convert_to_r_dataframe(counts.table)
    r_design = com.convert_to_r_dataframe(colData)

    runDESeq2 = R('''
    function(counts, design){

      suppressMessages(library(DESeq2))
      suppressMessages(library(vsn))
      suppressMessages(library(gplots))
      suppressMessages(library(RColorBrewer))

      dds <- DESeqDataSetFromMatrix(countData = counts, colData = design,
                                    design = ~ treatment)
      dds <- DESeq(dds)
      res <- results(dds, alpha = 0.05)
      resOrdered <- res[order(res$pvalue),]

      png("%(MA_plot_filename)s")
      plotMA(res, main="DESeq2", ylim=c(-3,3))
      dev.off()

      png("%(unshrunken_MA_plot_filename)s")
      resMLE <- results(dds, addMLE=TRUE)
      resMLE$log2FoldChange = resMLE$lfcMLE
      plotMA(resMLE, main="DESeq2", ylim=c(-3,3))
      dev.off()

      png("%(dispersion_plot_filename)s")
      plotDispEsts(dds)
      dev.off()

      rld <- rlogTransformation(dds, blind=TRUE)
      vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

      hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
      distsRL <- dist(t(assay(rld)))
      mat <- as.matrix(distsRL)
      rownames(mat) <- colnames(mat) <- dds$treatment
      hc <- hclust(distsRL)
      png("%(heatmap_plot_filename)s")
      heatmap.2(mat, Rowv=as.dendrogram(hc),
      symm=TRUE, trace="none",
      col = rev(hmcol), margin=c(13, 13))
      dev.off()

      res$gene_id = row.names(res)
      write.table(res, "%(outfile)s", sep="\t", quote=FALSE, row.names=FALSE)

      sink("%(outfile_summary)s")
      summary(res)
      sink()

    }''' % locals())

    runDESeq2(r_counts, r_design)


@cluster_runnable
def plotPCAsAndDendo(counts_inf, design_inf, outfile):
    ''' make basic PCA and dendogram plots from counts object '''

    plot_base = P.snip(outfile, ".pca_PC1_PC2.png")
    logfile = plot_base + ".log"
    dendogram_plot_filename = (plot_base + ".dendogram.png")
    coloured_dendogram_plot_filename = plot_base + ".coloured.dendogram.png"
    pca_variance_plot_filename = (plot_base + ".pca_variance_explained.png")
    pca_plot_filename = (plot_base + ".pca_PC1_PC2.png")
    pca_plot2_filename = (plot_base + ".pca_PC3_PC4.png")

    # create Counts object
    counts = Counts.Counts(pd.io.parsers.read_csv(
        counts_inf, sep="\t", index_col=0, comment="#"))

    # create Design object
    design = Expression.ExperimentalDesign(pd.read_csv(
        design_inf, sep="\t", index_col=0, comment="#"))

    # validate design against counts
    design.validate(counts)

    # restrict counts to samples in design table
    counts.restrict(design)

    design.revalidate(counts)
    counts.normalise(method="total-column")
    counts.log(base=2, pseudocount=0.1)

    outf = IOTools.openFile(logfile, "w")
    outf.write("pre-filter read data: %i observations for %i samples\n"
               % counts.table.shape)

    # remove gene models with low counts (less than 2/1000000)
    counts.removeObservationsFreq(min_counts_per_row=1)

    outf.write("post-filter read data: %i observations for %i samples\n"
               % counts.table.shape)

    # make PC plots
    counts.plotPCA(design, variance_plot_filename=pca_variance_plot_filename,
                   pca_plot_filename=pca_plot_filename)

    counts.plotPCA(design, variance_plot_filename=pca_variance_plot_filename,
                   pca_plot_filename=pca_plot2_filename,
                   x_axis="PC3", y_axis="PC4")
    # cluster
    counts.plotDendogram(plot_filename=dendogram_plot_filename,
                         distance_method="euclidean",
                         clustering_method="ward.D2")

    coloured_dendogram = R('''
      function(df){

        suppressMessages(library(dendextend))

        hc = hclust(dist(t(df), method = "euclidean"), method = "ward")
        dend <- as.dendrogram(hc)

        groupCodes <- c(rep("Saline",4), rep("Dex",4))
        colorCodes <- c(Saline="darkorange4", Dex="chartreuse4")
        labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

        png("%(coloured_dendogram_plot_filename)s")
        par(cex=1.5, mar = c(8,4,1,1))
        plot(dend)
        dev.off()

    }''' % locals())

    r_df = com.convert_to_r_dataframe(counts.table)

    coloured_dendogram(r_df)

    outf.close()


@cluster_runnable
def plotHeatmap(counts_inf, design_inf, outfile, short_names):
    ''' plot Heatmap
    Genes restricted to those with at least 5 counts
    '''

    counts = Counts.Counts(counts_inf)
    counts.removeObservationsFreq(5)
    retain_genes = counts.table.index

    design = Expression.ExperimentalDesign(design_inf)

    counts = Counts.Counts(counts_inf)
    counts.restrict(design)
    counts.transform(method="rlog", design=design, blind=False)
    counts.table = counts.table.ix[retain_genes]

    df = counts.table
    df.columns = [x.replace("Saline", "Veh") for x in df.columns]

    max_df = df.sum(1)
    max_df.sort()
    df = df.ix[max_df.index]

    if short_names:
        #remove generation portion from sample name
        df.columns = ["-".join(x.split("-")[1:]) for x in df.columns]

    plotHeatmap = R('''function(df){

    write.table(df, paste0("%(outfile)s", "_r.tsv"), sep="\t")

    library("Biobase")
    library("RColorBrewer")
    library("pvclust")
    library("ggplot2")
    library("gplots")

    result = pvclust(as.matrix(df), method.dist="correlation",
                     method.hclust="ward.D2", nboot=1000)

    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(99)

    treatment_number = sapply(
        colnames(df), FUN=function(x) ifelse(grepl("Veh",x), 1, 2))
    colcols = c("#9e66ab", "#f9a65a")[treatment_number]

    png("%(outfile)s", width = 7, height = 8, units = 'in', res = 590)

    colnames(df) <- gsub("\\\.", " ", colnames(df))

    heatmap.2(as.matrix(df),
          col = hmcol,
          key = FALSE,
          ColSideColors = colcols,
          scale="none",
          symbreaks=FALSE,
          trace="none",
          Colv = as.dendrogram(result$hclust),
          margin=c(16, 1),
          dendrogram="column",
          cexCol=3,
          labRow = "",
          Rowv=FALSE,
          hclustfun = function(x) hclust(x, method = 'ward.D2'))

    legend("topleft", legend=c("Veh", "Dex"),
        fill= c("#9e66ab", "#f9a65a"),
        border=FALSE, bty="n", y.intersp = 0.7, cex=1)

    dev.off()

    }''' % locals())

    r_df = com.convert_to_r_dataframe(df)

    plotHeatmap(r_df)


@cluster_runnable
def plotCorrelationHeatmap(counts_inf, design_inf, outfile, short_names):
    ''' plot Heatmap
    Genes restricted to those with at rlog >5
    '''

    counts = Counts.Counts(counts_inf)
    retain_genes = counts.table.index

    design = Expression.ExperimentalDesign(design_inf)

    counts = Counts.Counts(counts_inf)
    counts.restrict(design)
    counts.transform(method="rlog", design=design, blind=False)
    counts.table = counts.table.ix[retain_genes]

    df = counts.table
    df.columns = [x.replace("Saline", "Veh") for x in df.columns]

    df = df.ix[df.min(axis=1) > 5]

    if short_names:
        #remove generation portion from sample name
        df.columns = ["-".join(x.split("-")[1:]) for x in df.columns]

    if "miRNA" in counts_inf:
        pad = "0.003"
    elif "piR" in counts_inf:
        pad = "0.01"
    elif "tRNA" in counts_inf or "tRF" in counts_inf:
        pad = "0.002"

    plotHeatmap = R('''function(df){

    library(dendextend)
    library(ggplot2)
    library(reshape2)
    library(ggdendro)
    library(gplots)
    library(gridExtra)
    library(grid)
    library(lattice)
    library(pvclust)

    spearman <- function(x, ...) {
    x <- as.matrix(x)
    res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
    res <- as.dist(res)
    attr(res, "method") <- "spearman"
    return(res)
    }

    pv_result = pvclust(as.matrix(df),
    method.dist=spearman, method.hclust="ward.D2", nboot=1000)

    ddata <- dendro_data(as.dendrogram(pv_result), type = "rectangle")

    labels_df <- ddata$labels
    labels_df['id'] <- lapply(labels_df['label'],
                              FUN=function(x) ifelse(grepl("Veh", x), 0, 1))

    m=30
    m_txt = element_text(size=m)

    blank_theme = theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y = m_txt,
    axis.title.x = m_txt,
    axis.title.y = m_txt,
    legend.text = m_txt,
    legend.title = m_txt,
    legend.key.size = unit(2, "cm"),
    legend.position = "right"
    )

    p_dend <- ggplot(segment(ddata)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_rect(data = labels_df, aes(
            xmin=x-0.46,  xmax=x+0.46,
            ymin=y-%(pad)s, ymax=y,
            fill=factor(id), colour=factor(id))) +
    scale_fill_manual(values=c("#DA9870", "#AB98C8")) +
    scale_colour_manual(values=c("#DA9870", "#AB98C8")) +
    blank_theme +
    theme(
    axis.text.y = element_text(size=20),
    axis.text.x=element_blank(),
    legend.position="None",
    plot.margin = unit(c(0, 7.5, -1.7, 2.4), "cm")) +
    #plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    xlab("") + ylab("")

    ggsave("%(outfile)s_dend.png")

    cor_df = as.data.frame(cor(df, method="spearman"))
    cor_df$v1 <- factor(rownames(cor_df), levels=ddata$labels$label)
    cor_df = melt(cor_df, id_vars="v1")
    cor_df$variable <- factor(cor_df$variable, levels=ddata$labels$label)

    p = ggplot(cor_df, aes(v1, variable, fill=value)) +
    geom_tile() +
    blank_theme +
    scale_fill_gradient(
        low="white", high="skyblue4",
        limits=c(0.85, 1), breaks = seq(0.85, 1, 0.05),
        name="Spearman's\nrho\n") +
    theme(
    axis.text.x=element_text(size=m, angle=90, hjust=1, vjust=0.5),
    aspect.ratio=1,
    plot.margin = unit(c(0,0,0,0), "cm")) +
    xlab("") + ylab("")

    png("%(outfile)s", width = 42, height = 42, units = 'cm', res = 400)
    grid.arrange(p_dend, p, ncol=1, heights = c(1, 3))
    dev.off()

    }''' % locals())

    r_df = com.convert_to_r_dataframe(df)

    plotHeatmap(r_df)

@cluster_runnable
def spikeInDESeq2(counts, design, outfile_prefix):
    ''' run DESeq2 +/- spike ins and plot MA and power '''

    spike_counts = Counts.Counts(
        pd.read_table(counts, sep="\t", index_col=0, comment="#"))

    design = Expression.ExperimentalDesign(design)

    spike_counts.restrict(design)

    # run initial analysis without spike ins
    spike_counts.table = spike_counts.table.ix[
        ["spike" not in x for x in spike_counts.table.index]]
    experiment = Expression.DEExperiment_DESeq2()

    design.factors['group'] = design.table['group']

    # run DESeq2
    results = experiment.run(spike_counts,
                             design,
                             model=None,
                             contrasts=None,
                             outfile_prefix=outfile_prefix,
                             fdr=0.1,
                             ref_group="F1-Saline")
    # extract results
    results.getResults(fdr=0.1)
    results.summariseDEResults()

    outfile = outfile_prefix + "_MA_no_spike.png"

    plot_MA = R('''
    function(df){

    library(ggplot2)

    m=40
    m_txt = element_text(size=m)

    p = ggplot(df, aes(log(control_mean,2), l2fold)) +
    geom_point(size=1.5, colour="black") +
    ylab("Fold change (Log 2)") +
    xlab("Mean expression (Log 2)") +
    scale_colour_discrete(name="Significant") +
    theme_bw() +
    scale_y_continuous(breaks=seq(-8,8,2), limits=c(-8,8)) +
    scale_x_continuous(breaks=seq(-5,20,5), limits=c(-5,20)) +
    theme(
    aspect.ratio=1,
    axis.text.x=m_txt,
    axis.text.y=m_txt,
    axis.title.x=m_txt,
    axis.title.y=m_txt,
    legend.text=m_txt,
    legend.title=m_txt)
    guides(colour = guide_legend(override.aes = list(size=3)))
    ggsave("%(outfile)s", width=10, height=10)
    }''' % locals())

    plot_MA(com.convert_to_r_dataframe(results.table))

    # repeat with spike-ins
    spike_counts = Counts.Counts(
        pd.read_table(counts, sep="\t", index_col=0, comment="#"))

    spike_counts.restrict(design)
    experiment = Expression.DEExperiment_DESeq2()

    design.factors['group'] = design.table['group']

    # run DESeq2
    results = experiment.run(spike_counts,
                             design,
                             model=None,
                             contrasts=None,
                             outfile_prefix=outfile_prefix,
                             fdr=0.1,
                             ref_group="F1-Saline")
    # extract results
    results.getResults(fdr=0.1)
    results.summariseDEResults()

    df_spike = copy.deepcopy(results.table)
    df_spike.set_index("test_id", inplace=True)
    df_spike = df_spike.ix[[x for x in df_spike.index if "spike" in x]]
    df_spike['spike_fold_change'] = [x.split("_")[2] for x in df_spike.index]
    df_spike['initial'] = [x.split("_")[1] for x in df_spike.index]

    outfile = outfile_prefix + "_spike_in_power.png"

    plot_power = R('''
    function(df){
    library(ggplot2)
    library(grid)
    
    s=15
    m=20
    s_txt = element_text(size=s)
    m_txt = element_text(size=m)

    p = ggplot(df, aes(as.numeric(as.character(initial)),
                       significant,
                       colour=as.factor(spike_fold_change))) +
    stat_summary(fun.y = "mean", geom="line", size=2) +
    ylab("Power") +
    xlab("Mean expression (Log 2)") +
    scale_colour_discrete(name="Fold-change (Log 2)") +
    theme_bw() +
    scale_x_continuous(breaks=seq(1,9,2), limits=c(1,9)) +
    theme(
    aspect.ratio=1,
    axis.text.x=m_txt,
    axis.text.y=m_txt,
    axis.title.x=m_txt,
    axis.title.y=m_txt,
    legend.text=m_txt,
    legend.title=m_txt,
    legend.key.size=unit(1, "cm")) +
    guides(colour = guide_legend(override.aes = list(size=2)))

    ggsave("%(outfile)s", width=10, height=10)
    }''' % locals())

    plot_power(com.convert_to_r_dataframe(df_spike))

    df = results.table
    df['spike'] = ["spike" in x for x in df['test_id']]
    df['significant'] = [1 if df['significant'][x] == 1
                         and df['spike'][x] == True
                         else 0 for x in df.index]

    outfile = outfile_prefix + "_MA_with_spike.png"

    plot_MA_2 = R('''
    function(df, outfile){

    library(ggplot2)
    library(grid)

    m=20
    s=15
    m_txt = element_text(size=m)
    s_txt = element_text(size=s)

    p = ggplot(df, aes(log(control_mean,2), l2fold,
                       colour=interaction(spike, as.factor(significant)))) +
    geom_point(size=2) +
    ylab("Fold change (Log 2)") +
    xlab("Mean expression (Log 2)") +
    scale_colour_discrete(
        name="", label=c("sRNA", "Spike-in - Not Sig.", "Spike-in - Sig.")) +
    theme_bw() +
    scale_y_continuous(breaks=seq(-8,8,2), limits=c(-8,8)) +
    scale_x_continuous(breaks=seq(-5,20,5), limits=c(-5,20)) +
    theme(
    aspect.ratio=1,
    axis.text.x=m_txt,
    axis.text.y=m_txt,
    axis.title.x=m_txt,
    axis.title.y=m_txt,
    legend.text=s_txt,
    legend.title=m_txt,
    legend.key.size=unit(1, "cm")) +
    guides(colour = guide_legend(override.aes = list(size=3)))

    ggsave("%(outfile)s", width=10, height=10)
    }''' % locals())

    plot_MA_2(com.convert_to_r_dataframe(df))
