'''
PipelineRRBS.py - Utility functions for mapping short reads
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Mapping reads is a common task in pipelines. Different pipelines
combine different sources of input (:term:`fastq` files, :term:`sra` files)
of different data (single end, paired end) with different mapping
algorithms (bowtie, tophat, stampy). This module provides utility
functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra` archive
might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

The module currently is able to deal with:

   * tophat mapping against genome
   * bowtie mapping against transcriptome, genome and junctions
   * bwa against genome
   * stampy against genome

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

Code
----

'''

import os
import sys
import shutil
import glob
import collections
import re
import gzip
import itertools
import CGAT.Pipeline as P
import logging as L
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Fastq as Fastq
import CGAT.IndexedFasta as IndexedFasta
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelienMapping
import pysam

SequenceInformation = collections.namedtuple("SequenceInformation",
                                             """paired_end
                                                 filename_first
                                                 filename_second
                                                 readlength_first
                                                 readlength_second
                                                 is_colour""")


def getReadLengthFromFastq(filename):
    '''return readlength from a fasta/fastq file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.

    '''

    with IOTools.openFile(filename) as infile:
        record = iterate(infile).next()
        readlength = len(record.seq)
        return readlength


def getReadLengthFromBamfile(filename):
    '''return readlength from a bam file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.
    '''

    samfile = pysam.Samfile(filename, "rb")
    record = samfile.fetch().next()
    readlength = record.rlen
    samfile.close()
    return readlength


def getSequencingInformation(track):
    '''glean sequencing information from *track*.'''

    colour = False
    if os.path.exists("%s.fastq.gz" % track):
        first_pair = "%s.fastq.gz" % track
        second_pair = None
    elif os.path.exists("%s.fastq.1.gz" % track):
        first_pair = "%s.fastq.1.gz" % track
        second_pair = "%s.fastq.2.gz" % track
    elif os.path.exists("%s.csfasta.gz" % track):
        first_pair = "%s.csfasta.gz" % track
        second_pair = None
        colour = True

    second_length = None
    if second_pair:
        if not os.path.exists(second_pair):
            raise IOError("could not find second pair %s for %s" %
                          (second_pair, first_pair))
        second_length = getReadLength(second_pair)

    return SequenceInformation._make((second_pair is not None,
                                      first_pair, second_pair,
                                      getReadLength(first_pair),
                                      second_length,
                                      colour))



###############################################################################
############################### Classes #######################################
###############################################################################
class Trimmer(PipelineMapping.Mapper):

    '''trims reads.

    preprocesses the input data, calls trimmer and post-process the output data.

    All in a single statement to be send to the cluster.

    Preprocessing function is defined in PipelineMapping.py
    '''

    datatype = "fastq"

    # set to True if you want to preserve colour space files.
    # By default, they are converted to fastq.
    preserve_colourspace = False

    # compress temporary fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    # remove non-unique matches in a post-processing step.
    # Many aligners offer this option in the mapping stage
    # If only unique matches are required, it is better to
    # configure the aligner as removing in post-processing
    # adds to processing time.
    remove_non_unique = False

    def __init__(self):
        pass

    def trimmer(self, infiles, outfile):
        '''build trimming statement on infiles.
        '''
        return ""

    def build(self, infiles, outfile):
        '''run mapper.'''

        cmd_preprocess, trimfiles = self.preprocess(infiles, outfile)
        cmd_trimmer = self.trimmer(trimfiles, outfile)
        cmd_postprocess = self.postprocess(infiles, outfile)
        cmd_clean = self.cleanup(outfile)

        assert cmd_preprocess.strip().endswith(";")
        assert cmd_trimmer.strip().endswith(";")
        if cmd_postprocess:
            assert cmd_postprocess.strip().endswith(";")
        if cmd_clean:
            assert cmd_clean.strip().endswith(";")

        statement = " checkpoint; ".join((cmd_preprocess,
                                          cmd_trimmer,
                                          cmd_postprocess,
                                          cmd_clean))

        return statement


class TrimGalore(Trimmer):

    '''run trim galore on reads.'''

    compress = True

    def __init__(self, nogroup=False, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.nogroup = nogroup

    def trimmer(self, infiles, outfile):
        '''build trimmer statement on infiles.

        The output is created in exportdir
        '''

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = os.path.join(self.tmpdir_fastq + "trim_galore")
        statement = ["mkdir -p %s;" % tmpdir]
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        # note: not fully implemented
        data_options = ["%(trim_galore_options)s"]

        tmpdir_fastq = self.tmpdir_fastq

        track = P.snip(os.path.basename(outfile), ".fastq.1.gz")
        
        if nfiles == 1:
            infiles = ",".join([self.quoteFile(x[0]) for x in infiles])

            statement.append('''
            trim_galore --rrbs --phred33
            --quality %(quality)s
            --length %(length)s
            --stringency %(stringency)s 
            %(infiles)s
            > %(tmpdir)s/%(track)s.fastq.gz
            2>> %(outfile)s.trim.log
            ''' % locals())

        elif nfiles == 2:
            track1 = track + ".1"
            track2 = track + ".2"
            infiles1 = ",".join([self.quoteFile(x[0]) for x in infiles])
            infiles2 = ",".join([self.quoteFile(x[1]) for x in infiles])

            # statement needs to specify output for both paired ends
            statement.append('''
            trim_galore --rrbs --paired --phred33
            --quality %(quality)s
            --length %(length)s
            --stringency %(stringency)s 
            %(infiles1)s %(infiles2)s
            > %(tmpdir)s/%(track)s.fastq.gz
            2>> %(outfile)s.trim.log
            ''' % locals())
        else:
            raise ValueError(
                "unexpected number read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return " ".join(statement)


class Bismark(object):

    '''maps reads.

    preprocesses the input data, calls mapper and post-process the output data.

    All in a single statement to be send to the cluster.

    Preprocessing function is defined in PipelineMapping.py
    '''

    datatype = "fastq"

    # set to True if you want to preserve colour space files.
    # By default, they are converted to fastq.
    preserve_colourspace = False

    # compress temporary fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    # remove non-unique matches in a post-processing step.
    # Many aligners offer this option in the mapping stage
    # If only unique matches are required, it is better to
    # configure the aligner as removing in post-processing
    # adds to processing time.
    remove_non_unique = False

    def __init__(self):
        pass

    def mapper(self, infiles, outfile):
        '''build trimming statement on infiles.
        '''
        return ""

    def build(self, infiles, outfile):
        '''run mapper.'''

        cmd_mapper = self.mapper(infile, outfile)

        assert cmd_mapper.strip().endswith(";")

        statement = " checkpoint; ".join(cmd_mapper)

        return statement





def splitGeneSet(infile):
    ''' split a gtf file by the first column '''

    last = None
    outfile = None
    outprefix = P.snip(infile, ".gtf.gz")

    for line in IOTools.openFile(infile):

        this = line.split("\t")[0]

        if this == last:
            outfile.write(line)

        else:
            last = this
            if outfile is not None:
                outfile.close()

            outfile = IOTools.openFile("%s.%s.gtf.gz" % (outprefix, this), "w")
            outfile.write(line)
