==================
Gene counts tables
==================

Gene counts are provided in tabulated form below.
Samples are labelled as Tissue_Treatment_Replicate

Note: The F2 samples showed very variable distributions of sRNA sizes,
suggestive of extreme degredation in some sample. All analysis has
therefore been restricted to the F1 samples.

For all tables, the counts given are the estimated number of reads per sample
which align to each feature (gene). This is an estimate as some reads
will align to multiple features and some reads will have been ignored
for 

miRNA loci were obtained from Ensembl78
tRNA loci were obtained from UCSC
piRNA loci were obtained from piRBase

The gene counting was performing in two ways:

1. Genome alignment:
   Reads aligned to the genome, followed by counting of
   reads falling into genomic features

2. Iterative alignment to features:
   Iterative alignment of reads to features a la Oliver Rando

In addition, the genome alignment was performed with two mappers (bwa
and butter). Only the bwa counts tables are provided here as butter
alignment was sub-optimal.

Counts tables for the genome alignment and iterative alignment are
provided seperately.

tRNA loci show considerable sequence similarity which complicates
assignment of reads. For this reason, tRNA loci were
merged by codon for the tRNA-fragment (tRF) analysis when using the
first method. To quantify tRFs using the second method, tRFs were
merged based on the first 30 bases (from the 5' end). The gene name
given here is the 30 bp sequence.  A seperate file called
tRNA_sequences.tsv maps these sequences onto the tRNAs in which the
sequence is found.


Counts by genome alignment
==========================

Table of counts at tRNA fragment loci. Counts are given per codon as
explained above.

.. report:: tables.tRFCountsGenome
   :render: xls-table
   :force:

Table of counts at miRNA loci.

.. report:: tables.miRNACountsGenome
   :render: xls-table
   :force:

Table of counts at piRNA loci.

.. report:: tables.piRNACountsGenome
   :render: xls-table
   :force:


Counts by iterative alignment direct to features
================================================

Table of counts at tRNA fragment loci. Counts are given per 30 bp 5'
sequence as explained above.

.. report:: tables.tRFCountsIterative
   :render: xls-table
   :force:

Table which maps tRNA 30 bp 5' sequences to tRNA loci

.. report:: tables.tRFmapSequenceLoci
   :render: xls-table
   :force:

Table of counts at mature miRNA loci. Counts represent the maximum
coverage depth at the locus

.. report:: tables.matureMiRNACountsIterative
   :render: xls-table
   :force:

Table of counts at putative piRNA loci. Counts represent the maximum
coverage depth at the locus

.. report:: tables.piRNACountsIterative
   :render: xls-table
   :force:

All the tables above can be downloaded by right-clicking on the links below
