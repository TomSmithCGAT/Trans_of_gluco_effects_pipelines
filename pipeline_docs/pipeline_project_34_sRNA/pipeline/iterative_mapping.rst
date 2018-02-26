===============
Serial Mapping 
===============

Quantification of small RNA species abundance from sRNA-Seq can
broadly be performed in three ways, supervised semi-supervised or
unsupervised. In the semi-supervised and unsupervised approaches, the
reads are aligned back to the reference genome. Estimates of small RNA
species abundance is then performed either using the genomic loci of
known annotations (semi-supervised) or using *de-novo* annotations
(supervised). In the supervised approach, reads are aligned directly
to the sequences of the known annotations.

Alignment to the genome (semi-supervised/unsupervised) avoids biasing
the read alignments towards known annotations. The downside of aliging
to the genome is that reads may erroneously align to repeative
elements with very similar sequences to sRNAs.

The supervised approach makes the assumption that the vast majority of
the sRNA-Seq reads will have originated from the sRNA species and
avoids mistakenly aligning reads from sRNA species to elsewhere in the
genome. However, the supervised approach may artifically "force" reads
to align to the wrong sRNA species if the sRNA annotation is
incomplete. 

The initial analysis was performed in a semi-supervised manner using
annotations from ucsc (tRNAs, rRNAs), miRBase (miRNAs) and piRBank
(piRNAs). 

This section presents an identical analysis based on a supervised
quantification. The supervised approach was developed in conversation
with Oliver Rando.

Reads were serially aligned to rRNAs, tRNAs, miRNAs, piRNAs and other
small RNA species in turn. The alignment was also performed serially
with regards to the number of "best" matches found. In the first round
of alignment, only unique alignments were retained. Reads with 2
"best" alignments were then assigned to a sRNA based on the density of
the previous alignments. This process was then repeated for reads
aligning to 3 locations and so on. Reads that aligned to more than 30
locations were discarded.

Quantification of 5' tRNA fragments was performed by collapsing
individual tRNA sequences based on the 30nt 5' sequence. Most tRNAs
possess unique 5' sequences however there are 26 Cys(TGY) tRNAs which
share the same 30nt sequence. In these cases it is not possible to
determine which tRNA the reads originate from and quantification is
best expressed as the total from all tRNAs with the same 30nt 5'
sequence.

piRNA and miRNA expression was estimated from the number of reads
assigned to each gene. Mature miRNA expression was estimated by taking
the maximum coverage across the miRNA hairpin.

The sections below present the analysis of the small RNA
quantification estimates following serial alignment.

In all sections the analysis is performed seperately for miRNAs,
mature miRNAs, piRNAs and tRFs (tRNA fragments)

--------
Analyses
--------

.. toctree::
   :maxdepth: 3

   iterative_clustering.rst
   iterative_pca.rst
   iterative_deseq2.rst
