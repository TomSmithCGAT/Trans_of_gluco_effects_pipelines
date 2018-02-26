======================
Hierachical clustering
======================

This section presents an analysis to perform unsupervised clustering
of the samples using the distances between their expression profiles
at various classes of sRNA species. If the treatment has a clear
affect on the expression on the sRNA genes, we would expect samples to
cluster by treatment.

miRNA and tRNA annotations were derived from ENSEMBL. piRNA
annotations were obtained from piRBase. sRNA loci quanitifaction was
performed with featureCounts before normalisation by library size and
removal of low abundant sRNA loci.

Two aligners were used to map the sRNA reads onto the reference genome,
butter and bwa. Butter is an iterative mapper designed for sRNA-Seq
analysis which attempts to align reads which map to multiple locations
to their true single location. To do so it uses the
alignments of longer reads which are more likely to map to a single
genomie loci to inform the alignment of shorter reads. BWA on the other
hand will report a single location at random if more than one equally
good alignment is found. 

From the plots below, it is apparent that the samples do not cluster
by treatment type according to the quantification of expression at any
of the sRNA species examined here. This suggests the treatment does
not have a clear consistent impact of the expression of these sRNA species.

The following plots show the result of hierarchical clustering using
euclidean distance and the Ward clustering criterion. The prefix of
the filename above the plot indicates the sRNA species used for the
clustering. Clustering was performed sperately on all samples and the
F1 samples in isolation as specified in the filename ("F1" or
"all"). The filename also contains the name of the aligner ("butter"
or "bwa") used


.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: clustering.dir/*dendogram.png
	  
   High resolution plots can be downloaded using the links below
