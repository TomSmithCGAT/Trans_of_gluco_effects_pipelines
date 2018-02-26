======================
Hierachical clustering
======================

The following plots show the clustering of sRNA expression values
after performing a serial alignment and quantification. 

This section presents an analysis to perform unsupervised clustering
of the samples using the distances between their expression profiles
at various classes of sRNA species. If the treatment has a clear
affect on the expression on the sRNA genes, we would expect samples to
cluster by treatment.

miRNA and tRNA annotations were derived from ENSEMBL. piRNA
annotations were obtained from piRBase. sRNA loci quanitifaction was
performed by counting the reads aligned to each feature. miRNA
expression was quantified at two levels, the mature miRNA and the
whole miRNA loci.

From the plots below, it is apparent that the samples do not cluster
by treatment type according to the quantification of expression at any
of the sRNA species examined here. This suggests the treatment does
not have a clear consistent impact of the expression of these sRNA species.

The following plots show the result of hierarchical clustering using
euclidean distance and the Ward clustering criterion. The prefix of
the filename above the plot indicates the sRNA species used for the
clustering. Clustering is shown for the F1 samples only.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: iterative_mapping.dir3/plots.dir/*coloured.dendogram.png
	  
   High resolution plots can be downloaded using the links below
