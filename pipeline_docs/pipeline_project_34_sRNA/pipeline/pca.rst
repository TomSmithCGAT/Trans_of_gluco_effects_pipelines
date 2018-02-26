============================
Principal Component Analysis
============================

This section presents principal components analysis (PCA) of the
sRNA-Seq. Principal components are orthogonal transformations of the
original data. The principal components (PCs) are ordered such that
the first PC contains the largest possible amount variance and each
following PC and the largest possible variance under the constraint
that it is orthogonal to the previous PCs.

Two aligners were used to map the sRNA reads onto the reference
genome, butter and bwa. Butter is an iterative mapper designed for
sRNA-Seq analysis which attempts to align reads which map to multiple
locations to their true single location. To do so it uses the
alignments of longer reads which are more likely to map to a single
genomie loci to inform the alignment of shorter reads. BWA on the
other hand will report a single location at random if more than one
equally good alignment is found.

miRNA and tRNA annotations were derived from ENSEMBL. piRNA
annotations were obtained from piRBase. sRNA loci quanitifaction was
performed with featureCounts before normalisation by library size and
removal of low abundant sRNA loci.

For our purposes, it is hoped that one of the top PCs will seperate
the treatment groups, which would indicate that a large proportion of
the variance in the data is due to the differences between the
treatment groups. However, when comparing just the F1 samples, the
samples do not seperate by treatment on PC1-PC4 (as shown below), or
indeed any of the following PCs. Furthermore, in most instances the
amount of variance explained by each PC decreases steadily (see
variance line plots below), suggesting the PCs are seperating the
samples largely by noise. The one execption to this is the PCA for
tRNA loci when using the aligner butter in which PC1 explains more
than 2/3 of the variance. However PC1 still does not seperate the
treatment groups and the same analysis using the aligner bwa does not
show the same results. This suggests the observed seperation on PC1
for tRNA loci when using butter may be an artifact of the underlying
assumption of this aligner that multi-mapping sRNA reads originate
from a single genomic loci when in fact they originate from multiple
tRNA loci with very similar sequences.

When all samples are included in the analysis, the first two PCs
explain a much greater proportion of the variance and seperate two F2
DEX samples (1 and 4) from the rest of the samples. However, there is
no seperation of the treatment groups on these PCs or later PCs.

In conclusion, it appears that the vast majority of the variance in
the data in noise and that the treatment has little effect on the
expression of the sRNA species examine

The following plots show the PCA results. The prefix of the filename
above the plot indicates the sRNA species. PCA was performed sperately
on all samples and the F1 samples in isolation as specified in the
filename (“F1” or “all”). The filename also contains the name of the
aligner (“butter” or “bwa”) used.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: clustering.dir/*pca_PC*.png
	  
   High resolution plots can be downloaded using the links below

The following plots show the variance explained for each PC.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: clustering.dir/*pca_variance_explained.png
	  
   High resolution plots can be downloaded using the links below


