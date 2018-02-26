================================
Differential expression analysis
================================

The sections on clustering and PCA analysis suggest that the samples
do not cluster by treatment and that the vast majority of the variance
in the data in attriutable to noise, not treatment. 

Nethertheless, it is possible that there are a small number of
differentially expression (DE) sRNA loci in the data. Therefore DESeq2 was
used to identify DE sRNA loci between F1 Dex and F1 Saline. 

Across the sRNA loci tested there were a total of 3 loci identified as
differentially expression (<0.01%). These loci were only identified
for a single aligner. Upon inspection in a genome browser, it's clear
that the 2 loci differentially expressed in the samples aligned by
butter are false positives from binary expression values. Looking
in the F2 samples, the two genes are either highly expressed or not
expressed at all seemingly at random. By chance, they happen to be
expressed in all four F1 Saline samples and no F1 Dex samples. The
other DE gene only just passes the threshold for DE and has a low fold
change (< 2-fold). If the p-value adjustment was performed across all
loci tested rather than within individual sRNA loci, this gene would
not be signficantly DE.

In summary, there are no annotated miRNA, piRNA or tRNA genes which we
can confidently say are differentially expressed in the F1 samples.

The plots below show expression against log-fold change, where
the log-fold change for loci with low expression has been adjusted to
reflect the lower accuracy of quantification. The 3 genes identified
as DE by DESeq2 are highlighted in red.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: DESeq.dir/*_MA_plot.png
	  
   High resolution plots can be downloaded using the links below



