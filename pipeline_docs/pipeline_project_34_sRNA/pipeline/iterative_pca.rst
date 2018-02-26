============================
Principal Component Analysis
============================

The following section presents principal components analysis (PCA) of the
sRNA-Seq expression data after performing the serial alignment and quantification. 

For an explanation of PCA, please see the previous section describing
PCA following quantification from genome-wide alignment.

The samples do not seperate by treatment group in the first four
principal components for any of the sRNA species examined.
It appears then that the vast majority of the variance in the data in noise
and that the treatment has little or no effect on the expression of
the sRNA species examined

The following plots show the PCA results. The prefix of the filename
above the plot indicates the sRNA species. PCA was performed sperately
on the F1 samples 

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: iterative_mapping.dir3/plots.dir/*pca_PC*.png
	  
   High resolution plots can be downloaded using the links below

The following plots show the variance explained for each PC.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: iterative_mapping.dir3/plots.dir/*pca_variance_explained.png
	  
   High resolution plots can be downloaded using the links below


