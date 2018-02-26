===========================
H3 modification enrichments
===========================

This section presents an analysis of the enrichment of IP signal in
various annotated features. The raw data used in this analysis is the
context results section from the mapping pipeline.

Only features to which at least 1/10000 of the mapped reads aligned
are included in the analysis.

The first part of the analysis shows the enrichment/depletion of IP
signal for the three H3 tri-methylation marks across annotated
features in the rat genome. This is followed by a table of t-test
results for all features and all modifications. The intention here is
to identify genomic features at which the enrichment/depletion of IP
signal vs H3 Input is significantly different between Saline and Dex
samples. **There are no significant differences between Saline and Dex**.

The second part is a principal components analysis using the
enrichments scores in each sample. The intention here is to examine
how much of the variability in the enrichment scores can be attributed
to expected factors such as the difference between the H3
tri-methylation marks and the difference between Saline and Dex
treatment groups. **The majority of the variation appears to relate to
the litter number**. **A minor proportion of the variation appears to
relate to the difference between Saline and Dex**.

The third part is a blocked anova analysis to identify significantly
different enrichment scores for the Saline and Dex IPs, with litter
number as the blocking factor. This analysis is motivated by the
principal component analysis which indicates that the majority
of the variation in the enrichment scores originates from the litter
number or a latent variable which is confounded in the litter
number. The inclusion of the litter number as a blocking factor will
enable detection of significant differences after accounting for the
effect of litter number. **There are no signficant differences between
enrichment scores for Saline and Dex after accounting for variation
due to litter number at a 5% False Discovery Rate (FDR).However,
there are 3 significant features at a very relaxed 20% FDR**.


Enrichment vs H3 Input at annotated features
--------------------------------------------
The plot below shows the enrichment of IP signal vs H3 Input in
annotated features for the 3 tri-methylation modifications. The
features shown are restricted to include only those features in which
the absolute maximum log2-fold enrichment vs. H3 Input across all
samples is > 0.5. There is no clear difference between Saline and Dex
at any of the three features shown, although the enrichment at CpG
islands is slightly higher, and the depletion at SINE "ID" elements is
slightly lower for all 3 tri-methylation modifications.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment.png
	  
   High resolution plots can be downloaded using the links below


T-test for significant difference in enrichment
-----------------------------------------------
The table below shows the the mean enrichment values (e.g log2 H3K4me3
vs. H3 Input), the standard error (se) and the results of performing a
two-tailed student's t-test testing the null hypothesis that the
enrichment is the same in Saline and Dex. The test is performed for
each tri-methylation mark at each feature seperately. p-values were
then adjusted using the Benjamini & Hochberg proceedure to obtain the
False discovery rate (FDR). We can see from the table that there are
no significant differences after p-value adjustment (using 5% FDR
cut-off). Even before adjustment, there are only 3 p-values below 0.05
(all H3K27me3). Given we are performing 84 tests, this is actually
slightly less than we would expect by chance (4.2). Furthermore, all 3
features are very slightly depleted in the IP signal and the
difference between mean enrichment in Saline and Dex is very slight.



.. report:: enrichment.EnrichmentTTest
   :render: table
   :large: xls
   :force:

   The table can be downloaded using the links below

Principal components analysis
-----------------------------

The plots below show the results of a principal components analysis to
examine the variance within the enrichment scores at each feature for
each sample. Ideally, most of the variability should originate from
expected factors such as the difference between the three
tri-methylation modifications and the difference between the treatment
groups. 

The principal components (PC) represent orthogonal transformations
of the original set of observations (enrichment scores) ordered such that the
first PC has the largest possible variance. The first PC will
therefore align with the source of greatest variability between the
samples. We can visualise the results by plotting the samples in PC
space to show how they seperate on the PCs. 

The first plots show a PC analysis of all samples together. The
remaining plots show the same analysis but conducted on each
tri-methylation modfication seperately. The rationale for this is that
we aren't particularly neccessarily expecting differences between
Saline and Dex to be the same for all three tri-methylation
modifications and what we really want to see is a difference between
Saline and Dex at each modification seperately.

A summary of the results is included at the top of each section


All samples
-----------

In the first plot below we can see that the first and second PCs account
for 78.8% of the variance and appear to seperate the H3K4me3 samples
somewhat from the H3K9me3 and H3K27me3 samples. This indicates that
the major source of variance within the enrichment scores originates from the
difference between the tri-methylation modifications. This is of
course not suprising. 

Plotting the 3rd-6th PC we can see that no PC seperates the samples by
treatment group, suggesting Dex does not have a consistent impact on
the enrichment scores for the three tri-methylation
modifications. This again is perhaps to be expected. The next section
therefore examines each modification in turn.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_all_PC*.png

   High resolution plots can be downloaded using the links below	  

H3K4me3
-------

The first and second PCs account for 73.2 % of the variance in the
enrichment scores for H3K4me3 IP and appear to seperate the
samples by "litter number". This is suprising given that the litter
number for Saline and Dex samples does not link them in any way
biologically but was used to number the litters in turn during the
experiments. The litter number may therefore be a proxy for some other
experimental "batch effect" variable which is introducing
variability in the enrichment scores. This kind of variable is usually
referred to as a latent or hidden variable

The 3rd and 4th PCs account for 23.6% of the variance in the
enrichment scores for H3K4me3 samples and appear to seperate the
samples by treatment group (i.e Dex and Saline). This suggests that a
minor proportion of the variation is due to treatment but that this is
less than the variation due to the latent variable(s). Removal of
this latent variable(s) may therefore be required to identify
differences between Saline and Dex.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_H3K4PC1_PC2*.png

   PC1 and PC2.	 

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_H3K4PC3_PC4*.png

   PC3 and PC4.
   High resolution plots can be downloaded using the links below

H3K9me3
-------

The results for the PCA of H3K9me enrichment scores is similar to the
above H3K4me3 analysis.  The first PC accounts for 75.5% of the
variance in enrichment scores for the H3K9me3 IP and again appears to
seperate the samples by a "litter number". The 2nd PC accounts for
12.9% of the variance and seperates the samples by treatment group.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_H3K9PC1_PC2*.png

   PC1 and PC2.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_H3K9PC3_PC4*.png

   PC3 and PC4.
   High resolution plots can be downloaded using the links below

H3K27me3
--------
The results for the PCA of H3K9me enrichment scores is very similar to the
above H3K9me3 analysis.  The first PC accounts for 72.4% of the
variance in enrichment scores for the H3K9me3 IP and again appears to
seperate the samples by a "litter number". The 2nd PC accounts for
17.3% of the variance and seperates the samples by treatment group.


.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_H3K27PC1_PC2*.png

   PC1 and PC2.

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_pca_H3K27PC3_PC4*.png

   PC3 and PC4.
   High resolution plots can be downloaded using the links below


Blocked ANOVA
-------------
The table below shows the results of a blocked-design ANalysis Of
Variance Analysis (ANOVA) on the enrichment scores for the 3
tri-methylation marks at annotated features. The intention here is to
identify features with signficantly different enrichment in Saline and
Dex samples. The litter number is included as a blocking factor to
account for the previously noted considerable proportion of the
variation which appears to originate from differences between the
litter numbers.

After adjusting for multiple testing, there are no
significantly different enrichments at a 5% FDR. However, there are 3
features which are significant if we relax the FDR cut-off to
20%. These might be worth further consideration on the
understanding that they are quite likely to be false positives. 

Below the table are plots of the enrichment scores for these three
features at all tri-methylation marks. 

As we can see, the enrichment values for the DNA transposon family
hAT_Tip100 are close to zero and although the difference between
Saline and Dex is significant at the 20% FDR cut-off for H3K27me3, the
range of enrichment values encompasses zero for both Saline and
Dex. Similarly, the significant differences between Saline and Dex for
the transposon-like B2 repeats for H3K4me3 only relates to an
apparently slightly more depleted IP signal vs. H3 Input in the Dex
samples at these repeats. Neither of these results seem particularly
worth following up.

The significant difference between Saline and Dex (at 20% FDR) at ERV1
elements represents a difference in enrichment of IP signal which is
perhaps more likely to be validated. However, again the difference is
very slight and the enrichment of IP vs. H3 input is less than
1.3-fold (2^0.35) in all cases.

.. report:: enrichment.EnrichmentBlockedAnova
   :render: table
   :large: xls
   :force:

   The table can be downloaded using the links below


.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: enrichment.dir/enrichment_blockedANOVA_*png

   High resolution plots can be downloaded using the links below
