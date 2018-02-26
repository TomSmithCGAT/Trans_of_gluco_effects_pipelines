==================================
Differential histone modifications
==================================

So far, a single tool has been used to examine differential histone
modifications between Saline and Dex treatments.  `MMDiff
<http://www.biomedcentral.com/1471-2164/14/826>`_ (developed by Guido
Sanguinetti's group at the School of Informatics, University of
Edinburgh) identifies differences between ChIP-Seq peaks by examining
the shape of the peak profiles. The advantage of this approach is that
it can identify both changes in the number of aligned reads (a proxy
for binding strengh of a transcription factor, or the histone
modification frequency) and changes in peak profile which may also
represent functional differences.

The difference between two peak profiles is expressed as a distance
metric. This distance metric is calculated between all pairs (between
replicate and between groups). For each peak, the average distance is
calculate between all pairs across the groups tested (Saline and
Dex). This is then compared to the average distance between replicates
for peaks with a similar number of counts to estimate the
p-value. This p-value is then adjusted for multiple testing and all
peaks with an adjusted p-value <0.05 are identified as having
significantly different peak profiles between the two groups.

For this analysis I used all peaks identified in at least three samples
across the 6 samples for a single histone modification. Only the peaks
called by Sicer in narrow and broad mode were considered here as the
other two peak callers used were not consistent, as described in the
previous section.

MMDiff will then compare the profile of reads in the
3 Saline samples and the 3 Dex samples for the genomic loci
surrounding the peak

The plots below show the calculate distances within and between
groups. Significant differences are highlighted in red.

**There are no signficant differences between Saline and Dex**

sicer - narrow
--------------
.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: MMDiff.dir/*_sicer_narrow_mmdiff.png
	  
   High resolution plots can be downloaded using the links below


sicer - broad
--------------
.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: MMDiff.dir/*_sicer_broad_mmdiff.png
	  
   High resolution plots can be downloaded using the links below





