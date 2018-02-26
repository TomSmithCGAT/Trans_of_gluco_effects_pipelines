==========================
Peak calling intersections
==========================

Three peakcallers have been used to call peaks in the H3K4me3, H3K9me3
and H3K27me3 samples.

- `macs2 <http://liulab.dfci.harvard.edu/MACS/00README.html>`_ This
  peak caller has been optimised for narrow peak calling

- `sicer
  <http://www.genomatix.de/online_help/help_regionminer/sicer.html>`_
  This peak caller is particularly useful for histone modification
  peaks and has been run here in "narrow" and "broad" peak detection
  modes

- `peakranger <http://ranger.sourceforge.net/>`_ This peak caller has
  a "narrow" peak caller (ranger) and a "broad" peak caller
  (ccat). Both have been used here

To assess the reproducibility of peak calling across the three
replicates, the peak profiles were intersected. Where two peaks
overlapped by at least 1bp in two samples, the peak was considered to
be identified in both samples. The venn diagrams below show the
overlap of peaks between the three replicates

Sicer in narrow mode appears to call reasonably consistent peaks
for all marks. Sicer broad is a slight improvement in term of the
fraction of peaks which are reproduced across all 6 samples per
tri-methylation mark.

Peakranger (ccat and ranger) is much less consistent than sicer.

Macs2 is even worse for consistency.

**Sicer is the most suitable peakcaller for our IP**.

sicer - narrow
--------------


.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: intersections.dir/*_sicer_narrow_peak_venn.png
	  
   High resolution plots can be downloaded using the links below


sicer - broad
--------------
.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: intersections.dir/*_sicer_broad_peak_venn.png
	  
   High resolution plots can be downloaded using the links below


Peakranger - ccat
-----------------
.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: intersections.dir/*_ccat_region_peak_venn.png
	  
   High resolution plots can be downloaded using the links below


Peakranger - ranger
-------------------
.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: intersections.dir/*_ranger_region_peak_venn.png
	  
   High resolution plots can be downloaded using the links below

macs2
-----

From the venn diagrams it's clear that macs2 is calling far fewer
peaks and these are not reproducible across the replicates

.. report:: project34Report.imagesTracker
   :render: gallery-plot
   :glob: intersections.dir/*-macs2_peaks_venn.png
	  
   High resolution plots can be downloaded using the links below
