.. _pairwise_asymmetries:


====================
pairwise_asymmetries
====================


|


**This function performs strand asymmetry estimation between two strand-assigned motifs in proximity or overlapping each other.**

|


.. list-table:: **Required inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - motifsA
     - One or more BED-formatted files

   * - motifsB 
     - One or more BED-formatted files


|

.. list-table:: **Optional inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - --names_A
     - The name of each of the inputted BED-formatted motif A files.

   * - --names_B
     - The name of each of the inputted BED-formatted motif B files.

   * - --upstream_only
     - Only look for occurrences of motif A upstream of motif B. Incompattible with -downstream.

   * - --downstream_only
     - Only look for occurrences of motif A downstream of motif B. Incompatible with -upstream.

   * - --orientation
     - Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations.

   * - --expected_asym
     - The expected asymmetry bias between the pairs of motifs regarding same or opposite strand orientation. Default is 0.5.

   * - --expected_asym_conv_div
     - The expected convergent / divergent asymmetry bias between the pairs of motifs. Default is 0.5.

   * - --min_distance
     - Minimum distance to consider in the analysis. Default is 0.

   * - --max_distance
     - Maximum distance to consider in the analysis. Default is 100.

   * - --bins
     - Number of bins to subdivide the analysis in. Default is 1, which does not perform this analysis.

   * - --plots
     - Returns the associated plots of the asymmetries for each file.


|


**Outputs:**

* Table containing strand asymmetries of ++/+-/-+/-- orientations, same-strand, opposite strand, p-value, p-value with Binomial correction (corrected for multiple file queries). 

* Barplots of strand asymmetries for each pair of comparisons for same versus opposite and convergent versus divergent strand orientations.

* Barplot of strand-asymmetries across the bins for same versus opposite and convergent versus divergent strand orientations.


