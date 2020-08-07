.. _contained_asymmetries:
 
=====================
contained_asymmetries
=====================

|

**This function performs strand asymmetry estimation of strand-assigned nucleotide-sequence motifs within strand-assigned regions.**

|

.. list-table:: **Required inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - regions
     - One or more BED-formatted files, containing the regions within which to estimate motif asymmetries.

   * - motifs
     - One or more BED-formatted files, for each of which the asymmetries are calculated.


|


.. list-table:: **Optional inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - --names_regions
     - Assign a name for each of the BED-formatted region files to allow for more human readable output.

   * - --names_motifs
     - Assign a name for each of the BED-formatted motif files to allow for more human readable output.

   * - --orientation_region
     - Orient file(s) relative to annotated BED-formated region file(s) and perform the analysis for the un-annoated file with the new annotations.

   * - --orientation_motif 
     - Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations.

   * - --expected_asym
     - The expected asymmetry bias between the regions and the motifs regarding same or opposite strand orientation. Default is 0.5.

   * - --expected_asym_conv_div
     - The expected convergent / divergent asymmetry bias between the regions and the motifs. Default is 0.5.

   * - --score
     - Optional flag. If provided, assumes the last column of the region files is a scoring metric and uses it to subdivide the analysis into bins.

   * - --bins_score
     - Number of bins to subdivide the score column into. Only runs when --score is provided. Default value is 10.

   * - --plots
     - Returns the associated plots of the asymmetries for each file.


|

**Outputs:**

* Table containing strand asymmetries of ++/+-/-+/-- orientations, same-strand, opposite strand, p-value, p-value with Bonferoni correction.

* Barplot of strand asymmetries for each pair of comparisons.

* Barplot of strand-asymmetries across the score quantiles.

* Table of strand asymmetries for same / opposite and convergent / divergent asymmetries per bin, associated p-vaues and p-values Bonferoni correction.

|

.. note::

   The results from files consisting of palindromic motifs will not show strand asymmetries.


