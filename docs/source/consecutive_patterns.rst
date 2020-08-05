.. _consecutive_patterns:
  
====================
consecutive_patterns
====================

|

**This function enables the identification of patterns in the orientation of consecutive occurrences of a feature.**


|


.. list-table:: **Required inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - --paths
     -  paths: One or more BED-formatted files.



|



.. list-table:: **Optional inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - --names
     - Assign a name for each of the inputted BED-formatted files to allow for more human readable output.

   * - --min_distance
     - Minimum distance between consecutive occurrences to consider in the analysis. Default is 0.

   * - --max_distance
     - Maximum distance between consecutive occurrences to consider in the analysis. Default is 100.

   * - --patterns   
     - Patterns to search, comma separated. Default is same / opposite strand orientation analysis. 

   * - --orientation
     - Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations.

   * - --bins
     - Number of bins to subdivide the analysis in. Default is 1 which does not perform this analysis.

   * - --threshold
     - Threshold of p-value of consecutive patterns to save in new BED file.

   * - --plots
     - Returns the associated plots of the asymmetries for each file.


|

**Outputs:**

* Table of strand asymmetries for all patterns for each motif.
* BED files with the statistically significant consecutive regions. One file for each pattern for each file inputted.
* Barplots of expected and observed consecutive occurrences of each of the patterns.
* Plots showing the distribution of consecutive occurrences for each pattern.
* Heatmap of subdivisions of the signal in distances across all patterns if --bins is selected.
