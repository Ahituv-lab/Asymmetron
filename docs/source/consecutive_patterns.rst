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

   * - paths
     - One or more BED-formatted files.



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

   * - --simulations
     - Number of simulations to perform from which the empirical p-value is derived. Default is N=100.

   * - --plots
     - Returns the associated plots of the asymmetries for each file.


|

**Outputs:**

* Table of strand asymmetries for all patterns for each motif.
* BED files with the statistically significant consecutive regions. One file for each pattern for each file inputted.
* Barplots of expected and observed consecutive occurrences of each of the patterns.
* Plots showing the distribution of consecutive occurrences for each pattern.
* Heatmap of subdivisions of the signal in distances across all patterns if --bins is selected.

|

.. note::

   The results from files consisting of palindromic motifs will not show strand asymmetries.

| 

.. note::

   * The default strand asymmetry estimation is obtained as shown in the tutorial of this function. However, when the option for custom patterns is selected, then the following procedure is followed to estimate strand asymmetries.

   We define as p the probability that the motif is found in strand A\ :sub:`n`\, calculated as the number of appearances of the motif in strand A\ :sub:`n`\ divided by the total number of appearances of the motif in the genome. Similarly, we define qas the probability that the motif is found in strand B\ :sub:`n`. Since there are only two strands, p=1-q. the probability that any given pattern would emerge by chance is estimated based on the number of + or - signs it contains. A pattern that contains u “+” signs and v “-” signs, should appear with a probability p\ :sub:`u`\*q\ :sub:`v`\. We also apply the Bonferroni correction, based on the total number of motif occurrences. We can thus determine the p-value of any given pattern emerging. 
