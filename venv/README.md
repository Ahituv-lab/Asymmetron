===========
# Asymmetron
===========

Copyright 2020. All rights reserved.

## Installation


## Summary 

Asymmetron is a toolkit for the identifcation of asymmetry patterns in biological sequences. Asymmetron can identify strand asymmetries within consecutive occurrences of a single genomic element and for pairs of overlapping and non-overlapping genomic elements. It can also measure strand asymmetries of genomic elements relative to transcriptional and replicative orientations. Asymmetron can assign strand orientation to third features such as mutations, by orienting them relative to other genomic elements. 

It encompasses four functions:
1.	consecutive_patterns.py:	Estimates the asymmetry biases within consecutive occurrences of a single motif.
2.	contained_asymmetries.py	Estimates the asymmetry biases of a motif within an encompassing region.	
3.	pairwise_asymmetries.py		Estimates the asymmetry biases between two motifs.
4.	orientation.py			Orients an un-annotated BED file relative to overlapping instances of a second, annotated BED file.

## Instructions for each function

### 1. consecutive_patterns.py
#### Input requirements:
	-motifs: One or more BED-formatted files
#### Optional inputs:
	-names: The name or each of the inputted BED-formatted files
	-min_dist: Minimum distance of consecutive occurrences to consider in the analysis. Default is 0.
	-max_dist: Maximum distance of consecutive occurrences to consider in the analysis. Default is 100.
	-patterns: Patterns to search, comma separated. Default is ++,--,+-,-+.
	-orientation: Orient file(s) relative to annotated BED-formated file(s) and perform the analysis for the un-annotated file with the new annotations.
	-bins: Number of bins to subdivide the analysis in. Default is 1, which does not perform this analysis.
	-threshold: Threshold of p-value of consecutive patterns to save in new BED file.
	-plots: Returns the associated plots of the asymmetries for each file.
#### Outputs:
	-Table of strand asymmetries for all patterns for each motif.
	-BED files with the statistically significant consecutive regions. One file for each pattern for each file inputted.
	-Histogram of expected and observed distribution of consecutive patterns.

### 2. contained_asymmetries.py
#### Input requirements:
	-regions: One or more BED-formatted files, containing the regions within which to estimate motif asymmetries.
	-motifs: One or more BED-formatted files, for each of which the asymmetries are calculated.
#### Optional inputs:
	-names_A: The name of each of the inputted BED-formatted region files.
	-names_B: The name of each of the inputted BED-formatted motif files.
	-upstream: Only look for occurrences of motif A upstream of motif B
	-downstream: Only look for occurrences of motif A downstream of motif B
	-orientation_region: Orient file(s) relative to annotated BED-formated region file(s) and perform the analysis for the un-annoated file with the new annotations.
	-orientation_motif: Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations.
	-expected_asym: The expected asymmetry bias between the regions and the motifs. Default is 0.5.
	-expected_asym_conv_div: The expected convergent / divergent asymmetry bias between the regions and the motifs. Default is 0.5.
	-score: For region files, uses the last column to subdivide the analysis into quartiles. Default is not to perform this.
	-quartiles: Number of quartiles to subdivide the score into. Only runs when --score is provided. Default value is 10.
	-plots: Returns the associated plots of the asymmetries for each file.
#### Outputs:
	-Table containing strand asymmetries of ++/+-/-+/-- orientations, same-strand, opposite strand, p-value, p-value with Binomial correction.
	-Histograms of strand asymmetries for each pair of comparisons.
	-Barplot of strand-asymmetries across the quartiles.
	-Barplot of strand-asymmetries across the bins.

### 3. pairwise_asymmetries.py
#### Input requirements:
	-motifsA: One or more BED-formatted files
	-motifsB: One or more BED-formatted files
#### Optional inputs:
	-names_A: The name of each of the inputted BED-formatted motif A files.
	-names_B: The name of each of the inputted BED-formatted motif B files.
	-orientation: Orient file(s) relative to annotated BED-formated motif file(s) and perform the analysis for the un-annoated file with the new annotations.
	-expected_asym: The expected asymmetry bias between the pairs of motifs. Default is 0.5.
        -expected_asym_conv_div: The expected convergent / divergent asymmetry bias between the pairs of motifs. Default is 0.5.
	-min_dist: Minimum distance to consider in the analysis. Default is 0.
	-max_dist: Maximum distance to consider in the analysis. Default is 100.
	-bins: Number of bins to subdivide the analysis in. Default is 1, which does not perform this analysis.
	-plots: Returns the associated plots of the asymmetries for each file.
#### Outputs:
	-Table containing strand asymmetries of ++/+-/-+/-- orientations, same-strand, opposite strand, p-value, p-value with Binomial correction.
	-Histograms of strand asymmetries for each pair of comparisons.
	-Barplot of strand-asymmetries across the bins.

### 4. orientation.py 	
#### Input requirements:
	-motif_no_annotation
	-motif_annotation
#### Output:
	-Un-annotated BED formatted motif with annotation based on -motif_annotation file.
