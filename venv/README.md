===========
# Asymmetron
===========

Copyright 2020. All rights reserved.
Repository:

## Installation

##Summary: Asymmetron is a toolkit for the identifcation of asymmetry patterns in biological sequences.
It encompasses four functions:
1.	consecutive_patterns.py:	Estimates the asymmetry biases within consecutive occurrences of a single motif.
2.	contained_asymmetries.py	Estimates the asymmetry biases of a motif within an encompassing region.	
3.	pairwise_asymmetries.py		Estimates the asymmetry biases between two motifs.
4.	orientation.py			Orients an unnanotated BED file relative to overlapping instances of a second, annotated BED file.


### 1. consecutive_patterns.py
Input requirements:\
	-motifs: One or more BED-formatted files\
Optional inputs:\
	-names: The name or each of the inputted BED-formatted files\
	-min_dist: Minimum distance of consecutive occurrences to consider in the analysis. Default is 0.\
	-max_dist: Maximum distance of consecutive occurrences to consider in the analysis. Default is 100.\
	-patterns: Patterns to search, comma separated. Default is ++,--,+-,-+.\
	-bins: Number of bins to subdivide the analysis in. Default is 1, which does not perform this analysis.\
	-plots: Returns the associated plots of the asymmetries for each file.\
	

### 2. contained_asymmetries.py
Input requirements:\
	-regions: One or more BED-formatted files, containing the regions within which to estimate motif asymmetries.\
	-motifs: One or more BED-formatted files, for each of which the asymmetries are calculated.\
Optional inputs:\
	--names_A: The name of each of the inputted BED-formatted region files.\
	--names_B: The name of each of the inputted BED-formatted motif files.\
	--expected_asym: The expected asymmetry bias between the regions and the factors. Default is 0.5.\
	--expected_asym_conv_div: The expected convergent / divergent asymmetry bias between the regions and the factors. Default is 0.5.\
	--score: For region files, uses the last column to subdivide the analysis into quartiles. Default is not to perform this.\
	--quartiles: Number of quartiles to subdivide the score into. Only runs when --score is provided. Default value is 10.\
	--plots: Returns the associated plots of the asymmetries for each file.\

### 3. pairwise_asymmetries.py
Input requirements:\
	-motifsA: One or more BED-formatted files\
	-motifsB: One or more BED-formatted files\
Optional inputs:\
	



### 4. orientation.py 	
Input requirements:


