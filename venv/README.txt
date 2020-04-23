===========
Asymmetron
===========

Copyright 2020. All rights reserved.
Repository:

Summary: Asymmetron is a toolkit for the identifcation of asymmetry patterns in biological sequences.
It encompasses four functions:
1.	consecutive_patterns.py:	Estimates the asymmetry biases within consecutive occurrences of a single motif.
2.	contained_asymmetries.py	Estimates the asymmetry biases of a motif within an encompassing region.	
3.	pairwise_asymmetries.py		Estimates the asymmetry biases between two motifs.
4.	orientation.py			Orients an unnanotated BED file relative to overlapping instances of a second, annotated BED file.


1. consecutive_patterns.py
Input requirements:
	-motifs: One or more BED-formatted files
Optional inputs:
	-names: The name or each of the inputted BED-formatted files
	-min_dist: Minimum distance of consecutive occurrences to consider in the analysis. Default is 0.
	-max_dist: Maximum distance of consecutive occurrences to consider in the analysis. Default is 100.
	-patterns: Patterns to search, comma separated. Default is ++,--,+-,-+.
4.	
Input requirements:


