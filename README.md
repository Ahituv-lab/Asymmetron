===========
# Asymmetron
===========

## Installation

1.- Clone github repository

	git clone https://github.com/Ahituv-lab/Asymmetron

2.- Create a conda virtual enviroment with the necesary dependencies 

	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge

	conda create --name asymmetron pybedtools python=3.5 seaborn numpy

3.- Activate virtual enviroment
	
	conda activate asymmetron

## Introduction
Even though the DNA double helix is a symmetric structure, many biological processes such as replication, transcription and transcription factor binding are directional. The directionality of such processes results in the inhomogeneous distribution of genomic sequences relative to the two complementary DNA strands. Reflecting the directionality biases, strong compositional strand asymmetries have been observed across the entire tree of life, ranging all the way from viral to eukaryotic genomes.

DNA mutations can be oriented relative to transcription and replication, using as reference the template/non-template and leading/lagging strands, respectively. If the reference nucleotide or motif at the site of the mutation is found more frequently in one strand relative to the other, following correction for background strand preferences, it indicates a mutational strand asymmetry. This mutational strand imbalance can have a major impact on disease, development and evolution.

Sequences that are non palindromic can be oriented relative to one another. An exemplar case reflecting the contribution of strand asymmetries in biological processes is the transcription factor CCCTC-binding factor (CTCF). The orientation of the CTCF motif dictates chromatin looping and three dimensional genome topology. Other cases of strand asymmetries include endogenous repetitive elements with preferences in their orientation relative to each other and relative to transcriptional and replicative direction, which can influence their jumping activity.

By studying systematically strand asymmetries we can identify novel DNA elements, improve our comprehension regarding their interactions with one another and advance our understanding regarding the contribution of underlying processes in mutagenesis and evolution. To date, there is no versatile tool to perform analysis of strand asymmetries across biological problems.

## Quick Summary 
Asymmetron is a toolkit for the identifcation of asymmetry patterns in biological sequences. Asymmetron can identify strand asymmetries within consecutive occurrences of a single genomic element and for pairs of overlapping and non-overlapping genomic elements. It can also measure strand asymmetries of genomic elements relative to transcriptional and replicative orientations. Asymmetron can assign strand orientation to third features such as mutations, by orienting them relative to other genomic elements. 

![Schematic_Asymmetron](Resources/Schematic_Asymmetron.png)

## Documentation

An in depth documentation, including a tutorial with examples can be found at:
https://asymmetron.readthedocs.io/en/latest/index.html

## Asymmetron is licensed under the Apache 2.0 License.

## For questions, ideas, feature requests and potential bug reports please contact ilias.georgakopoulossoares@ucsf.edu 
