#!/usr/bin/env bash

#slurms settings






#Rationale: pipeline for Variant to Gene mapping analysis. Output from each analysis: list of genes


#1.Variant annotation
Rscript src/Variant_annotation_FAVOR.Rationale


#2.Quantitative trait loci:
#Genomic boundaries: based on recombination hotspots as found in European ancestry 1000 Genome Project Phase I by ldetect. Source of file:
#https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-all.bed
#rename: input/boundaries_recomb_hotspots_ldetect.txt
Rscript src/boundaries_sugg_sentinels.R
