# snakevirome

## Authors
Etienne Loire (CIRAD) & Antony Exbrayat (CIRAD)

## Purpose:

Metagenomic analysis of viral samples

## Steps:

* cleaning
* merging
* filtering
* assembly
* taxonomic annotation 

## Usage:

* edit config.yaml to precise dataset and dependencies path
* edit snakefile to accomodate read files. Currently set to {sample_name}_1.fastq.gz and {sample_name}_2.fastq.gz

## Dependencies

Oh boy... So many
* cutadapt
* bwa
* flash
* megahit
* cap3
* python3
* samtools
* picard tools
* ncbi-blast
* ETE toolkit

