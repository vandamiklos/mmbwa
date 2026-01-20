# mmbwa

A command line tool that improves alignment accuracy of highly fragmented variants while minimizing runtime 
using minimap2 and bwa-mem.

## Table of contents
* [Overview](#overview)
* [Install](#install)
* [Quick usage](#quick-usage)
* [Requirements](#requirements)
* [Detailed usage](#detailed-usage)
* [Options](#options)
* [Benchmark](#benchmark)

## Overview
### Start from FASTQ files
Reads are first aligned with minimap2, the output is piped into a filter script that selects 
reads with long soft-clips and directs them to bwa-mem.
### Start from minimap2 alignment BAM/SAM file
Optionally minimap2 alignment BAM/SAM files can be used as input, in this case the filtering and bwa-mem alignment
steps are run.

## Install

### PIP
```bash
pip install -e .
```

### Conda

in progress

## Requirements
Python 3.7+

samtools, bedtools, pysam, minimap2, bwa-mem

## Quick usage
Default settings correspond to Nanopore ONT reads.
### Align with minimap2, filter, and re-align with bwa-mem
```bash
mmbwa ref_genome --input-fq /path/to/file.fq --output /path/to/outdir
```

### Use pre-generated minimap2 alignment (SAM/BAM), filter and realign with bwa-mem
```bash
mmbwa ref_genome --input-aln /path/to/file.bam --output /path/to/outdir
```

## Detailed usage
### Soft-clip threshold
The `--threshold` flag controls which primary alignments are filtered and sent to re-alignment.
It is calculated by the fraction of soft-clip length and primary alignment length.
The higher this fraction is the less likely it is that a read is re-aligned.

### Change minimap2 parameters

### Change bwa-mem parameters

### Exclude genomic regions

### Re-align unmapped reads

### Keep temporary files

### Sort and index final output




## Options

* ref: Reference genome FASTA file, required argument.
* --input-fq: Input FASTQ file, required if no minimap2 alignment is available
* --input-aln: Input SAM/BAM file if already aligned by minimap2
* --threads: Number of threads
* --mm-args: Arguments for minimap2 (default: -ax map-ont)

example usage for PacBio HiFi/CCS genomic reads:
```bash
mmbwa ref_genome --input-aln /path/to/file.bam --output /path/to/outdir --mm-args "-ax map-hifi"
```

* --bwa_args: Arguments for bwa-mem (default: None)
* --output: Path to output directory
* --regions: Option to exclude specified genomic regions from bwa-mem alignment. E.g. mapping in centromeric regions 
considerably slows down the alignment
* --regions_bed: BED file with genomic intervals
* --regions_overlap: Min overlap threshold with regions. If the primary alignment overlaps a region but is below this 
threshold that primary will not be excluded from bwa-mem alignment
* --keep_temp: Keep temporary files
  * minimap2.sam: SAM file generated after minimap2 alignment
  * filtered.bam: BAM file containing reads that are below the soft-clip threshold
  * softclipped.bam: BAM file containing reads that are above the soft-clip threshold - these are fed to bwa-mem
* --threshold: If the fraction of the soft-clip length/read length is above this threshold, the read is re-aligned with 
bwa-mem, default = 0.1
* --sort: Sort the final output BAM
* --index: Index the final output BAM
* --unmapped: Re-align reads unmapped by minimap2 (not tested yet)

The final BAM output is created by merging reads that were below the soft-clipping threshold with the reads re-aligned by bwa-mem.


## Benchmark