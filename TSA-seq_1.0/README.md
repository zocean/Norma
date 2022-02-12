## TSA-seq 1.0 processing code

## Introduction

This repo hosts the pipeline for processing and analyzing genome-wide proximity mapping data generated by TSA-seq technique. Users can apply script to finish the following task:

- Normalized and processed the TSA-seq data using a sliding window approach.
- Partition the genome into SON-enriched and SON-depleted region.
- Analyze the correlations of SON TSA-seq with other Genomic features.
- Generated data that can be loaded in UCSC Genome Browser tracks for visualization.

The codes have been tested under:

- Python (tested in Python 2.7)
- argparse
- numpy
- scipy
- bx-python
- R (tested in have any requirement of read mapping procedure

## Reads mapping

The pipeline does not have any requirement of read mapping procedure. It begins with bam file. For read mapping, any well-known mapping tools such as Bowtie, BWA and SOAP2 can be used. We do suggest to remove PCR duplicates after the mapping process, and bam file should be sorted and indexed. For example, in the paper we used the following code to map reads to human genome and generate bam files

```shell
# Since K562 is derived from a female, we remove chromosome Y from human reference genome hg19 and named it hg19F
bowtie2 -p 8 -x hg19F -U SON_TSA-seq_pulldown.fastq -S SON_TSA-seq_pulldown.sam

# Convert Sam to Bam
samtools view -bS SON_TSA-seq_pulldown.sam > SON_TSA-seq_pulldown.bam

# You can also combine the above two steps if you do not want to save sam file on the disk using
# bowtie2 -p 8 -x hg19F -U SON_TSA-seq_pulldown.fastq | samtools view -bS - > SON_TSA-seq_pulldown.bam

# Sort bam file
samtools sort SON_TSA-seq_pulldown.bam SON_TSA-seq_pulldown_sort

# remove pcr duplicate from bam file
samtools rmdup SON_TSA-seq_pulldown_sort SON_TSA-seq_pulldown_rmdup.bam

# Index bam file
samtools index SON_TSA-seq_pulldown_rmdup.bam
```
## Normalization

We then normalized and processed the TSA-seq data using a sliding window approaches. The basic idea is to calculate the fold change ratio in each sliding window between pulldown sample and input sample normalized by the total number of mapped reads. Details of method equation can be found in the paper. In the paper, we used window size of 20kb with sliding window step 100bp:

```shell
# Normalize SON TSA-seq pulldown with matched input
python Norma_TSA-seq_1.0.py -N 20000 -r 100 -l 100 -o LMNB_on_input_20k -e SON_TSA-seq_pulldown_rmdup.bam -c SON_TSA-seq_pulldown_input_rmdump.bam

# the meaning the of parameters are shown below:
# -N => window size in base pair (bp)
# -r => sliding window step (bp)
# -l => read length (bp)
# -e => bam file for pulldown sample
# -c => bam file for matched input sample
# -o => output file prefix
```
The normalization step will generate three output files. One for the normalized TSA-seq score. The other two wig files show the signal profiles for pulldown sample and input sample under the resolution defined by -N and -r.

## About

The code for processing TSA-seq data is authored and maintained by [Yang Zhang](mailto:yangz6@andrew.cmu.edu)

## Cite

To cite TSA-seq 1.0 paper, use the following information.
```
@article{chen2018mapping,
  title={Mapping 3D genome organization relative to nuclear compartments using TSA-Seq as a cytological ruler},
  author={Chen, Yu and Zhang, Yang and Wang, Yuchuan and Zhang, Liguo and Brinkman, Eva K and Adam, Stephen A and Goldman, Robert and Van Steensel, Bas and Ma, Jian and Belmont, Andrew S},
  journal={Journal of Cell Biology},
  volume={217},
  number={11},
  pages={4025--4048},
  year={2018},
  publisher={Rockefeller University Press}
}
```