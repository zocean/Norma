#!/bin/bash

module load Bowtie2/2.1.0-IGB-gcc-4.9.4
module load SAMtools/1.5-IGB-gcc-4.9.4
module load Python/2.7.12-IGB-gcc-4.9.4

# python biocluster_pipeline.py --conf TSA-seq_conf.txt --dry_run
python biocluster_pipeline.py --conf TSA-seq_conf.txt
