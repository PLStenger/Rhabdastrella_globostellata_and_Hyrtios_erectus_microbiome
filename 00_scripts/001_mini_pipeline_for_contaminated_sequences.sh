#!/usr/bin/env bash

WORKING_DIRECTORY=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/99_contamination
OUTPUT=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/99_contamination

# You need a "manifest" for explaining and importing your data :
# In the fastq manifest formats, a manifest file maps sample identifiers to fastq.gz or fastq absolute filepaths that contain sequence and quality data for the sample, and indicates the direction of the reads in each fastq.gz / fastq absolute filepath. The manifest file will generally be created by you, and it is designed to be a simple format that doesn’t put restrictions on the naming of the demultiplexed fastq.gz / fastq files, since there is no broadly used naming convention for these files. There are no restrictions on the name of the manifest file.
# See https://docs.qiime2.org/2018.8/tutorials/importing/

###############################################################
### For importing your data in a Qiime2 format
###############################################################

cd $WORKING_DIRECTORY

eval "$(conda shell.bash hook)"
conda activate qiime2-2021.4


qiime tools import \
  --input-path contamination_seq.fasta \
  --output-path contamination_seq.qza \
  --type 'FeatureData[Sequence]'
