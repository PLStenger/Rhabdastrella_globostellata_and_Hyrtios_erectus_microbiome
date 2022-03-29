#!/usr/bin/env bash

WORKING_DIRECTORY=/scratch_vol1/fungi/Rhabdastrella_globostellata_microbiome/03_cleaned_data
OUTPUT=/scratch_vol1/fungi/Rhabdastrella_globostellata_microbiome/03_cleaned_data

# Here the WORKING_DIRECTORY and the OUTPUT is the same
# I only try to deal with the V1V3 - V3V4 combinaison

cd $WORKING_DIRECTORY

# By cat technic (for "normal" sequences and for "+" sequences)

cat BRG1V1_R1.paired.fastq.gz BRG1V3_R1.paired.fastq.gz > BRG1_R1_cat.fastq.gz
cat BRG1V3_R2.paired.fastq.gz BRG1V4_R2.paired.fastq.gz > BRG1_R2_cat.fastq.gz

cat BRG2V1_R1.paired.fastq.gz BRG2V3_R1.paired.fastq.gz > BRG2_R1_cat.fastq.gz
cat BRG2V3_R2.paired.fastq.gz BRG2V4_R2.paired.fastq.gz > BRG2_R2_cat.fastq.gz

cat BRG3V1_R1.paired.fastq.gz BRG3V3_R1.paired.fastq.gz > BRG3_R1_cat.fastq.gz
cat BRG3V3_R2.paired.fastq.gz BRG3V4_R2.paired.fastq.gz > BRG3_R2_cat.fastq.gz

#cat BRG1V1_R1.paired.fastq.gz BRG1V3+_R1.paired.fastq.gz > BRG1+_R1_cat.fastq.gz
#cat BRG1V3_R2.paired.fastq.gz BRG1V4+_R2.paired.fastq.gz > BRG1+_R2_cat.fastq.gz

#cat BRG2V1_R1.paired.fastq.gz BRG2V3+_R1.paired.fastq.gz > BRG2+_R1_cat.fastq.gz
#cat BRG2V3_R2.paired.fastq.gz BRG2V4+_R2.paired.fastq.gz > BRG2+_R2_cat.fastq.gz

#cat BRG3V1_R1.paired.fastq.gz BRG3V3+_R1.paired.fastq.gz > BRG3+_R1_cat.fastq.gz
#cat BRG3V3_R2.paired.fastq.gz BRG3V4+_R2.paired.fastq.gz > BRG3+_R2_cat.fastq.gz

