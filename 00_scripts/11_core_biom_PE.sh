#!/usr/bin/env bash

WORKING_DIRECTORY=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/05_QIIME2
OUTPUT=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/05_QIIME2/visual

DATABASE=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/98_database_files

###############################################################
### For Bacteria
###############################################################

cd $WORKING_DIRECTORY

eval "$(conda shell.bash hook)"
conda activate qiime2-2021.4

# Make the directory (mkdir) only if not existe already(-p)
mkdir -p subtables
mkdir -p export/subtables

# Aim: Filter sample from table based on a feature table or metadata

qiime feature-table filter-samples \
        --i-table core/RarTable.qza \
        --m-metadata-file $DATABASE/sample-metadata.tsv \
        --p-where "[#SampleID] IN ('Bourake_Rhabdastrella_globostellata_2019_1_V1V3', 'Bourake_Rhabdastrella_globostellata_2019_2_V1V3', 'Bourake_Rhabdastrella_globostellata_2019_3_V1V3', 'Bourake_Dysidea_2019_1_V1V3', 'Bourake_Dysidea_2019_2_V1V3', 'Bourake_Dysidea_2019_3_V1V3', 'Bourake_SeaWater_downstream_0_2micro_V1V3', 'Bourake_SeaWater_downstream_3micro_V1V3', 'Bourake_SeaWater_upstream_0_2micro_V1V3', 'Bourake_SeaWater_upstream_3micro_V1V3')"  \
        --o-filtered-table subtables/RarTable-all.qza
 
# Aim: Identify "core" features, which are features observed,
     # in a user-defined fraction of the samples
        
qiime feature-table core-features \
        --i-table subtables/RarTable-all.qza \
        --p-min-fraction 0.1 \
        --p-max-fraction 1.0 \
        --p-steps 10 \
        --o-visualization visual/CoreBiom-all.qzv  
        
qiime tools export --input-path subtables/RarTable-all.qza --output-path export/subtables/RarTable-all    
qiime tools export --input-path visual/CoreBiom-all.qzv --output-path export/visual/CoreBiom-all
biom convert -i export/subtables/RarTable-all/feature-table.biom -o export/subtables/RarTable-all/table-from-biom.tsv --to-tsv
sed '1d ; s/\#OTU ID/ASV_ID/' export/subtables/RarTable-all/table-from-biom.tsv > export/subtables/RarTable-all/ASV.tsv
