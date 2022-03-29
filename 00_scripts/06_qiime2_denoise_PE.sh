#!/usr/bin/env bash

WORKING_DIRECTORY=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/05_QIIME2
OUTPUT=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/05_QIIME2/visual

METADATA=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/98_database_files/sample-metadata.tsv
# negative control sample :
#NEG_CONTROL=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/98_database_files/Negative_control_Sample_RepSeq_V4.qza
#NEG_CONTROL=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/98_database_files/SILVA-138-SSURef-Full-Seqs.qza
NEG_CONTROL=/Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Dysidea_microbiome/Rhabdastrella_globostellata_and_Dysidea_microbiome/99_contamination/contamination_seq.qza

# https://chmi-sops.github.io/mydoc_qiime2.html

# https://docs.qiime2.org/2021.2/plugins/available/dada2/denoise-single/
# Aim: denoises single-end sequences, dereplicates them, and filters
# chimeras and singletons sequences
# Use: qiime dada2 denoise-single [OPTIONS]

# DADA2 method

cd $WORKING_DIRECTORY

eval "$(conda shell.bash hook)"
conda activate qiime2-2021.4

# dada2_denoise :
#################

# Aim: denoises paired-end sequences, dereplicates them, and filters
#      chimeras and singletons sequences

# https://github.com/benjjneb/dada2/issues/477

qiime dada2 denoise-paired --i-demultiplexed-seqs core/demux.qza \
--o-table core/Table.qza  \
--o-representative-sequences core/RepSeq.qza \
--o-denoising-stats core/Stats.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 4  

# sequence_contamination_filter :
#################################

# Aim: aligns feature sequences to a set of reference sequences
#      to identify sequences that hit/miss the reference
#      Use: qiime quality-control exclude-seqs [OPTIONS]

# Here --i-reference-sequences correspond to the negative control sample (if you don't have any, like here, take another one from an old project, the one here is from the same sequencing line (but not same project))

qiime quality-control exclude-seqs --i-query-sequences core/RepSeq.qza \
      					     --i-reference-sequences $NEG_CONTROL \
      					     --p-method vsearch \
      					     --p-threads 6 \
      					     --p-perc-identity 1.00 \
      					     --p-perc-query-aligned 1.00 \
      					     --o-sequence-hits core/HitNegCtrl.qza \
      					     --o-sequence-misses core/NegRepSeq.qza
                                         

# table_contamination_filter :
##############################

# Aim: filter features from table based on frequency and/or metadata
#      Use: qiime feature-table filter-features [OPTIONS]

qiime feature-table filter-features --i-table core/Table.qza \
     					      --m-metadata-file core/HitNegCtrl.qza \
     					      --o-filtered-table core/NegTable.qza \
     					      --p-exclude-ids

# table_contingency_filter :
############################

# Aim: filter features that show up in only one samples, based on
#      the suspicion that these may not represent real biological diversity
#      but rather PCR or sequencing errors (such as PCR chimeras)
#      Use: qiime feature-table filter-features [OPTIONS]

# contingency:
    # min_obs: 2  # Remove features that are present in only a single sample !
    # min_freq: 0 # Remove features with a total abundance (summed across all samples) of less than 0 !


qiime feature-table filter-features  --i-table core/NegTable.qza \
        					       --p-min-samples 2 \
        					       --p-min-frequency 0 \
        					       --o-filtered-table core/ConTable.qza


# sequence_contingency_filter :
###############################

# Aim: Filter features from sequence based on table and/or metadata
       # Use: qiime feature-table filter-seqs [OPTIONS]

qiime feature-table filter-seqs --i-data core/NegRepSeq.qza \
      					  --i-table core/ConTable.qza \
      					  --o-filtered-data core/ConRepSeq.qza


# sequence_summarize :
######################

# Aim: Generate tabular view of feature identifier to sequence mapping
       # Use: qiime feature-table tabulate-seqs [OPTIONS]

qiime feature-table summarize --i-table core/Table.qza --m-sample-metadata-file $METADATA --o-visualization visual/Table.qzv
qiime feature-table summarize --i-table core/ConTable.qza --m-sample-metadata-file $METADATA --o-visualization visual/ConTable.qzv
qiime feature-table summarize --i-table core/NegTable.qza --m-sample-metadata-file $METADATA --o-visualization visual/NegTable.qzv
qiime feature-table tabulate-seqs --i-data core/NegRepSeq.qza --o-visualization visual/NegRepSeq.qzv
qiime feature-table tabulate-seqs --i-data core/RepSeq.qza --o-visualization visual/RepSeq.qzv
qiime feature-table tabulate-seqs --i-data core/HitNegCtrl.qza --o-visualization visual/HitNegCtrl.qzv

mkdir -p export/core
mkdir -p export/visual

qiime tools export --input-path core/Table.qza --output-path export/core/Table
qiime tools export --input-path core/ConTable.qza --output-path export/core/ConTable
qiime tools export --input-path core/NegTable.qza --output-path export/core/NegTable
qiime tools export --input-path core/NegRepSeq.qza --output-path export/core/NegRepSeq
qiime tools export --input-path core/RepSeq.qza --output-path export/core/RepSeq
qiime tools export --input-path core/HitNegCtrl.qza --output-path export/core/HitNegCtrl
qiime tools export --input-path core/ConRepSeq.qza --output-path export/core/ConRepSeq
qiime tools export --input-path core/Stats.qza  --output-path export/core/Stats

qiime tools export --input-path visual/NegTable.qzv --output-path export/visual/NegTable
qiime tools export --input-path visual/ConTable.qzv --output-path export/visual/ConTable
qiime tools export --input-path visual/Table.qzv --output-path export/visual/Table
qiime tools export --input-path visual/HitNegCtrl.qzv --output-path export/visual/HitNegCtrl
qiime tools export --input-path visual/RepSeq.qzv --output-path export/visual/RepSeq
qiime tools export --input-path visual/NegRepSeq.qzv --output-path export/visual/NegRepSeq



#########!/usr/bin/env bash
########
########WORKING_DIRECTORY=/scratch_vol1/fungi/Rhabdastrella_globostellata_microbiome/05_QIIME2
########OUTPUT=/scratch_vol1/fungi/Rhabdastrella_globostellata_microbiome/05_QIIME2/visual
########
########METADATA=/scratch_vol1/fungi/Rhabdastrella_globostellata_microbiome/98_database_files/sample-metadata.tsv
######### negative control sample :
########NEG_CONTROL=/scratch_vol1/fungi/Rhabdastrella_globostellata_microbiome/98_database_files/Negative_control_Sample_RepSeq_V4.qza
########
######### https://chmi-sops.github.io/mydoc_qiime2.html
########
######### https://docs.qiime2.org/2021.2/plugins/available/dada2/denoise-single/
######### Aim: denoises single-end sequences, dereplicates them, and filters
######### chimeras and singletons sequences
######### Use: qiime dada2 denoise-single [OPTIONS]
########
######### DADA2 method
########
########cd $WORKING_DIRECTORY
########
########eval "$(conda shell.bash hook)"
########conda activate qiime2-2021.4
########
######### dada2_denoise :
#########################
########
######### Aim: denoises paired-end sequences, dereplicates them, and filters
#########      chimeras and singletons sequences
########
######### https://github.com/benjjneb/dada2/issues/477
########
######### https://docs.qiime2.org/2020.8/plugins/available/dada2/denoise-paired/
######### p-trunc : Position at which forward read sequences should be truncated due to decrease in quality. This truncates the 3' end of the of the input sequences, which will be the bases that were sequenced in the last cycles. Reads that are shorter than this value will be discarded. After this parameter is applied there must still be at least a 12 nucleotide overlap between the forward and reverse reads. If 0 is provided, no truncation or length filtering will be performed
######### --p-trim-left-f : Position at which forward read sequences should be trimmed due to low quality. This trims the 5' end of the input sequences, which will be the bases that were sequenced in the first cycles.
#########--p-trim-left-r : Position at which reverse read sequences should be trimmed due to low quality. This trims the 5' end of the input sequences, which will be the bases that were sequenced in the first cycles.
########
#########qiime dada2 denoise-paired --i-demultiplexed-seqs core/demux.qza \
#########--o-table core/Table.qza  \
#########--o-representative-sequences core/RepSeq.qza \
#########--o-denoising-stats core/Stats.qza \
#########--p-trim-left-f 0 \
#########--p-trim-left-r 0 \
#########--p-trunc-len-f 0 \
#########--p-trunc-len-r 0 \
#########--p-n-threads 4                        
########                           
######### For seeing the SampleData file (convert from .qza to .qzv)
########qiime metadata tabulate \
########  --m-input-file core/Stats.qza \
########  --o-visualization visual/Stats.qzv
########
########## For seeing the Table file (convert from .qza to .qzv)
#########qiime feature-table summarize --i-table core/Table.qza --o-visualization visual/Table.qzv
#########
########## For seeing the RedSeq file (convert from .qza to .qzv)
#########qiime feature-table tabulate-seqs --i-data core/RepSeq.qza --o-visualization visual/RepSeq.qzv
########
########
######### sequence_contamination_filter :
#########################################
########
######### Aim: aligns feature sequences to a set of reference sequences
#########      to identify sequences that hit/miss the reference
#########      Use: qiime quality-control exclude-seqs [OPTIONS]
########
######### 0.97 is the default (see https://docs.qiime2.org/2020.2/plugins/available/quality-control/exclude-seqs/)
######### Whereas N. Fernandez put 1.00 but didn't work forward (" All features were filtered out of the data." error to the step "qiime feature-table filter-seqs")
########
######### Here --i-reference-sequences correspond to the negative control sample (if you don't have any, like here, take another one from an old project, the one here is from the same sequencing line (but not same project))
########
#########    --i-reference-sequences /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/05_QIIME2/old_analysis/Negative_control/V4/Negative_control_Sample_RepSeq_V4.qza \
########
########
########qiime quality-control exclude-seqs --i-query-sequences core/RepSeq.qza \
########      					     --i-reference-sequences $NEG_CONTROL \
########      					     --p-method vsearch \
########      					     --p-threads 4 \
########      					     --p-perc-identity 1.00 \
########      					     --p-perc-query-aligned 1.00 \
########      					     --o-sequence-hits core/HitNegCtrl.qza \
########      					     --o-sequence-misses core/NegRepSeq.qza
########               
########               
######### For seeing the NegRepSeq file (convert from .qza to .qzv)
########qiime feature-table tabulate-seqs --i-data core/NegRepSeq.qza --o-visualization visual/NegRepSeq.qzv
########        
######### For seeing the HitNegCtrl file (convert from .qza to .qzv)        
########qiime feature-table tabulate-seqs --i-data core/HitNegCtrl.qza --o-visualization visual/HitNegCtrl.qzv
########        
########                     
########
######### table_contamination_filter :
######################################
########
######### Aim: filter features from table based on frequency and/or metadata
#########      Use: qiime feature-table filter-features [OPTIONS]
########
######### --p-exclude-ids: --p-no-exclude-ids If true, the samples selected by `metadata` or `where` parameters will be excluded from the filtered table instead of being retained. [default: False]:
########
########qiime feature-table filter-features --i-table core/Table.qza \
########     					      --m-metadata-file core/HitNegCtrl.qza \
########     					      --o-filtered-table core/NegTable.qza \
########     					      --p-exclude-ids
########
######### For seeing the HitNegCtrl file (convert from .qza to .qzv)
########qiime feature-table tabulate-seqs --i-data core/HitNegCtrl.qza --o-visualization visual/HitNegCtrl.qzv
########
######### For seeing the NegTable file (convert from .qza to .qzv)
########qiime feature-table summarize --i-table core/NegTable.qza --o-visualization visual/NegTable.qzv
########
########
########
########
######### table_contingency_filter :
####################################
########
######### Aim: filter features that show up in only one samples, based on
#########      the suspicion that these may not represent real biological diversity
#########      but rather PCR or sequencing errors (such as PCR chimeras)
#########      Use: qiime feature-table filter-features [OPTIONS]
########
######### contingency:
########    # min_obs: 2  # Remove features that are present in only a single sample !
########    # min_freq: 0 # Remove features with a total abundance (summed across all samples) of less than 0 !
########
########
########qiime feature-table filter-features  --i-table core/NegTable.qza \
########        					       --p-min-samples 2 \
########        					       --p-min-frequency 0 \
########        					       --o-filtered-table core/ConTable.qza
########
######### For seeing the ConTable file (convert from .qza to .qzv)
########qiime feature-table summarize --i-table core/ConTable.qza --o-visualization visual/ConTable.qzv
########
########
######### sequence_contingency_filter :
#######################################
########
######### Aim: Filter features from sequence based on table and/or metadata
########       # Use: qiime feature-table filter-seqs [OPTIONS]
########
########qiime feature-table filter-seqs --i-data core/NegRepSeq.qza \
########      					  --i-table core/ConTable.qza \
########      					  --o-filtered-data core/ConRepSeq.qza
########
######### For seeing the ConRepSeq file (convert from .qza to .qzv)
########qiime feature-table tabulate-seqs --i-data core/ConRepSeq.qza --o-visualization visual/ConRepSeq.qzv
########
######### sequence_summarize :
##############################
########
######### Aim: Generate tabular view of feature identifier to sequence mapping
########       # Use: qiime feature-table tabulate-seqs [OPTIONS]
########
########qiime feature-table summarize --i-table core/Table.qza --m-sample-metadata-file $METADATA --o-visualization visual/Table.qzv
########qiime feature-table summarize --i-table core/ConTable.qza --m-sample-metadata-file $METADATA --o-visualization visual/ConTable.qzv
########qiime feature-table summarize --i-table core/NegTable.qza --m-sample-metadata-file $METADATA --o-visualization visual/NegTable.qzv
########qiime feature-table tabulate-seqs --i-data core/NegRepSeq.qza --o-visualization visual/NegRepSeq.qzv
########qiime feature-table tabulate-seqs --i-data core/RepSeq.qza --o-visualization visual/RepSeq.qzv
########qiime feature-table tabulate-seqs --i-data core/HitNegCtrl.qza --o-visualization visual/HitNegCtrl.qzv
########
########mkdir -p export/core
########mkdir -p export/visual
########
########qiime tools export --input-path core/Table.qza --output-path export/core/Table
########qiime tools export --input-path core/ConTable.qza --output-path export/core/ConTable
########qiime tools export --input-path core/NegTable.qza --output-path export/core/NegTable
########qiime tools export --input-path core/NegRepSeq.qza --output-path export/core/NegRepSeq
########qiime tools export --input-path core/RepSeq.qza --output-path export/core/RepSeq
########qiime tools export --input-path core/HitNegCtrl.qza --output-path export/core/HitNegCtrl
########qiime tools export --input-path core/ConRepSeq.qza --output-path export/core/ConRepSeq
########qiime tools export --input-path core/Stats.qza  --output-path export/core/Stats
########
########qiime tools export --input-path visual/NegTable.qzv --output-path export/visual/NegTable
########qiime tools export --input-path visual/ConTable.qzv --output-path export/visual/ConTable
########qiime tools export --input-path visual/Table.qzv --output-path export/visual/Table
########qiime tools export --input-path visual/HitNegCtrl.qzv --output-path export/visual/HitNegCtrl
########qiime tools export --input-path visual/RepSeq.qzv --output-path export/visual/RepSeq
########qiime tools export --input-path visual/NegRepSeq.qzv --output-path export/visual/NegRepSeq
