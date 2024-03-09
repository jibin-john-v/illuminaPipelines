import os

##https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)
##https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/030-tutorial-downstream/020-taxonomy.html


#Filter out rare ASVs
def qiime_filter_rare(feature_table,min_frequency):
    '''Filter out rare ASVs'''
    os.system(f'''qiime feature-table filter-features \
             --i-table {feature_table} \
             --p-min-frequency {min_frequency} \
             --p-min-samples 1 \
             --o-filtered-table feature-table-1.qza''')

def qiime_filter_contaminant(feature_table,taxonomy):
    '''Filter out contaminant and unclassified ASVs'''
    os.system(f'''qiime taxa filter-table \
         --i-table {feature_table} \
         --i-taxonomy {taxonomy} \
         --p-mode contains \
         --p-include p__ \
         --p-exclude 'p__;,Chloroplast,Mitochondria' \
         --o-filtered-table feature-table-2.qza''')

def qiime_filter_low_sequence_counts(feature_table,min_depth):
    '''Filtering samples with low sequence counts'''
    os.system(f'''qiime feature-table filter-samples \
                --i-table {feature_table}\
                --p-min-frequency {min_depth} \
                --o-filtered-table filtered-table-qc-passed.qza''')

def qiime_filter_asv_sequences(feature_sequence,feature_table):
    '''remove ASV sequences that are no longer represented in our table'''
    os.system(f'''qiime feature-table filter-seqs \
                 --i-data {feature_sequence} \
                 --i-table {feature_table} \
                 --o-filtered-data filtered-sequences-qc-passed.qza''')


def qiime_filter_rare_fraction(feature_table,max_depth):
    '''rarefaction'''
    os.system(f'''qiime diversity alpha-rarefaction \
   --i-table {feature_table} \
   --p-max-depth {max_depth} \
   --p-steps 20 \
   --p-metrics 'observed_features' \
   --o-visualization rarefaction_curves_test.qzv''')




'''
qiime feature-table filter-features \
   --i-table deblur_output/table.qza \
   --p-min-frequency X \
   --p-min-samples 1 \
   --o-filtered-table deblur_output/deblur_table_filt.qza


#Filter out contaminant and unclassified ASVs
qiime taxa filter-table \
   --i-table deblur_output/deblur_table_filt.qza \
   --i-taxonomy taxa/classification.qza \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table deblur_output/deblur_table_filt_contam.qza


qiime taxa filter-table \
   --i-table deblur_output/deblur_table_filt.qza \
   --i-taxonomy taxa/classification.qza \
   --p-include p__ \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table deblur_output/deblur_table_filt_contam.qza


#Exclude low-depth samples
qiime feature-table summarize \
   --i-table deblur_output/deblur_table_filt_contam.qza \
   --o-visualization deblur_output/deblur_table_filt_contam_summary.qzv

qiime diversity alpha-rarefaction \
   --i-table deblur_output/deblur_table_filt_contam.qza \
   --p-max-depth X \
   --p-steps 20 \
   --p-metrics 'observed_features' \
   --o-visualization rarefaction_curves_test.qzv

qiime feature-table filter-samples \
   --i-table deblur_output/deblur_table_filt_contam.qza \
   --p-min-frequency SET_CUTOFF \
   --o-filtered-table deblur_output/deblur_table_final.qza


#Subset and summarize filtered table

qiime feature-table filter-seqs \
   --i-data deblur_output/representative_sequences.qza \
   --i-table deblur_output/deblur_table_final.qza \
   --o-filtered-data deblur_output/rep_seqs_final.qza

qiime feature-table summarize \
   --i-table deblur_output/deblur_table_final.qza \
   --o-visualization deblur_output/deblur_table_final_summary.qzv
'''