# pylint: disable=missing-module-docstring
import os

##Reference https://github.com/beiko-lab/CBW2021_Module2_16S_Analysis/wiki/MIC-Module-2-tutorial
##                      https://github.com/qiime2/docs/blob/master/source/tutorials/read-joining.rst#id1 
##                      https://docs.qiime2.org/2023.9/tutorials/fmt/
##                      https://github.com/LangilleLab/microbiome_helper/wiki/QIIME2-DADA2-Quick-Reference
##                      https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/020-tutorial-upstream/040-denoising.html

## Database reference : https://github.com/qiime2/docs/blob/master/source/data-resources.rst

#Clustering 16S sequences

def dada2_denoise_paired(fastq_imported,n_threads,trunc_len_f,trunc_len_r):
    '''Performing sequence quality control (i.e., denoising) for paired end data'''
    os.system(f'''
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {fastq_imported} \
            --p-trunc-len-f {trunc_len_f} \
            --p-trunc-len-r {trunc_len_r} \
            --p-trim-left-f 0 \
            --p-trim-left-r 0 \
            --p-max-ee-f 2 \
            --p-max-ee-r 3 \
            --o-table feature-table-0.qza \
            --p-n-threads {n_threads} \
            --o-representative-sequences asv-sequences-0.qza \
            --o-denoising-stats dada2-stats.qza ''')

def dada2_denoise_single(fastq_imported,n_threads,trunc_len):
    '''Performing sequence quality control (i.e., denoising) for single end data'''
    os.system(f'''
        qiime dada2 denoise-single \
            --i-demultiplexed-seqs {fastq_imported} \
            --p-trunc-len {trunc_len} \
            --p-trim-left 0 \
            --p-max-ee 2 \
            --o-table feature-table-0.qza \
            --p-n-threads {n_threads} \
            --o-representative-sequences asv-sequences-0.qza \
            --o-denoising-stats dada2-stats.qza ''')


def dada2_stat_visualize(dada2_stats):
    '''generate a viewable summary'''
    os.system(f'''qiime metadata tabulate \
                --m-input-file {dada2_stats} \
                --o-visualization dada2-stats-summ.qzv''')


def qiime2_denoisesummary(unfiltered_table_qza,repseqs,metadata):
    '''Generating summaries of the feature table and feature data'''
    os.system(f'''qiime feature-table summarize \
                    --i-table {unfiltered_table_qza} \
                    --m-sample-metadata-file {metadata} \
                    --o-visualization feature-table-0-summary.qzv''')
    os.system(f'''qiime feature-table tabulate-seqs \
                    --i-data {repseqs} \
                    --o-visualization asv-sequences-0-summary.qzv''')