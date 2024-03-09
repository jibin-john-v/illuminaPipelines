import sys,os,glob,subprocess
from pathlib import Path
import pandas as pd
import concurrent.futures


sys.path.insert(0, '../../master_codes/') 
from fastq_qc.fastq_qc import falco_run,fastp_run
from utility.utility import fastq_files_location_df,fastq_files_location_df
from postalignment.markduplicate_bqsr  import sambamba_unique_mdup_unmap_run
from alignment.alignment import bowtie2_aligner_run
from postalignment.postalignment  import bamindex_run
from chipseq.chipseq import macs2_narrowpeak_call,macs2_broadpeak_call,remove_blacklist_region,bigWig_filecreation,bigWig_comparison_filecreation


'''References 
            https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/05_filtering_BAM_files.md
            https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/schedule/2-day.md
            https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/chipseq_tutorials.html
            https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html #For QC
            https://bioinformaticsdotca.github.io/EPI_2021
            https://nbisweden.github.io/workshop-archive/workshop-ChIP-seq/2018-11-07/labs/lab-processing.html
            https://www.sciencedirect.com/science/article/pii/S1046202320300591
            https://www.frontiersin.org/articles/10.3389/fcell.2023.1167111/full
''''


genome_index="/Users/JJOHN41/Desktop/projects/test_data/reference_files/GRCm39.primary_assembly.genome"
threads=10
binSize=20
fastq_path = "/Users/JJOHN41/Desktop/projects/test_data/raw_fastq_files/"
effective_genome_size=2913022398
meta_df="/Users/JJOHN41/Desktop/projects/test_data/raw_fastq_files/single_end_fastq_files_locations.csv"
blacklist_bed="/Users/JJOHN41/Desktop/projects/test_data/reference_files/mm10-blacklist.v2.bed"
paired_end, single_end = fastq_files_location_df(fastq_path)


#run_fastqc_parallel("falco_run",single_end, n_parallel=5,n_threads=3,out_name="output",qc_cutoff=20)
for index in range(single_end.shape[0]):
    falco_run(single_end.loc[index,"SampleName"],[single_end.loc[index,"R1"]],"falco")


#run_fastqc_parallel("fastp_run",single_end, n_parallel=5,n_threads=3,out_name="output",qc_cutoff=20)
for index in range(single_end.shape[0]):
    fastp_run(single_end.loc[index,"SampleName"],[single_end.loc[index,"R1"]],5,20)

#Alignment using bowtie2
for index in range(single_end.shape[0]):
    s_name=single_end.loc[index,"SampleName"]
    bowtie2_aligner_run(s_name,[f'Results/cleaned_fastQ_files/{s_name}_trimmed.fq.gz'],genome_index,threads)

#potalignment qc using sambamba
for index in range(single_end.shape[0]):
    s_name=single_end.loc[index,"SampleName"]
    sambamba_unique_mdup_unmap_run(s_name, f'Results/Bam_Files/{s_name}.bam', threads)
    bamindex_run(f'Results/Bam_Files/{s_name}_final.bam')


meta_df=pd.read_csv(meta_df) #meta_df.columns =['sample', 'R1', 'R2', 'group', 'assay', 'chip_sample']
meta_df_pivot = meta_df.pivot(index="chip_sample", columns="assay", values=["sample", "group"])
meta_df_pivot = meta_df_pivot.reset_index()
meta_df_pivot.columns=meta_df_pivot.columns.droplevel(0)
meta_df_pivot=meta_df_pivot.iloc[:,1:4]
meta_df_pivot.columns=["chip","input","group"]

#peak calling
for index in range(meta_df_pivot.shape[0]):
    chip_name=meta_df_pivot.loc[index,"chip"]
    input_name=meta_df_pivot.loc[index,"input"]
    macs2_narrowpeak_call(chip_name, 
                               f"Results/Bam_Files/{chip_name}_final.bam", 
                               f"Results/Bam_Files/{input_name}_final.bam",
                               effective_genome_size)
    macs2_broadpeak_call(chip_name, 
                               f"Results/Bam_Files/{chip_name}_final.bam", 
                               f"Results/Bam_Files/{input_name}_final.bam",
                               effective_genome_size)

##remove blacklisted region
for index in range(meta_df_pivot.shape[0]):
    chip_name=meta_df_pivot.loc[index,"chip"]
    remove_blacklist_region(blacklist_bed, 
                            f'Results/broad_peak_calling/{chip_name}_peaks.broadPeak',
                            'Results/broad_peak_calling/',
                            f'{chip_name}_peaks_filtered.broadPeak')
    remove_blacklist_region(blacklist_bed, 
                            f'Results/narrow_peak_calling/{chip_name}_peaks.narrowPeak',
                            'Results/narrow_peak_calling/',
                            f'{chip_name}_peaks_filtered.narrowPeak')


input_bamfiles=[]
chip_bamfiles=[]
input_names=[]
chip_names=[]

for index in range(meta_df_pivot.shape[0]):
    chip_name=meta_df_pivot.loc[index,"chip"]
    input_name=meta_df_pivot.loc[index,"input"]
    input_bamfiles.append(f"Results/Bam_Files/{input_name}_final.bam")
    chip_bamfiles.append(f"Results/Bam_Files/{chip_name}_final.bam")
    input_names.append(f'{input_name}_i')
    chip_names.append(f'{chip_name}_c')

bam_files=" ".join(chip_bamfiles+input_bamfiles)
s_names=" ".join(chip_names+input_names)


os.system(f'''multiBamSummary bins --bamfiles {bam_files} \
        --outFileName multiBamArray_dT201_preproc_bam_chr12.npz \
        --binSize=5000 \
        --extendReads=110 --labels {s_names}''')


os.system(f'''plotCorrelation \
                        --corData multiBamArray_dT201_preproc_bam_chr12.npz \
                        --plotFile REST_bam_correlation_bin.pdf \
                        --outFileCorMatrix corr_matrix_bin.txt \
                        --whatToPlot heatmap \
                        --corMethod spearman''')



#bigWig file, creation
for index in range(single_end.shape[0]):
    s_name=single_end.loc[index,"SampleName"]
    bigWig_filecreation(s_name,f'Results/Bam_Files/{s_name}_final.bam',binSize,"Results/bigWig-file")

for index in range(meta_df_pivot.shape[0]):
    chip_name=meta_df_pivot.loc[index,"chip"]
    input_name=meta_df_pivot.loc[index,"input"]
    chip_bam=meta_df_pivot.loc[index,"chip"]
    bigWig_comparison_filecreation(chip_name,
                                   input_name,
                                   f"Results/Bam_Files/{chip_name}_final.bam",
                                   f"Results/Bam_Files/{input_name}_final.bam",
                                   binSize,
                                   "Results/bigWig-file")

##create meta data file for ChIPQC
chipqc_df=meta_df_pivot.copy() ##example https://raw.githubusercontent.com/hbctraining/Intro-to-ChIPseq/master/samplesheet_chr12.csv
chipqc_df=chipqc_df.rename(columns={"chip":"SampleID","input":"ControlID","group":"Factor"})
num_factors = chipqc_df['Factor'].nunique()
chipqc_df['Replicate']=chipqc_df.groupby('Factor').cumcount() // (num_factors // 2) + 1
chipqc_df['bamReads']="Results/Bam_Files/"+chipqc_df['SampleID']+'_final.bam'
chipqc_df['bamControl']="Results/Bam_Files/"+chipqc_df['ControlID']+'_final.bam'
chipqc_df['Peaks']="Results/peak_calling/"+chipqc_df['SampleID']+"_peaks.narrowPeak"
chipqc_df['PeakCaller']="narrow"
chipqc_df['Tissue']="NA"
chipqc_df['Condition']=chipqc_df["Factor"]
chipqc_df=chipqc_df[["SampleID","Factor","Replicate","bamReads","ControlID","bamControl","Peaks","PeakCaller","Tissue","Condition"]]
os.system("mkdir -p Results/ChIPQC")
chipqc_df.to_csv("Results/ChIPQC/metadata_for_ChIPQC.csv",index=None)



#https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html
os.system(f'Rscript ../../master_codes/chipseq/chipqc.R --metadata Results/ChIPQC/metadata_for_ChIPQC.csv')