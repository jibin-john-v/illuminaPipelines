#References: https://bioinfo-dirty-jobs.github.io/rana2//lectures/02.linux/
            #https://www.nature.com/articles/nprot.2016.095
            #https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/
            #https://rpubs.com/selveyad/
            #https://www.waccbip.org/research/scientific-technology-systems/high-performance-computing


os.system("export PATH=/Users/JJOHN41/software/samtools/bin:$PATH")

import sys,os,glob,subprocess
from pathlib import Path
import pandas as pd
import concurrent.futures


sys.path.insert(0, '../../illumina_mastercode/') 
from fastq_qc.fastq_qc import falco_run,fastp_run
from utility.utility import fastq_files_location_df,fastq_files_location_df
from alignment.alignment import hisat2aligner_run

from postalignment.postalignment  import bamindex_run
from rnaseq_specific.qunatification import quant_stringtie_run,merge_stringtie_run,gffcompare_run,quant_stringtie_run2




genome_index="/Users/JJOHN41/Desktop/projects/test_data/rnaseq/chrX_data/indexes/chrX_tran"
threads=10
binSize=20
fastq_path = "/Users/JJOHN41/Desktop/projects/test_data/rnaseq/chrX_data/samples/"
gtf_file="/Users/JJOHN41/Desktop/projects/test_data/rnaseq/chrX_data/genes/chrX.gtf"
paired_end, single_end = fastq_files_location_df(fastq_path)


fastq_df=paired_end.copy()


#run_fastqc_parallel("falco_run",single_end, n_parallel=5,n_threads=3,out_name="output",qc_cutoff=20)
for index in range(fastq_df.shape[0]):
    falco_run(fastq_df.loc[index,"SampleName"],[fastq_df.loc[index,"R1"],fastq_df.loc[index,"R2"]],"falco")


#run_fastqc_parallel("fastp_run",single_end, n_parallel=5,n_threads=3,out_name="output",qc_cutoff=20)
for index in range(fastq_df.shape[0]):
    fastp_run(fastq_df.loc[index,"SampleName"],[fastq_df.loc[index,"R1"],fastq_df.loc[index,"R2"]],5,20)


#Alignment using  hisat2
for index in range(single_end.shape[0]):
    s_name=single_end.loc[index,"SampleName"]
    hisat2aligner_run(fastq_df.loc[index,"SampleName"],
                        [fastq_df.loc[index,"R1"],fastq_df.loc[index,"R2"]],threads,genome_index,"FR")

#Assemble and quantify expressed genes and transcripts 
for index in range(fastq_df.shape[0]):
    s_name=fastq_df.loc[index,"SampleName"]
    quant_stringtie_run(s_name,gtf_file,threads,f"Results/bamfiles/{s_name}.bam")

##Merge gtf files
merge_list=pd.DataFrame(glob.glob("Results/quantification/*gtf"))
merge_list.to_csv("merge_list.txt",index=None)
merge_stringtie_run(gtf_file,threads,"merge_list.txt")

##Compare reference and merged gtf files
gffcompare_run(gtf_file,"Results/quantification/stringtie_merged.gtf")


##Run the stringtie using stringtie_merged.gtf file and bam file
for index in range(fastq_df.shape[0]):
    s_name=fastq_df.loc[index,"SampleName"]
    quant_stringtie_run2(s_name,"Results/quantification/stringtie_merged.gtf",threads,f"Results/bamfiles/{s_name}.bam")

