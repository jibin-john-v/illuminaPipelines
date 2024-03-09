import os

##Alignment
        #Cell Ranger https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials
        #STARsolo https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
        #Alevin  https://salmon.readthedocs.io/en/latest/alevin.html 
        #UMI-tools  https://github.com/CGATOxford/UMI-tools
        #zUMIs  https://github.com/sdparekh/zUMIs
##Reff : https://bioinformatics-core-shared-training.github.io/SingleCell_RNASeq_Sept23/UnivCambridge_ScRnaSeqIntro_Base/Markdowns/03_CellRanger.html

##create index for cellranger
def cellranger_index(ref_fasta,gtf_file, threads,cellranger_index):
    os.makedirs("Results/log_files/", exist_ok=True)
    os.system(f'''cellranger mkref \
                --fasta={ref_fasta} \
                --genes={gtf_file} \
                --genome={cellranger_index} \
                --nthreads={threads}''')


##Align using bismark using bowtie2
def cellranger_aligner(fasq_folder, sample_name, fastq_sample_name,cellranger_index, threads):
    os.makedirs("Results/log_files/", exist_ok=True)
    os.system(f'''cellranger count \
                        --id={sample_name} \
                        --transcriptome={cellranger_index} \
                        --fastqs={fasq_folder} \
                        --sample={fastq_sample_name} \
                        --localcores={threads} \
                        --localmem=24''')