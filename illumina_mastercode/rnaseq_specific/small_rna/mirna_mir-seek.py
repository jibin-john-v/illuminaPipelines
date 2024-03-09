import os 


def mirna_mir_seek(fastq_folder,fast_suffix,output_folder):
    os.system(f'mkdir -p {output_folder}')
    os.system(f'''/mnt/disks/sdc/Software/mir-seek/mir-seek run \
                    --input {fastq_folder}/*{fast_suffix} \
                  --output {output_folder} \
                  --genome hg38 \
                  --mode local''')


mirna_mir_seek("/mnt/disks/sdc/mirna_data/Analysis/Results/CleanedFastQ_Files/",
                "_trimmed.fq.gz",
                "/mnt/disks/sdc/mirna_data/Analysis/Results/CleanedFastQ_Files/Result_mir_seek")