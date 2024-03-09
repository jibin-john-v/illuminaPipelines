import os
import subprocess

##Generate Methylation Profiles
def bismark_methylation(bismark_index,bam_file,read_type):
    bam_dir = "Results/bamfiles/"
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs("Results/log_files/", exist_ok=True)
    os.system(f' mkdir -p Results/Bam_Files')
    if read_type=="single":
        command=f'''bismark_methylation_extractor \
                            --cytosine_report \
                            --single-end \
                            --comprehensive \
                            --output_dir \
                            --bedGraph  \
                            --cytosine_report \
                            {bam_file} \
                            --genome_folder {bismark_index}''')
    elif read_type=="paired":
        command=f'''bismark_methylation_extractor \
                            --cytosine_report \
                            --paired-end \
                            --comprehensive \
                            --output_dir \
                            --bedGraph  \
                            --cytosine_report \
                            {bam_file} \
                            --genome_folder {bismark_index}''')
    else:
        print("Input options are not correct read_type should be either single or paired")
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command bismark_methylation executed successfully. for the sample {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Command bismark_methylation execution failed with error: {e} for the sample {bam_file}")

