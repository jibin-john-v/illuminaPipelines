

#Reference: STAR: https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I ;
#                 https://www.reneshbedre.com/blog/star-aligner-twopass-mode.html

#Good tutorials ; https://nbisweden.github.io/workshop-RNAseq/2011/lab_smallrna.html 
#                 https://faqs.lexogen.com/faq/how-can-i-analyze-my-small-rna-sequencing-data

import os

#STAR Alignment
def staraligner_run(sample_name, reads, threads, genomedir):
    bam_dir="Results/bamfiles/"
    os.system(f"mkdir -p {bam_dir}")
    os.system(f"mkdir -p Results/log_files/")
    if len(reads) > 1:
        os.system(f'''STAR --genomeDir {genomedir} \
                   --readFilesIn {reads[0]} {reads[1]} \
                   --readFilesCommand zcat \
                   --outSAMunmapped Within \
                   --outFileNamePrefix {bam_dir}/{sample_name}.\
                   --runThreadN {threads} > Results/log_files/alignment_star_{sample_name}.log 2>&1''')
    else:
        os.system(f'''STAR --genomeDir {genomedir} \
                   --readFilesIn {reads[0]} \
                   --readFilesCommand zcat \
                   --outSAMunmapped Within \
                   --outFileNamePrefix {bam_dir}/{sample_name}.\
                   --runThreadN {threads} > Results/log_files/alignment_star_{sample_name}.log 2>&1''')        


##Alignment using BWA
def bwaaligner_run(sample_name, reads, threads, fasta):
    bam_dir="Results/bamfiles/"
    os.system(f"mkdir -p Results/log_files/")
    os.system(f"mkdir -p {bam_dir}")
    ID = f"HTYYCBBXX.1.{sample_name}"
    SM = sample_name
    LB = "Lib-1"
    PU = "HTYYCBBXX.1"
    PL = "ILLUMINA"
    if len(reads) > 1:
        os.system(f"""bwa mem \
                -R "@RG\\tID:{ID}\\tSM:{SM}\\tLB:{LB}\\tPU:{PU}\\tPL:{PL}" \
                -M \
                -t {threads} {fasta} {reads[0]} {reads[1]} >Results/log_files/alignment_bwamem_{sample_name}.log \
                | samtools sort -@ {threads} \
                -o Results/bamfiles/{sample_name}.bam """)
    else:
        os.system(f"""bwa mem \
                -R "@RG\\tID:{ID}\\tSM:{SM}\\tLB:{LB}\\tPU:{PU}\\tPL:{PL}" \
                -M \
                -t {threads} {fasta} {reads[0]} >Results/log_files/alignment_bwamem_{sample_name}.log  \
                | samtools sort -@ {threads} \
                -o Results/bamfiles/{sample_name}.bam """)


##Alignment using BWA-mem2 ; https://github.com/bwa-mem2/bwa-mem2
def bwa_mem2ligner_run(sample_name, reads, threads, fasta):
    bam_dir="Results/bamfiles/"
    os.system(f"mkdir -p {bam_dir}")
    os.system(f"mkdir -p Results/log_files/")
    ID = f"HTYYCBBXX.1.{sample_name}"
    SM = sample_name
    LB = "Lib-1"
    PU = "HTYYCBBXX.1"
    PL = "ILLUMINA"
    if len(reads) > 1:
        os.system(f"""bwa-mem2 mem \
                -R "@RG\\tID:{ID}\\tSM:{SM}\\tLB:{LB}\\tPU:{PU}\\tPL:{PL}" \
                -M \
                -t {threads} {fasta} {reads[0]} {reads[1]} 2>Results/log_files/alignment_bwamem2_{sample_name}.log \
                | samtools sort -@ {threads} \
                -o Results/bamfiles/{sample_name}.bam """)
    else:
        os.system(f"""bwa-mem2 mem \
                -R "@RG\\tID:{ID}\\tSM:{SM}\\tLB:{LB}\\tPU:{PU}\\tPL:{PL}" \
                -M \
                -t {threads} {fasta} {reads[0]} 2> Results/log_files/alignment_bwamem2_{sample_name}.log \
                | samtools sort -@ {threads} \
                -o Results/bamfiles/{sample_name}.bam """)



##Alignment using bowtie1
def bowtie1_srna_aligner_run(sample_name,reads,index,threads):
    os.system(' mkdir -p Results/Bam_Files')
    os.system(f"mkdir -p Results/log_files/")
    if len(reads) > 1:
        os.system(f'''bowtie \
                        -q \
                        -v 1 \
                        -k 10 \
                        -S \
                        -t \
                        -p {threads} \
                        --seedlen 16 \
                        -x {index} \
                        -1 {reads[0]} -2 {reads[1]} 2> Results/Bam_Files/bowtie1_alignment_{sample_name}.log  \
                        |samtools sort \
                        --output-fmt BAM \
                        -o Results/Bam_Files/{sample_name}.bam ''')
    else:
        os.system(f'''bowtie \
                -q \
                -v 1 \
                -k 10 \
                -S \
                -t \
                -p {threads} \
                --seedlen 16 \
                -x {index} {reads[0]} 2> Results/Bam_Files/bowtie1_alignment_{sample_name}.log \
                |samtools sort \
                --output-fmt BAM  \
                -o Results/Bam_Files/{sample_name}.bam ''')

##Alignment using hisat2 ; https://daehwankimlab.github.io/hisat2/ 
def hisat2aligner_run(sample_name, reads, threads, fasta,strand="FR"):
    bam_dir="Results/bamfiles/"
    os.system(f"mkdir -p {bam_dir}")
    os.system(f"mkdir -p Results/log_files/")
    ID = f"HTYYCBBXX.1.{sample_name}"
    SM = sample_name
    LB = "Lib-1"
    PU = "HTYYCBBXX.1"
    PL = "ILLUMINA"
    if len(reads) > 1:
        os.system(f"""hisat2 \
                        -p {threads} \
                        --rg-id={ID} --rg SM:{SM} --rg LB:{LB} --rg PL:{PL} --rg PU:{PU} \
                        -x {fasta} \
                        --dta \
                        --rna-strandness {strand} \
                        -1 {reads[0]} \
                        -2 {reads[1]} \
                         2> Results/log_files/alignment_hisat2_{sample_name}.log \
                        | samtools sort -O BAM \
                        | tee Results/bamfiles/{sample_name}.bam \
                        | samtools index - Results/bamfiles/{sample_name}.bai """)
    else:
        os.system(f"""hisat2 \
                        -p {threads} \
                        --rg-id={ID} --rg SM:{SM} --rg LB:{LB} --rg PL:{PL} --rg PU:{PU} \
                        -x {fasta} \
                        --dta \
                        --rna-strandness {strand} \
                        -1 {reads[0]} \
                         2> Results/log_files/alignment_hisat2_{sample_name}.log \
                        | samtools sort -O BAM \
                        | tee Results/bamfiles/{sample_name}.bam \
                        | samtools index - Results/bamfiles/{sample_name}.bai""")






##Alignment using bowtie2
def bowtie2_aligner_run(sample_name,reads,index,threads):
    os.system(' mkdir -p Results/Bam_Files')
    ID = f"HTYYCBBXX.1.{sample_name}"
    SM = sample_name
    LB = "Lib-1"
    PU = "HTYYCBBXX.1"
    PL = "ILLUMINA"
    if len(reads) > 1:
        os.system(f'''bowtie \
                        -q \
                        --local \
                       --rg-id={ID} --rg SM:{sample_name} --rg LB:{LB} --rg PL:{PL} --rg PU:{PU} \
                        -p {threads} \
                        -x {index} \
                        -1 {reads[0]} -2 {reads[1]} 2> Results/log_files/bowtie2_alignment_{sample_name}.log  \
                        |samtools sort \
                        --output-fmt BAM \
                        -o Results/Bam_Files/{sample_name}.bam ''')
    else:
        os.system(f'''bowtie2 \
                -q \
                --local \
                --rg-id={ID} --rg SM:{sample_name} --rg LB:{LB} --rg PL:{PL} --rg PU:{PU} \
                --threads {threads} \
                -x {index}  \
                -U {reads[0]} 2> Results/log_files/bowtie2_alignment_{sample_name}.log \
                |samtools sort \
                --output-fmt BAM  \
                -o Results/Bam_Files/{sample_name}.bam ''')



##Alignment using mirdeep2 ; https://daehwankimlab.github.io/hisat2/

def mirdeep2_aligner_run(sample_name, reads, bowtie_index, minimum_readlength=17):
    bam_dir = "Results/bamfiles/"
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs("Results/log_files/", exist_ok=True)
    global compressed
    
    def uncompress_if_gzipped(file_path):
        def is_gzipped(file_path):
            with open(file_path, 'rb') as f:
                return f.read(2) == b'\x1f\x8b'
        global compressed 
        if is_gzipped(file_path):
            output_file_path = file_path.rstrip('.gz')
            os.system(f"gunzip {file_path}")
            compressed = "gzipped"
            return output_file_path
        else:
            compressed = "notgzipped"
            return file_path
    reads = uncompress_if_gzipped(reads)
    
    os.system(f'''/mnt/disks/sdc/Software/mirdeep2/bin/mapper.pl \
                  {reads} \
                  -e \
                  -h \
                  -i \
                  -j \
                  -l {minimum_readlength} \
                  -m \
                  -p {bowtie_index} \
                  -s {sample_name}_reads_collapsed.fa \
                  -t {sample_name}_reads_vs_refdb.arf \
                  -v -o 4 2>Results/log_files/alignment_mirdeep2_{sample_name}.log''')
    if compressed == "gzipped":
        os.system(f"gzip {reads}")



##Align using bismark using bowtie2
def bismark_bowtie2aligner(reads, bismark_index):
    bam_dir = "Results/bamfiles/"
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs("Results/log_files/", exist_ok=True)
    os.system(f'''bismark \
                        --multicore 4 \
                        --bowtie2 \
                        {bismark_index} \
                        -1 {reads[0]} \
                        -2 {reads[1]} \
                        --bam ''')
    #for single end
    #os.system(f'''bismark --genome {bismark_index} {reads[0]}''')


#sortmerna ; to remove rRNA sequence