##Software needed mirdeep2_patch ; https://github.com/Drmirdeep/mirdeep2_patch
import gzip,os

#prepare the inputs form bowtie index preparation
os.system("perl -plane 's/\s+.+$//' /mnt/disks/sdc/GATK_Resource/Homo_sapiens_assembly38.fasta > /mnt/disks/sdc/GATK_Resource/Homo_sapiens_assembly38_modified.fasta") ##no whit space allowed
os.system("sed -i '/^>/!s/[RYKMSWBDHVN]/N/g'  /mnt/disks/sdc/GATK_Resource/Homo_sapiens_assembly38_modified.fasta")

##bowtie index 
bowtie1_index_run("/mnt/disks/sdc/GATK_Resource/Homo_sapiens_assembly38_modified.fasta",
                   "Homo_sapiens_assembly38_modified",
                   "/mnt/disks/sdc/GATK_Resource/",12)

##Extract mirna fasta file of specific species
def extract_mirna_fasta(inputfasta,species,output_fasta):
    os.system(f'perl /mnt/disks/sdc/Software/mirdeep2_patch/extract_miRNAs.pl \
                    {inputfasta} {species} > {output_fasta}')


##Extract hairpin mirna fasta file of specific species
extract_mirna_fasta("/mnt/disks/sdc/GATK_Resource/mirbase_hairpin.fa","hsa","/mnt/disks/sdc/GATK_Resource/hsa_mirbase_hairpin.fa")
extract_mirna_fasta("/mnt/disks/sdc/GATK_Resource/mirbase_mirbase.fa","hsa","/mnt/disks/sdc/GATK_Resource/hsa_mirbase.fa")
extract_mirna_fasta("/mnt/disks/sdc/GATK_Resource/mirbase_mirbase.fa","mmu,ptr","/mnt/disks/sdc/GATK_Resource/other_mirbase.fa")

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

##Qunatification and novel miRNA identification using mirdeep2 ; https://daehwankimlab.github.io/hisat2/
def mirdeep2_quantification(sample_name,bowtie_index,mature_mirna_fa,other_sp_mature_mirna_fa,hairpin_fa):
    os.system(f'/mnt/disks/sdc/Software/mirdeep2/bin/miRDeep2.pl \
            {sample_name}_reads_collapsed.fa  \
            {bowtie_index} \
            {sample_name}_reads_vs_refdb.arf \
            {mature_mirna_fa} \
            {other_sp_mature_mirna_fa} \
            {hairpin_fa} \
            -t hsa 2>Results/log_files/r{sample_name}_report.log')





mirdeep2_aligner_run("SRR17237669", 
                      "/mnt/disks/sdc/mirna_data/Analysis/Results/CleanedFastQ_Files/SRR17237669_trimmed.fq.gz", 
                      '/mnt/disks/sdc/GATK_Resource/Homo_sapiens_assembly38_modified',minimum_readlength=18)

mirdeep2_quantification("SRR17237669",'/mnt/disks/sdc/GATK_Resource/Homo_sapiens_assembly38_modified.fasta',
                           "/mnt/disks/sdc/GATK_Resource/hsa_mirbase.fa",
                           "/mnt/disks/sdc/GATK_Resource/other_mirbase.fa",
                           "/mnt/disks/sdc/GATK_Resource/hsa_mirbase_hairpin.fa")

