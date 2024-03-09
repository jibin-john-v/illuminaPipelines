#Reference: STAR: https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I ;
#                 https://www.reneshbedre.com/blog/star-aligner-twopass-mode.html

#Good tutorials ; https://nbisweden.github.io/workshop-RNAseq/2011/lab_smallrna.html 
#                 https://faqs.lexogen.com/faq/how-can-i-analyze-my-small-rna-sequencing-data

import os
'''
genome_fasta="~/db/hg38/hg38.fa"
genomeDir="~/db/hg38/"
gtf_file="~/db/hg38/hg38.gtf"
read_length=100
sjdbOverhang=read_length-1
n_threads=5
'''

"""
#1) Indexing genome with annotations
os.system(f'''STAR --runMode genomeGenerate \
                  --genomeDir {genomeDir} \
                  --genomeFastaFiles {genome_fasta} \
                  --sjdbGTFfile {gtf_file} \
                  --runThreadN 30 \
                  --sjdbOverhang {sjdbOverhang}''')


#3) Indexing genome with annotations and SJ.out.tab files
os.system(f"cat *.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ_filtered.tab")
os.system(f'''STAR --runMode genomeGenerate \
                   --genomeDir ~/db/hg38/SJ_Index/ \
                   --genomeFastaFiles {genome_fasta} \
                   --sjdbGTFfile {gtf_file} \
                   --runThreadN 30 \
                   --sjdbOverhang 89 \
                   --sjdbFileChrStartEnd {SJ_filtered.tab} ''')
"""

##bowtie index 
def bowtie1_index_run(fasta,index_prefix,out_dir,threads):
    os.system(f'''mkdir -p {out_dir}''')
    os.system(f'''bowtie-build {fasta} \
                        {out_dir}/{index_prefix} \
                         --threads {threads}''')

##bowtie2 index 
def bowtie2_index_run(fasta,index_prefix,out_dir,threads):
    os.system(f'''mkdir -p {out_dir}''')
    os.system(f'''bowtie2-build {fasta} \
                        {out_dir}/{index_prefix} \
                         --threads {threads}''')

#hasat2 index
def hasat2_index_run2(gtf,fasta,out_dir):
    os.system(f" mkdir -p {out_dir}")
    gtf_out=gtf.split("/")[-1].replace(gtf,"")
    os.system(f"extract_splice_sites.py {gtf} > {out_dir}/{gtf_out}.ss")
    os.system(f"extract_exons.py {gtf} > {out_dir}/{gtf_out}.exon")
    os.system(f"hisat2-build \
                        --ss {out_dir}/{gtf_out}.ss \
                        --exon {out_dir}/{gtf_out}.exon \
                        {fasta}")