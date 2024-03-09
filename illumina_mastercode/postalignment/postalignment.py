import os

def sam_tobam_run(sample_name, samfile, threads, outdir):
    os.system(f"mkdir -p {outdir}")
    os.system(f"samtools view -@ {threads} -o {outdir}/{sample_name}.bam {samfile}")
    os.system(f'samtools index {sample_name}.bam {sample_name}.bai')

def sort_bam_run(sample_name, bamfile, threads, outdir):
    os.system(f"mkdir -p {outdir}")
    os.system(f"samtools sort -@ {threads} -o {outdir}/{sample_name}_sorted.bam {bamfile}")
    os.system(f'samtools index {sample_name}_sorted.bam {sample_name}_sorted.bai')

def bamindex_run(bamfile):
    bam_index=bamfile.replace(".bam",".bai")
    os.system(f'samtools index {bamfile} {bam_index}')

def qualimap_rna_run(sample_name:str,bam_file,gtf_file):
    qualimap_dir=f"Results/QC/post_alignment/qualimap"
    qualimap_dir = f"{qualimap_dir}/{sample_name}"
    os.system(f"mkdir -p {qualimap_dir} ")
    os.system(f'''qualimap rnaseq \
                        -bam  {bam_file} \
                        -gtf {gtf_file} \
                        -outfile {sample_name} \
                        -outdir {qualimap_dir} \
                        --java-mem-size=8G > Results/log_files/post_alignment-qc_qualimap_{sample_name}.log 2>&1 ''')


def qualimap_genome_run(sample_name:str,bam_file,bedfile,threads):
    qualimap_dir=f"Results/QC/post_alignment/qualimap"
    qualimap_dir = f"{qualimap_dir}/{sample_name}"
    os.system(f"mkdir -p {qualimap_dir} ")
    os.system(f'''qualimap bamqc \
                        -bam BamFiles/{bam_file} \
                        --java-mem-size=8G
                        -nt ${threads} \
                        -outdir {qualimap_dir} \
                        -outfile ${sample_name} \
                        --feature-file {bedfile}''')


##samtools based flagstat & idxstats
def idx_flagstat_run(sample_name:str,bam_file,threads):
    flagidx_dir=f"Results/QC/post_alignment/flag_idxstat"
    flagidx_dir = f"{flagidx_dir}/{sample_name}"
    os.system(f"mkdir -p {flagidx_dir}")

    os.system(f'''samtools flagstat \
                        BamFiles/{bam_file} \
                        --threads {threads} \
                        --output-fmt tsv \
                        > {flagidx_dir}/{sample_name}_flagstat.tsv \
                         2> Results/log_files/post_alignment-qc_flagstat_{sample_name}.log''')
    
    os.system(f'''samtools idxstats \
                        BamFiles/{bam_file} \
                        --threads {threads} \
                        > {flagidx_dir}/{sample_name}_idxstats.tsv \
                        2> Results/log_files/post_alignment-qc_idxstats_{sample_name}.log''')



#https://rseqc.sourceforge.net/#
def rseqc_run(sample_name:str,bam_file,gtf_bed):
    rseqc_dir=f"Results/QC/post_alignment/rseqc"
    rseqc_dir = f"{rseqc_dir}/{sample_name}"
    
    os.system(f'''infer_experiment.py \
                  -r {gtf_bed} \
                  -i {bam_file} \
                  > {rseqc_dir}/{sample_name}_infer_experiment.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_infer_experiment_{sample_name}.log''')
    
    os.system(f'''geneBody_coverage.py \
                  -r {gtf_bed} \
                  -i {bam_file} \
                  -o {rseqc_dir}/{sample_name}_geneBody_coverage.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_geneBody_coverage_{sample_name}.log''')
    
    os.system(f'''junction_annotation.py \
                  -r {gtf_bed} \
                  -i {bam_file} \
                  -o {rseqc_dir}/{sample_name}_junction_annotation.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_junction_annotation_{sample_name}.log''')
    
    os.system(f'''read_distribution.py \
                  -r {gtf_bed} \
                  -i {bam_file} \
                  > {rseqc_dir}/{sample_name}_read_distribution.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_read_distribution_{sample_name}.log''')
    
    os.system(f'''read_duplication.py \
                  -i {bam_file} \
                  -o {rseqc_dir}/{sample_name}_read_duplication.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_read_duplication_{sample_name}.log''')
    
    os.system(f'''inner_distance.py \
                  -r {gtf_bed} \
                  -i {bam_file} \
                  -o {rseqc_dir}/{sample_name}_inner_distance.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_inner_distance_{sample_name}.log''')

    """os.system(f'''RPKM_saturation.py \
                  -r {gtf_bed} \
                  -i {bam_file} \
                  -d '1++,1--,2+-,2-+' \
                  -o {rseqc_dir}/{sample_name}_inner_distance.tsv 
                  2> Results/log_files/post_alignment-qc_rseqc_RPKM_saturation_{sample_name}.log''')"""




#dupRadar should introduce ;https://www.bioconductor.org/packages/release/bioc/vignettes/dupRadar/inst/doc/dupRadar.html


