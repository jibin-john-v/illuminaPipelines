import os
import subprocess




########Mark duplicate unmapped,unequally mapped reads selection https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
def sambamba_unique_mdup_unmap_run(sample_name, bam_file, threads):
    os.system(f' mkdir -p Results/Bam_Files')
    os.system(f' mkdir -p Results/log_files')
    command =f'''sambamba view \
                -h \
                -t {threads} \
                -f bam \
                -F "[XS] == null and not unmapped and not duplicate" {bam_file} \
                > Results/Bam_Files/{sample_name}_final.bam 2> Results/log_files/post_alignment-qc_sambamba__dupuniquemultimap_{sample_name}.log'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command sambamba_unique_mdup_unmap_run executed successfully. for the sample {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Command execution failed with error: {e} for the sample {bam_file}")


##Main reference: https://hpc.nih.gov/training/gatk_tutorial/
#######Mark duplicate  https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
def gatk_markduplicate_run(sample_name:str,bam_file):
    os.system(f'''gatk \
                        --java-options "-Djava.io.tmpdir=`pwd` -Xms16G -Xmx16G" MarkDuplicates \
                        --INPUT BamFiles/{bam_file} \
                        --OUTPUT BamFiles/${sample_name}_Mduplicates.bam \
                        --METRICS_FILE BamFiles/{sample_name}_marked_dup_metrics.txt \
                        --CREATE_INDEX true --TMP_DIR `pwd`\
                        2> Results/log_files/post_alignment-qc_gatk_marked_duplicate_{sample_name}.log''')

##https://felixkrueger.github.io/Bismark/bismark/deduplication/
def deduplicate_bismark_run(sample_name:str,bam_file,read_type):
    os.system(f' mkdir -p Results/Bam_Files')
    os.system(f' mkdir -p Results/log_files')
    if read_type=="single":
          command=f'''deduplicate_bismark \
                            --single \
                            --outfile \
                            --output_dir \
                            --bam {bam_file} 2> Results/log_files/deduplicate_bismark_{sample_name}.log'''
    elif read_type=="paired":
        command =f'''deduplicate_bismark \
                            --paired \
                            --outfile \
                            --output_dir \
                            --bam {bam_file} 2> Results/log_files/deduplicate_bismark_{sample_name}.log'''
    else:
        print("Input options are not correct read_type should be either single or paired")
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command deduplicate_bismark_run executed successfully. for the sample {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Command deduplicate_bismark_run execution failed with error: {e} for the sample {bam_file}")


########-----------------------https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-------------------
### Base (Quality Score) Recalibration ; https://gatk.broadinstitute.org/hc/en-us/articles/360056969412-BaseRecalibrator" /https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator

def gatk_bqsr1_run(sample_name:str,bam_file,threads,dbsnp,GoldIndels,known_indels,fasta):
    os.system(f'mkdir -p Results/Recal_Files')
    os.system(f'''gatk \
                     --java-options "-Djava.io.tmpdir=`pwd` -Xms4G -Xmx4G -XX:ParallelGCThreads={threads}" \
                    BaseRecalibrator \
                    -I Results/BamFiles/{bam_file}\
                    -R {fasta} \
                    --known-sites {dbsnp} \
                    --known-sites {GoldIndels} \
                    --known-sites {known_indels} \
                    -O Results/Recal_Files/{sample_name}_recal_data.table \
                    2> Results/log_files/post_alignment-qc_gatk_bqsr_step1_{sample_name}.log''')


#ApplyBQSR ; https://gatk.broadinstitute.org/hc/en-us/articles/360056968652-ApplyBQSR"
def gatk_bqsr2_run(sample_name:str,bam_file,threads,fasta):
    os.system(f'mkdir -p Results/Recal_Files')
    os.system(f'''gatk \
                     --java-options "-Djava.io.tmpdir=`pwd` -Xms2G -Xmx2G -XX:ParallelGCThreads={threads}" \
                    ApplyBQSR \
                    -I Results/BamFiles/{bam_file}\
                    -R {fasta} \
                    --bqsr-recal-file Recal_Files/{sample_name}_recal_data.table \
                    -O Results/BamFiles/${sample_name}_recal.bam \
                    --create-output-bam-index true \
                    2> Results/log_files/post_alignment-qc_gatk_bqsr_step2_{sample_name}.log''')

#Create Post recal table #https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/"
def gatk_bqsr2_run(sample_name:str,threads,dbsnp,GoldIndels,known_indels,fasta):
    os.system(f'mkdir -p Results/Recal_Files')
    os.system(f'''gatk \
                     --java-options "-Djava.io.tmpdir=`pwd` -Xms2G -Xmx2G -XX:ParallelGCThreads={threads}" \
                    BaseRecalibrator \
                    -I Results/BamFiles/{sample_name}_recal.bam \
                    -R {fasta} \
                    --known-sites {dbsnp} \
                    --known-sites {GoldIndels} \
                    --known-sites {known_indels} \
                    --output  Results/Recal_Files/{sample_name}_post_recal_data.table
                    2> Results/log_files/post_alignment-qc_gatk_bqsr_step3_{sample_name}.log''')

#AnalyzeCovariates ; https://gatk.broadinstitute.org/hc/en-us/articles/360056967752-AnalyzeCovariates"
def gatk_bqsr3_run(sample_name:str,threads):
    os.system(f'mkdir -p Results/Recal_Files')
    os.system(f'''gatk \
                        --java-options "-Djava.io.tmpdir=`pwd` -Xms2G -Xmx2G -XX:ParallelGCThreads={threads}" \
                        AnalyzeCovariates \
                        -before Results/Recal_Files/{sample_name}_recal_data.table \
                        -after Results/Recal_Files/{sample_name}_post_recal_data.table \
                        -plots Results/Recal_Files/{sample_name}_AnalyzeCovariates.pdf \
                        2> Results/log_files/post_alignment-qc_gatk_bqsr_step4_{sample_name}.log ''')