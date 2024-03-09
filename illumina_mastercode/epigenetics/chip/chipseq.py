# pylint: disable=missing-function-docstring
# Reference https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/schedule/links-to-lessons.md
#          #https://nbisweden.github.io/workshop-archive/workshop-ChIP-seq/2018-11-07/labs/lab-processing.html
           #https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/schedule/2-day.md

import os,subprocess
######## MACS2 based narrow peak calling https://github.com/macs3-project/MACS
def macs2_narrowpeak_call(sample_name, chip_bam, input_bam,effective_genome_size):
    os.system(f' mkdir -p Results/narrow_peak_calling')
    os.system(f' mkdir -p Results/log_files')
    command =f'''macs2 callpeak \
                        --treatment {chip_bam} \
                        --control {input_bam} \
                        --format BAM \
                        --gsize {effective_genome_size} \
                        --name {sample_name} \
                        --outdir Results/narrow_peak_calling 2>  Results/log_files/peak_calling_macs2_narrowpeak_{sample_name}.log'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command macs2_narrowpeak_call executed successfully. for the sample {sample_name}")
    except subprocess.CalledProcessError as e:
        print(f"Command macs2_narrowpeak_call execution failed with error: {e} for the sample {sample_name}")


######## MACS2 based broad peak calling https://github.com/macs3-project/MACS
def macs2_broadpeak_call(sample_name, chip_bam, input_bam,effective_genome_size):
    os.system(f' mkdir -p Results/broad_peak_calling')
    os.system(f' mkdir -p Results/log_files')
    command =f'''macs2 callpeak \
                        --treatment {chip_bam} \
                        --control {input_bam} \
                        --format BAM \
                        --gsize {effective_genome_size} \
                        --name {sample_name} \
                        --broad \
                        --outdir Results/broad_peak_calling 2>  Results/log_files/peak_calling_macs2_broadpeak{sample_name}.log'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command macs2_peak_call executed successfully. for the sample {sample_name}")
    except subprocess.CalledProcessError as e:
        print(f"Command macs2_broadpeak_call execution failed with error: {e} for the sample {sample_name}")

#remove blacklist regions
def remove_blacklist_region(blacklist_bed,peakfile,output_folder,out_name):
    os.system(f' mkdir -p Results/{output_folder}')
    os.system(f' mkdir -p Results/log_files')
    command =f'''bedtools intersect \
                            -v \
                            -a {peakfile} \
                            -b {blacklist_bed} \
                            > {output_folder}/{out_name} 2>  Results/log_files/remove_blacklist_region_{out_name}.log'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command  remove_blacklist_region executed successfully. for the input {peakfile}")
    except subprocess.CalledProcessError as e:
        print(f"Command remove_blacklist_region execution failed with error: {e} for the input {peakfile}")


#overlaping peak
def overlap_region(peak1, peak2,output_file,output_folder):
    os.system(f' mkdir -p Results/{output_folder}')
    os.system(f' mkdir -p Results/log_files')
    command =f'''bedtools intersect \
                    -wo -f 0.3 -r \
                    -a {peak1} \
                    -b {peak2} \
                    >  Results/{output_folder}/{output_file}'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command  remove_blacklist_regiion executed successfully. for the input {peak1} & {peak2}")
    except subprocess.CalledProcessError as e:
        print(f"Command remove_blacklist_regiion execution failed with error: {e} for the input {peak1} & {peak2}")

##bigWig files creation
def bigWig_filecreation(sample_name,bam_file,binSize,output_folder):
    os.system(f' mkdir -p Results/{output_folder}')
    os.system(f' mkdir -p Results/log_files')
    command =f'''bamCoverage \
            -b {bam_file} \
            -o Results/{output_folder}{sample_name}.bw \
            --binSize {binSize} 2>  Results/log_files/peak_calling_macs2_bamCoverage{sample_name}.log'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command  bigWig_filecreation executed successfully. for the input {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Command bigWig_filecreation execution failed with error: {e} for the input {bam_file}")

def bigWig_comparison_filecreation(chipsname,contsname,chip_bam,input_bam,binSize,output_folder):
    os.system(f' mkdir -p Results/{output_folder}')
    os.system(f' mkdir -p Results/log_files')
    command =f'''bamCompare \
            -b1 {chip_bam} \
            -b2 {input_bam} \
            -o Results/{output_folder}{chipsname}_{contsname}.bw \
            --binSize {binSize} 2>  Results/log_files/peak_calling_macs2_bamCompare_{chipsname}_{contsname}.log'''
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command  bigWig_comparison_filecreation executed successfully. for the input {chipsname}_{contsname}")
    except subprocess.CalledProcessError as e:
        print(f"Command bigWig_comparison_filecreation execution failed with error: {e} for the input {chipsname}_{contsname}")


#Create the matrix
'''computeMatrix reference-point --referencePoint center \
-b 4000 -a 4000 \
-R ~/chipseq_workshop/results/macs2/wt_peaks_final.bed \
-S visualization/bigWig/wt_sample1_chip.bw visualization/bigWig/wt_sample2_chip.bw \
--skipZeros \
-o ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-p 6'''


#Drawing the profile plot
'''
plotProfile -m ~/chipseq_workshop/results/visualization/wt_matrix.gz \
-out ~/chipseq_workshop/results/visualization/figures/plot1_wt_replicates.png \
--regionsLabel "" \
--perGroup \
--colors red blue \
--samplesLabel "WT_replicate1" "WT_replicate2" \
--refPointLabel "PRDM16 binding sites"
'''



##visualization   ;https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html
                #https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html






#Handling replicates using the Irreproducibility Discovery Rate (IDR) framework ;https://github.com/karmel/homer-idr


##Other peak callers

#HOMER: a suite of tools for peak calling and motif discovery. http://homer.ucsd.edu/homer/index.html
#SPP: an R package, that is implemented in the ENCODE processing pipeline. Best for narrow peak calling. Uses a sliding window to calculate scores based on fragment counts from up- and downstream flanking windows.
        #https://hbctraining.github.io/Intro-to-ChIPseq/lessons/peak_calling_spp.html
        #https://github.com/hms-dbmi/spp

#epic2: Ideal for broad peak calling (a re-implementation of an older tool called SICER)
        #https://github.com/biocore-ntnu/epic2  

#haystack bio: Epigenetic Variability and Motif Analysis Pipeline
        #https://github.com/pinellolab/haystack_bio
