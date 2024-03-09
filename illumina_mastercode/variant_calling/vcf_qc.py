
###-------------------------------Variant Quality Score Recalibration (VQSR)-----------------------------------------
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-



##VariantRecalibrator 
    #https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator

gatk \
    --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 \
    -tranche 99.5 -tranche 99.0 -tranche 97.0 \
    -tranche 96.0 -tranche 95.0 -tranche 94.0 \
    -tranche 93.5 -tranche 93.0 -tranche 92.0 \
    -tranche 91.0 -tranche 90.0 \
    -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
    -V merged.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /fdb/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz  \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 /fdb/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 /fdb/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -an QD \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an FS -an SOR  \
    -mode SNP \
    -O merged_SNP1.recal \
    --tranches-file output_SNP1.tranches \
    --rscript-file output_SNP1.plots.R


gatk \
    --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 \
    -tranche 99.5 -tranche 99.0 -tranche 97.0 \
    -tranche 96.0 -tranche 95.0 -tranche 94.0 \
    -tranche 93.5 -tranche 93.0 -tranche 92.0 \
    -tranche 91.0 -tranche 90.0 \
    -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
    -V merged.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 /fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz \
    -an QD \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an FS \
    -an SOR \
    -an DP \
    -mode INDEL \
    -O merged_indel1.recal \
    --tranches-file output_indel1.tranches \
    --rscript-file output_indel1.plots.R


##Allele-specific version of the SNP filtering (beta)
    #https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR

gatk \
    --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V merged.vcf.gz \
    --recal-file merged_SNP1.recal \
    -mode SNP \
    --tranches-file output_SNP1.tranches \
    --truth-sensitivity-filter-level 99.9 \
    --create-output-variant-index true \
    -O SNP.recalibrated_99.9.vcf.gz


gatk \ 
    --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V SNP.recalibrated_99.9.vcf.gz \
    -mode INDEL \
    --recal-file merged_indel1.recal \
    --tranches-file output_indel1.tranches \
    --truth-sensitivity-filter-level 99.9 \
    --create-output-variant-index true \
    -O indel.SNP.recalibrated_99.9.vcf.gz

