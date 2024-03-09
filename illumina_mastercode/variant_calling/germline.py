##---------------------------------------Variant calling-----------------------------------------------------------------------"

                            #https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller ; 
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360057439091-StrandBiasBySample
                            #When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360035890551-Allele-specific-annotation-and-filtering-of-germline-short-variants
mkdir -p Gvcf_Files
gatk --java-options "-Djava.io.tmpdir=`pwd` -Xms20G -Xmx20G -XX:ParallelGCThreads=2" HaplotypeCaller \
  -R ${Fasta} \
  -I BamFiles/${sample_name}_recal.bam \
  -O Gvcf_Files/${sample_name}.g.vcf.gz \
  -ERC GVCF \
  -A StrandBiasBySample \
  --native-pair-hmm-threads ${threads} \
  -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation


for j in {1..22} X Y M; do ; \

gatk \
    --java-options "-Djava.io.tmpdir=`pwd` -Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenomicsDBImport \
    --genomicsdb-workspace-path genomicsDBI_Database/chr${j}_gdb \
    -R ${Fasta} \
    -V NA12878.g.vcf.gz \
    -V NA12891.g.vcf.gz \
    -V NA12892.g.vcf.gz \
    --tmp-dir `pwd` \
    --max-num-intervals-to-import-in-parallel 3 \
    --intervals chr${j} 
done 


for j in {1..22} X Y M; do ; \
cwd=`pwd`
gatk \
    --java-options "-Djava.io.tmpdir=`pwd` -Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    -R ${Fasta} \
    -V gendb:///${cwd}genomicsDBI_Database/chr${j}_gdb \
    -A StrandBiasBySample \
    -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
    -O Vcf_Files/{chr}.vcf
done



##deepvariant https://github.com/google/deepvariant

docker run \
  -v "YOUR_INPUT_DIR":"/input" \
  -v "YOUR_OUTPUT_DIR:/output" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \ **Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
  --ref=/input/YOUR_REF \
  --reads=/input/YOUR_BAM \
  --output_vcf=/output/YOUR_OUTPUT_VCF \
  --output_gvcf=/output/YOUR_OUTPUT_GVCF \
  --num_shards=$(nproc) \ **This will use all your cores to run make_examples. Feel free to change.**
  --logging_dir=/output/logs \ **Optional. This saves the log output for each stage separately.
  --haploid_contigs="chrX,chrY" \ **Optional. Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
  --par_regions_bed="/input/GRCh3X_par.bed" \ **Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
  --dry_run=false **Default is false. If set to true, commands will be printed out but not executed.


#https://github.com/google/deepvariant/blob/r1.6/docs/trio-merge-case-study.md
time sudo docker run \
  -v "${DIR}":"/data" \
  quay.io/mlin/glnexus:v1.2.7 \
  /usr/local/bin/glnexus_cli \
  --config DeepVariantWES \
  --bed "/data/${CAPTURE_BED}" \
  /data/HG004.g.vcf.gz /data/HG003.g.vcf.gz /data/HG002.g.vcf.gz \
  | sudo docker run -i google/deepvariant:${VERSION} bcftools view - \
  | sudo docker run -i google/deepvariant:${VERSION} bgzip -c \
  > ${DIR}/deepvariant.cohort.vcf.gz




#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/freebayes-dnaseq-workflow.html#gsc.tab=0
ls *-md-rg.bam > bam.fofn
freebayes-parallel \
   <(fasta_generate_regions.py ${ref}.fai 100000) 16 \
   --fasta-reference ${ref} \
   --bam-list bam.fofn  > output.vcf