

#Mutect2 ; https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2
##Good tutorials https://github.com/gatk-workflows/gatk4-jupyter-notebook-tutorials/tree/master/notebooks/Day3-Somatic
##https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/

#Tumor with matched normal
     gatk Mutect2 \
     -R reference.fa \
     -I tumor.bam \
     -I normal.bam \
     -normal normal_sample_name \
     --germline-resource af-only-gnomad.vcf.gz \
     --panel-of-normals pon.vcf.gz \
     -O somatic.vcf.gz



#


##CreateSomaticPanelOfNormals

#Step 1. Run Mutect2 in tumor-only mode for each normal sample.
 gatk Mutect2 \
   -R reference.fa \
   -I normal1.bam \
   -tumor normal1_sample_name \
   --germline-resource af-only-gnomad.vcf.gz \
   -O normal1_for_pon.vcf.gz


#Step 2. Create a file ending with .args or .list extension with the paths to the VCFs from step 1, one per line.
 gatk CreateSomaticPanelOfNormals \
   -vcfs normals_for_pon_vcf.args \
   -O pon.vcf.gz


 gatk CreateSomaticPanelOfNormals \
   -vcfs normal1_for_pon_vcf.gz \
   -vcfs normal2_for_pon_vcf.gz \
   -vcfs normal3_for_pon_vcf.gz \
   -O pon.vcf.gz