##Import required packages
import os
import pandas as pd


##Filtering parameters
vaf_upper=0.60
vaf_lower= 0.40
minimum_depth=10

biotype_to_include=['protein_coding']
impact_to_include=['MODERATE','HIGH']

##Identify the number of threads to be used for the pipeline run
n_threads=4

##Identify input/output files and path to be used during the pipeline workflow
vcf_file="Tumour_32_a.mutect2.filtered.vcf.gz"
vcf_path="/Users/awenhasan/Desktop/Pipeline/Raw_Data/"
out_path="/Users/awenhasan/Desktop/Pipeline/Analysis/"


vep_path="/Users/awenhasan/miniconda3/envs/WES-pipeline/bin/vep"
genome_fasta="/Users/awenhasan/Desktop/WES-Code/Refrences/Homo_sapiens_assembly38.fasta"

###required input file from cosmic 
census_gene="/Users/awenhasan/Desktop/WES-Code/Refrences/Cosmic_CancerGeneCensus_v99_GRCh38.tsv"
Hallmarks_gene="/Users/awenhasan/Desktop/WES-Code/Refrences/Cosmic_CancerGeneCensusHallmarksOfCancer_v99_GRCh38.tsv"
civic_genes="/Users/awenhasan/Desktop/WES-Code/Refrences/civic_01-Feb-2024-GeneSummaries.tsv"


##------------------------------------create cosmic and civic gene dataframe---------------------------------------------
#census gene reading using pd.read_csv and selecting the required columns and removing duplicate rows
census_df=pd.read_csv(f'{census_gene}',sep="\t")
census_df=census_df[['GENE_SYMBOL','SOMATIC','GERMLINE','MOLECULAR_GENETICS','ROLE_IN_CANCER']].drop_duplicates()

#opening halmark data frame using panda to read, and removing duplicates from the data frame
h_df=pd.read_csv(Hallmarks_gene,sep="\t")
h_df=h_df[['GENE_SYMBOL']].drop_duplicates()
h_df["Cosmic_Hallmark"]="Yes"

#Merging two data frame from gene census and halmark based on gene symbol column
cosmic_genes_df=pd.merge(census_df,h_df,on="GENE_SYMBOL",how="outer")#the outer option enables inclusion of genes from both data set , "left"would take from census, "right" would take from halmark
#merging vep data frame and cosmic data frame based on column name "Nearest" from vep and "gene symbol" from Cosmic in accordance to vep datafarme

##civic
civic_df=pd.read_csv(civic_genes,sep="\t")
civic_df=civic_df[["name","gene_civic_url"]].drop_duplicates()
civic_df=civic_df.rename(columns={"name":"GENE_SYMBOL"}) #replacing name from civic "name" with "Gene_symbol"

##Merge civic and cosmic data frame 
civic_cosmic_genes_df=pd.merge(cosmic_genes_df,civic_df,on="GENE_SYMBOL",how="outer")

##############------------------------Sample wise analysis -------------------------------------------------------------

output_file_name=vcf_file.replace(".vcf.gz","").split("/")[-1]

##Split multi allelic variants to bialleic variants 
os.system(f'''/Users/awenhasan/Desktop/WES-Code/software/bcftools/bcftools norm \
        --atomize \
        --multiallelics -any \
        --rm-dup all \
        --threads 5 {vcf_path}{vcf_file} | bgzip -c > {out_path}{output_file_name}_normalise.vcf.gz''')

##annotate normalized vcf file
os.system(f"""
            docker run -t -i -v $HOME/vep_data:/data \
            -v /Users/awenhasan/:/Users/awenhasan/ ensemblorg/ensembl-vep vep \
            --input_file {out_path}{output_file_name}_normalise.vcf.gz \
            --output_file {out_path}{output_file_name}_normalise_vep.vcf \
            --vcf \
            --species homo_sapiens \
            --assembly GRCh38 \
            --force_overwrite \
            --stats_file {output_file_name}.html \
            --warning_file {output_file_name}.error \
            --fork {n_threads} \
            --cache \
            --dir_cache $HOME/.vep/ \
            --fasta {genome_fasta} \
            --sift b --polyphen b \
            --nearest symbol \
            --overlaps \
            --gene_phenotype \
            --total_length \
            --numbers \
            --hgvs \
            --hgvsg \
            --protein \
            --symbol \
            --ccds \
            --uniprot \
            --tsl \
            --appris \
            --canonical \
            --mane \
            --mane_select \
            --biotype \
            --domains \
            --check_existing \
            --check_svs \
            --clin_sig_allele 1 \
            --no_check_alleles \
            --af \
            --max_af \
            --af_1kg \
            --af_gnomade \
            --af_gnomadg \
            --pubmed \
            --var_synonyms \
            --flag_pick_allele_gene \
            --failed 1 \
            --verbose \
            --buffer_size 500 \
            --custom file=/Users/awenhasan/Desktop/WES-Code/Refrences/01-Feb-2024-civic_accepted_hg38.vcf.gz,short_name=CIVIC,format=vcf,type=exact,coords=0,fields=GN%VT  \
            --custom file=/Users/awenhasan/Desktop/WES-Code/Refrences/Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz,short_name=COSMIC,format=vcf,type=exact,coords=0,fields=LEGACY_ID%SAMPLE_COUNT%TIER
              """)


#convert annotated vep vcf to tsv format
os.system(f'''vep-annotation-reporter {out_path}{output_file_name}_normalise_vep.vcf Allele Consequence IMPACT SYMBOL Gene Feature_type Feature BIOTYPE EXON INTRON HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation DISTANCE STRAND FLAGS PICK SYMBOL_SOURCE HGNC_ID CANONICAL MANE_SELECT MANE_PLUS_CLINICAL TSL APPRIS CCDS ENSP SWISSPROT TREMBL UNIPARC UNIPROT_ISOFORM SOURCE GENE_PHENO NEAREST SIFT PolyPhen DOMAINS HGVS_OFFSET HGVSg AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF gnomADe_AF gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF gnomADe_NFE_AF gnomADe_OTH_AF gnomADe_SAS_AF gnomADg_AF gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF gnomADg_OTH_AF gnomADg_SAS_AF MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC PHENO PUBMED SV OverlapBP OverlapPC VAR_SYNONYMS COSMIC COSMIC_LEGACY_ID COSMIC_SAMPLE_COUNT COSMIC_TIER \
          -o {out_path}{output_file_name}_normalise_vep.tab''')

####-----------------Read vcf file and select the required columns and create VAF,Depth,Ref_Depth and Alt_Depth----------------
os.system(f'''grep -v "##" {out_path}{output_file_name}_normalise_vep.vcf > {out_path}{output_file_name}_normalise_vep.vcf.tsv ''')
vcf_df=pd.read_csv(f'{out_path}{output_file_name}_normalise_vep.vcf.tsv',sep="\t")
vcf_df2=vcf_df[['#CHROM', 'POS','REF', 'ALT',"FILTER","FORMAT",vcf_df.columns[-1]]]
vcf_df2['VAF']=vcf_df2[vcf_df2.columns[-1]].str.split(":",expand=True)[2]
vcf_df2['Depth']=vcf_df[vcf_df.columns[-1]].str.split(":",expand=True)[3]
vcf_df2[["Ref_Depth","Alt_Depth"]]=vcf_df[vcf_df.columns[-1]].str.split(":",expand=True)[1].str.split(",",expand=True)
vcf_df2=vcf_df2.rename(columns={"#CHROM":"CHROM"})


##Read vep fle
#reading the test.tsv using .read_csv from pandas and saving variable into vep_df
vep_df=pd.read_csv(f'''{out_path}{output_file_name}_normalise_vep.tab''',sep="\t")

##merge vep annotated dataframe with vcf columns
vep_df=pd.merge(vep_df,vcf_df2, on=['CHROM', 'POS', 'REF', 'ALT'],how="outer")


####### merge vep annotated data frame with cosmic civic dataframe 
vep_cosmic_df=pd.merge(vep_df,civic_cosmic_genes_df,left_on='NEAREST',right_on='GENE_SYMBOL',how="left")
vep_cosmic_df.to_csv(f'{out_path}{output_file_name}_normalise_vep_cosmicgenes.csv',index=None)


##Filtered variants by VAF and depth
vep_cosmic_vaf_filt_df=vep_cosmic_df[ (vep_cosmic_df['VAF']<=vaf_upper) & 
                                         (vep_cosmic_df['VAF']>=vaf_lower) & 
                                         (vep_cosmic_df['Depth']>=minimum_depth) ]
vep_cosmic_vaf_filt_df.to_csv(f'{out_path}{output_file_name}_normalise_vep_cosmicgenes_vafDP_Filtered.csv',index=None)



###Filtered variants by variant impact/functional class
vep_cosmic_vaf_function_filt_df=vep_cosmic_vaf_filt_df[ (vep_cosmic_vaf_filt_df['BIOTYPE'].isin(biotype_to_include) ) &
                        (vep_cosmic_vaf_filt_df['IMPACT'].isin(impact_to_include) )]

vep_cosmic_vaf_function_filt_df.to_csv(f'{out_path}{output_file_name}_normalise_vep_cosmicgenes_vafDP_Function_Filtered.csv',index=None)

