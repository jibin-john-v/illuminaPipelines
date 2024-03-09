
import os
import pandas as pd 

#df=pd.read_csv("/Users/JJOHN41/Desktop/projects/microbiome/data/manifest_file.csv")
#df.to_csv("manifest_file.tsv",index=None,sep="\t")
##Import reference ; https://docs.qiime2.org/2023.9/tutorials/importing/

def qiime2_fastq_import(manifest_file,readtype):
    m_paired_columns=['sample-id','forward-absolute-filepath','reverse-absolute-filepath']
    m_singleend_columns=['sample-id','absolute-filepath']
    manifiest_df=pd.read_csv(manifest_file,sep="\t")
    
    if readtype=="paired":
        if all(item in list(manifiest_df.columns) for item in m_paired_columns):
            os.system(f'''
                    qiime tools import \
                    --type 'SampleData[PairedEndSequencesWithQuality]' \
                    --input-path {manifest_file} \
                    --output-path fastq_imported.qza \
                    --input-format PairedEndFastqManifestPhred33V2''')
        else:
            return 'The manifest file columns names are wrong'
    
    elif readtype=="singleend":
        if all(item in list(manifiest_df.columns) for item in m_singleend_columns):
            os.system(f'''
                        qiime tools import \
                        --type 'SampleData[SequencesWithQuality]' \
                        --input-path {manifest_file} \
                        --output-path fastq_imported.qza \
                        --input-format SingleEndFastqManifestPhred33V2''')
        else:
            return 'The manifest file columns names are wrong'
        



#qiime2_fastq_import("manifest_file.tsv","paired")
#qiime2_fastq_summarize("fastq_imported.qza")
