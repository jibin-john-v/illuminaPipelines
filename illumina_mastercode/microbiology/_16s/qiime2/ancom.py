import os

def qiime2_ancom_1(feature_table):
    os.system(f'''qiime composition add-pseudocount \
                    --i-table {feature_table} \
                    --p-pseudocount 1 \
                    --o-composition-table table_final_pseudocount.qza''')


def qiime2_ancom_1(metadata,):
    meta_df = pd.read_csv(f'{metadata}', sep="\t")
    groups = list(meta_df.columns[1:])
    
    for category in groups:
        os.system(f'''qqiime composition ancom \
              --i-table table_final_pseudocount.qza \
              --m-metadata-file {metadata}\
              --m-metadata-column {category} \
              --output-dir ancom_output''')