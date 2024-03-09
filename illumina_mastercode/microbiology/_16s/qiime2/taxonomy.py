import os

##Reference https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)

#Run taxonomic classification

def qiime_taxonomic_assignment(rep_sequences,reference,threads):
   '''Taxonomy assignment'''
   os.system(f'''qiime feature-classifier classify-sklearn \
                  --i-reads {rep_sequences} \
                  --i-classifier {reference} \
                  --p-n-jobs {threads} \
                  --o-classification taxonomy.qza''')


def qiime_taxa_export(taxa,rep_sequences):
   os.system(f'''qiime tools export \
                     --input-path {taxa} --output-path taxa''')
   os.system(f'''qiime metadata tabulate \
               --m-input-file {taxa} \
               --o-visualization taxonomy.qzv ''')
   os.system(f'''qiime feature-table tabulate-seqs \
                        --i-data {rep_sequences} \
                        --o-visualization representative_sequences.qzv''')



#classify-consensus-blast
#classify-consensus-vsearch

#