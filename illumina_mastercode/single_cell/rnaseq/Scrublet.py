
import scrublet as scr
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

##Reference: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html
##           https://github.com/swolock/scrublet/tree/master

input_dir = '/Users/JJOHN41/Desktop/projects/test_data/single_cell_rnaseq/data/ctrl_raw_feature_bc_matrix_copy/'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))


scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)



b_df=pd.read_csv(f'{input_dir}barcodes.tsv',sep="\t",header=None)
data_values = b_df[0].tolist()

df=pd.DataFrame([data_values,doublet_scores,predicted_doublets]).T
df.set_index(0,inplace=True)
df.columns=["Doublet_score","Is_doublet"]
df.to_csv('/Users/JJOHN41/Desktop/projects/test_data/single_cell_rnaseq/data/ctrl_raw_feature_bc_matrix/scrublet_calls.tsv',sep="\t")