import os
import pandas as pd


#https://docs.qiime2.org/2023.9/tutorials/moving-pictures/
##Reference based tree using SEPP 

def qiime2_taxa_barplot(feature_table,taxonomy,metadata):
   '''taxonomic composition barplots'''
   os.system(f'''qiime taxa barplot \
                    --i-table {feature_table} \
                    --i-taxonomy {taxonomy} \
                    --m-metadata-file {metadata} \
                    --o-visualization taxa-bar-plots-1.qzv''')


##de novo tree creation
#Generate a tree for phylogenetic diversity analyses
def qiime2_mafft_fasttree(repseqs):
    os.system(f'''qiime phylogeny align-to-tree-mafft-fasttree \
                  --i-sequences {repseqs} \
                  --o-alignment aligned-rep-seqs.qza \
                  --o-masked-alignment masked-aligned-rep-seqs.qza \
                  --o-tree unrooted-tree.qza \
                  --o-rooted-tree rooted-tree.qza''')

#Alpha and beta diversity analysis¶
def qiime2_alpha_beta(rootedtree,tableqza,sdepth,metadata):
    os.system(f'''qiime diversity core-metrics-phylogenetic \
                --i-phylogeny {rootedtree} \
                --i-table {tableqza} \
                --p-sampling-depth {sdepth} \
                --m-metadata-file {metadata} \
                --output-dir core-metrics-results''')


##Alpha rarefaction plotting
def qiime2_alpha_rarefaction(feature_table,rooted_tree,metadata,maxdepth):
   os.system(f'''qiime diversity alpha-rarefaction \
        --i-table {feature_table} \
        --i-phylogeny {rooted_tree} \
        --p-max-depth {maxdepth} \
        --m-metadata-file {metadata} \
        --o-visualization alpha-rarefaction.qzv''')
   
   os.system(f'''qiime diversity alpha-rarefaction \
      --i-table {feature_table} \
      --p-max-depth {maxdepth} \
      --p-steps 20 \
      --i-phylogeny {rooted_tree} \
      --o-visualization rarefaction_curves_eachsample.qzv''')


#stacked barchart of taxa relative abundances ; https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)

def qiime2_barchart(tableqza, metadata, taxa):
    os.system(f'''qiime taxa barplot \
                 --i-table {tableqza} \
                 --i-taxonomy {taxa} \
                 --m-metadata-file {metadata} \
                 --o-visualization taxa_barplot.qzv''')
    meta_df = pd.read_csv(f'{metadata}', sep="\t")
    groups = list(meta_df.columns[1:])
    
    for category in groups:
        os.system(f'''qiime feature-table group \
                       --i-table {tableqza} \
                       --p-axis sample \
                       --p-mode sum \
                       --m-metadata-file {metadata} \
                       --m-metadata-column {category} \
                       --o-grouped-table taxa_{category}.qza''')
        
        os.system(f'''qiime taxa barplot \
                       --i-table taxa_{category}.qza \
                       --i-taxonomy {taxa} \
                       --o-visualization taxa_barplot_{category}.qzv''')



##Alpha group significance
def qiime2_alpha_comparison(alpha_path,metadata):
    alpha_diversities = ['evenness_vector.qza','faith_pd_vector.qza','observed_features_vector.qza','shannon_vector.qza']
    
    for alpha_diversity in alpha_diversities:
        out_name=alpha_diversity.replace(".qza","_alpha_diversity.qzv")
        os.system(f'''qiime diversity alpha-group-significance \
                        --i-alpha-diversity {alpha_path}{alpha_diversity} \
                        --m-metadata-file {metadata} \
                        --o-visualization {out_name}''')
    

def qiime2_beta_comparison(beta_path,metadata,feature_table,taxonomy):
    beta_diversities = ['bray_curtis_distance_matrix.qza','unweighted_unifrac_distance_matrix.qza',
                        'jaccard_distance_matrix.qza','weighted_unifrac_distance_matrix.qza']
    
    for beta_diversity in beta_diversities:
        out_name=beta_diversity.replace(".qza","_beta_diversity_umap.qza")
        os.system(f'''qiime diversity umap \
                        --i-distance-matrix {beta_path}{beta_diversity} \
                        --o-umap {out_name}''')
    
    os.system(f'''qiime metadata tabulate \
              --m-input-file {metadata} unweighted_unifrac_distance_matrix_beta_diversity_umap.qza {beta_path}/faith_pd_vector.qza {beta_path}/evenness_vector.qza {beta_path}/shannon_vector.qza \
              --o-visualization expanded-metadata-summ.qzv''')
    
    os.system(f'''qiime taxa barplot \
            --i-table {feature_table} \
            --i-taxonomy {taxonomy} \
            --m-metadata-file {metadata} unweighted_unifrac_distance_matrix_beta_diversity_umap.qza \
            {beta_path}/faith_pd_vector.qza {beta_path}/evenness_vector.qza {beta_path}/shannon_vector.qza \
            --o-visualization taxa-bar-plots-2.qzv''')







'''
#microbial composition of the samples in the context of the sample metadata.
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv


##sample composition in the context of categorical metadata using PERMANOVA 
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column body-site \
  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv \
  --p-pairwise


##explore principal coordinates (PCoA) plots in the context of sample metadata
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv


##Taxonomic analysis
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza




##Differential abundance testing with ANCOM
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where "[body-site]='gut'" \
  --o-filtered-table gut-table.qza

qiime composition add-pseudocount \
  --i-table gut-table.qza \
  --o-composition-table comp-gut-table.qza

qiime composition ancom \
  --i-table comp-gut-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --o-visualization ancom-subject.qzv


qiime taxa collapse \
  --i-table gut-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table gut-table-l6.qza

qiime composition add-pseudocount \
  --i-table gut-table-l6.qza \
  --o-composition-table comp-gut-table-l6.qza

qiime composition ancom \
  --i-table comp-gut-table-l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --o-visualization l6-ancom-subject.qzv

'''