library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)


pheno_data = read.csv("/Users/JJOHN41/Desktop/projects/test_data/rnaseq/chrX_data/geuvadis_phenodata.csv")


bg_chrX = ballgown(dataDir = "/Users/JJOHN41/Desktop/projects/illumina/pipe_lines/rnaseq/Results/ballgown/", 
                    samplePattern = "ERR", 
                    pData=pheno_data)

#Filter to remove low abundance genes.
bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) > 1",genomesubset=TRUE)


#Identify transcripts that show statistically significant differences between groups. 
 results_transcripts = stattest(bg_chrX_filt, 
                            feature="transcript", 
                            covariate="sex", 
                            adjustvars = c("population"), 
                            getFC=TRUE, meas="FPKM")



 results_genes = stattest(bg_chrX_filt, 
                            feature="gene", 
                            covariate="sex", 
                            adjustvars = c("population"), 
                            getFC=TRUE, meas="FPKM")


#Add gene names and gene IDs to the results_transcripts data frame:
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), 
                                geneIDs=ballgown::geneIDs(bg_chrX_filt), 
                                results_transcripts)


##Sort the results from smallest p-value to largest: > 
results_transcripts = arrange(results_transcripts,pval) 
results_genes = arrange(results_genes,pval)


##Write the results to a CSV 
write.csv(results_transcripts, "chrX_transcript_results.csv", row.names=FALSE) 
write.csv(results_genes, "chrX_gene_results.csv", row.names=FALSE)