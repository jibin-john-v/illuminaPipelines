library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(sf)
library(glue)
library(celldex)
library(SingleR)
library(singleCellTK)
library(SingleCellExperiment)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)


file_path="/Users/JJOHN41/Desktop/projects/test_data/single_cell_rnaseq/data/"
project_name<-"pbmc3k"
out_path="/Users/JJOHN41/Desktop/projects/test_analysis/single_cell/"

files <- files <- list.files(file_path, recursive = FALSE, full.names = FALSE)

##Identify all the samples in the folder
files <- files[grep("\\matrix$", files)]
dir.create(glue("{out_path}data"), showWarnings = FALSE)
dir.create(glue("{out_path}Plots"), showWarnings = FALSE)


create_seurat_object <- function(files, file_path) {
        cell_ids <- character(0)  # Initialize an empty character vector
        seurat_obj_list <- list()  # Initialize an empty list to store Seurat objects

        # Create a Seurat object for each sample
        for (file in files) {
                seurat_data <- Read10X(data.dir = file.path(file_path, file))
                seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                                min.features = 100, 
                                                project = file)
                

                # Extract counts
                counts <- GetAssayData(object = seurat_obj, slot = "counts")
                # Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
                nonzero <- counts > 0
                # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
                keep_genes <- Matrix::rowSums(nonzero) >= 10
                # Only keeping those genes expressed in more than 10 cells
                filtered_counts <- counts[keep_genes, ]
                # Reassign to filtered Seurat object
                seurat_obj <- CreateSeuratObject(filtered_counts, meta.data = seurat_obj@meta.data)
                assign(file, seurat_obj)

                # Extract cell ID from the file name
                cell_id <- strsplit(file, "_")[[1]][1]

                # Append the cell ID to the vector
                cell_ids <- append(cell_ids, cell_id)

                # Append the Seurat object to the list
                seurat_obj_list <- append(seurat_obj_list, list(seurat_obj))
                }

        # Merge Seurat objects
        merged_seurat <- merge(x = seurat_obj_list[[1]], 
                                y = seurat_obj_list[-1], 
                                add.cell.id = cell_ids)

        return(list(seurat_obj = merged_seurat, cell_ids = cell_ids))

        }


##calculate additional meta data information
calculate_metrics <- function(files,cell_ids,seurat_object,out_path) {
        # Novelty score: Add number of genes per UMI for each cell to metadata
        seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)

        # Mitochondrial Ratio
        seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
        seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100

        # Ribosomal RNA percentage
        seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")

        # Ig genes percentage
        seurat_object[["percent.igg"]] <- PercentageFeatureSet(seurat_object, pattern = "^HLA-|IG[HJKL]")

        # Create metadata dataframe
        metadata <- seurat_object@meta.data
        # Add cell IDs to metadata
        metadata$cells <- rownames(metadata)

        # Create sample column
        metadata$sample <- NA
        for (index in seq(length(files) )) {
            metadata$sample[which(str_detect(metadata$orig.ident, files[index]))] <- cell_ids[index]
            }
        
        # Rename columns
        metadata <- tryCatch({
                metadata %>%
                        dplyr::rename(seq_folder = orig.ident,
                                nUMI = nCount_RNA,
                                nGene = nFeature_RNA)

                }, error = function(e) {
                # Handle the error
                print(paste("An error occurred:", e$message))
                # Return a default value or do some recovery action
                return(NA)
                })


        # Add metadata back to Seurat object
        #seurat_object@meta.data <- metadata
                                
        # Create .RData object to load at any time
        #save(seurat_object, file=glue("{out_path}data/merged_filtered_seurat.RData"))
        return(metadata=metadata)

        }



visualize_seurat_data <- function(metadata, out_path, prefix) {
        # Visualize the number of cell counts per sample
        png(file = file.path(glue("{out_path}/Plots/{prefix}_number_of_cell_counts_per_sample.png")), 
                                units = "in", width = 12, height = 12, res = 800) 
            print(metadata %>% 
            ggplot(aes(x=sample, fill=sample)) + 
            geom_bar() +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            theme(plot.title = element_text(hjust=0.5, face="bold")) +
            ggtitle("NCells"))
        dev.off()

        # Visualize the number UMIs/transcripts per cell
        png(file = file.path(glue("{out_path}/Plots/{prefix}_number_of_UMIs_counts_per_transcripts_percell_sample.png")), 
        units = "in", width = 12, height = 12, res = 800) 
            print(metadata %>% 
            ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
            geom_density(alpha = 0.2) + 
            scale_x_log10() + 
            theme_classic() +
            ylab("Cell density") +
            geom_vline(xintercept = 500))
        dev.off()

        # Visualize the distribution of genes detected per cell via histogram
        png(file = file.path(glue("{out_path}/Plots/{prefix}_distribution_of_genes_detected_per_cell.png")), 
        units = "in", width = 12, height = 12, res = 800) 
            print(metadata %>% 
            ggplot(aes(color=sample, x=nGene, fill= sample)) + 
            geom_density(alpha = 0.2) + 
            theme_classic() +
            scale_x_log10() + 
            geom_vline(xintercept = 300))
        dev.off()

        # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
        png(file = file.path(glue("{out_path}/Plots/{prefix}_genes_detected_per_UMI.png")), 
        units = "in", width = 12, height = 12, res = 800) 
            print(metadata %>%
            ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
            geom_density(alpha = 0.2) +
            theme_classic() +
            geom_vline(xintercept = 0.8))
        dev.off()

        # Visualize the distribution of mitochondrial gene expression detected per cell
        png(file = file.path(glue("{out_path}/Plots/{prefix}_Mitochondrial_counts_ratio.png")), 
        units = "in", width = 12, height = 12, res = 800) 
            print(metadata %>% 
            ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
            geom_density(alpha = 0.2) + 
            scale_x_log10() + 
            theme_classic() +
            geom_vline(xintercept = 0.2))
        dev.off()

        # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
        png(file = file.path(glue("{out_path}/Plots/{prefix}_correlation_between_genes_detected_and_number_of_UMIs.png")), 
        units = "in", width = 12, height = 12, res = 800) 
            print(metadata %>% 
            ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
            geom_point() + 
            scale_colour_gradient(low = "gray90", high = "black") +
            stat_smooth(method=lm) +
            scale_x_log10() + 
            scale_y_log10() + 
            theme_classic() +
            geom_vline(xintercept = 500) +
            geom_hline(yintercept = 250) +
            facet_wrap(~sample))
        dev.off()

}


create_plots <- function(seurat_obj, out_path, prefix, feature_pairs) {
        # Create directory if it doesn't exist # nolint
        if (!dir.exists(glue("{out_path}/Plots"))) 
                            dir.create(glue("{out_path}/Plots"), recursive = TRUE) # nolint

        # Save violin plot as PNG # nolint
        png(file = file.path(glue("{out_path}/Plots/{prefix}_Violinplots.png")), units = "in", width = 12, height = 12, res = 800) 
        print(VlnPlot(seurat_obj, 
                features = c("nGene", "nUMI", "mitoRatio", "percent.ribo","log10GenesPerUMI"),
                ncol = 5, pt.size = 0.1) + theme(plot.title = element_text(size = 10))) 
        dev.off()

        for (pair in feature_pairs) {
                filename <- glue("{out_path}/Plots/{prefix}_{pair[1]}_vs_{pair[2]}_scatter.png") 
                # Create scatter plot and save as PNG 
                png(file = filename, units = "in", width = 12, height = 12, res = 800) 
                print(FeatureScatter(seurat_obj, feature1 = pair[1], feature2 = pair[2])) 
                dev.off()
         }
        } 


# Assuming files is a vector of file names and file_path is the path to the directory containing the files
seurat_obj <- create_seurat_object(files, file_path)
merged_seurat<-seurat_obj$seurat_obj
cell_ids<-seurat_obj$cell_ids

# calculate_metrics
metadata<- calculate_metrics(files,cell_ids,merged_seurat,out_path) 
merged_seurat@meta.data<-metadata

save(merged_seurat, file=glue("{out_path}data/merged_filtered_seurat.RData"))

# visualize_seurat_data
visualize_seurat_data(metadata, out_path, "Rawdata")

feature_pairs <- list(
        c("nUMI", "mitoRatio"),c("nUMI", "nGene"),c("nUMI", "percent.ribo"),
        c("percent.ribo", "mitoRatio"),c("nGene","mitoRatio"))

create_plots(merged_seurat, out_path, "Rawdata", feature_pairs)


# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
			                          

# visualize filtered seurat_data
metadata2<- as.data.frame(filtered_seurat@meta.data)
visualize_seurat_data(metadata2, out_path, "Filtered")
feature_pairs <- list(
        c("nUMI", "mitoRatio"),c("nUMI", "nGene"),c("nUMI", "percent.ribo"),
        c("percent.ribo", "mitoRatio"),c("nGene","mitoRatio"))

create_plots(filtered_seurat, out_path, "Filtered", feature_pairs)
save(filtered_seurat, file=glue("{out_path}data/merged_filtered_QC_seurat.RData"))



# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

seurat_phase<- JoinLayers(seurat_phase)

# Score cells for cell cycle
# segregate this list into markers of G2/M phase and markers of S phase
s_genes <-cc.genes.updated.2019$s.genes
g2m_genes <-cc.genes.updated.2019$g2m.genes
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

ggplot(aes(Phase)) + geom_bar()


ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))

  
#Apply sctransform normalization https://satijalab.org/seurat/articles/sctransform_vignette ;single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
seurat_obj_filt <- SCTransform(seurat_phase,
                               vars.to.regress = c('mitoRatio','S.Score', 'G2M.Score','nUMI','nGene','log10GenesPerUMI'),
                                verbose = FALSE)






##Assign Cell-Cycle Scores




#add the doublet annotation generated by scrublet
#doublets <- read.table(glue('{file_path}scrublet_calls.tsv'),header = T,row.names = 1)
#seurat_obj <- AddMetaData(seurat_obj,doublets)
#head(seurat_obj[[]])









create_plots(merged_seurat, out_path, "Rawdata", feature_pairs)


minimum_nfeature=200
maximum_nfeature=2000
percent_mt=5
seurat_obj_filt <- subset(seurat_obj, subset = nFeature_RNA > minimum_nfeature & nFeature_RNA < maximum_nfeature & percent.mt < percent_mt)
create_plots(seurat_obj_filt, "Results", "Filtered", feature_pairs)


##Duplite base filtering
seurat_obj[['QC']] <- ifelse(seurat_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
seurat_obj[['QC']] <- ifelse(seurat_obj@meta.data$nFeature_RNA < 500 & seurat_obj@meta.data$QC == 'Pass','Low_nFeature',seurat_obj@meta.data$QC)
seurat_obj[['QC']] <- ifelse(seurat_obj@meta.data$nFeature_RNA < 500 & seurat_obj@meta.data$QC != 'Pass' & seurat_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seurat_obj@meta.data$QC,sep = ','),seurat_obj@meta.data$QC)
seurat_obj[['QC']] <- ifelse(seurat_obj@meta.data$percent.mt > 15 & seurat_obj@meta.data$QC == 'Pass','High_MT',seurat_obj@meta.data$QC)
seurat_obj[['QC']] <- ifelse(seurat_obj@meta.data$nFeature_RNA < 500 & seurat_obj@meta.data$QC != 'Pass' & seurat_obj@meta.data$QC != 'High_MT',paste('High_MT',seurat_obj@meta.data$QC,sep = ','),seurat_obj@meta.data$QC)
table(seurat_obj[['QC']])



#Apply sctransform normalization https://satijalab.org/seurat/articles/sctransform_vignette ;single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().

seurat_obj_filt <- SCTransform(seurat_obj_filt, vars.to.regress = "percent.mt", verbose = FALSE)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj_filt), 10)

# plot variable features with and without labels
png(file = file.path(glue("{out_path}/Plots/Filtereddata_Top10_variablefeatres.png")), units = "in", width = 12, height = 12, res = 800) 
plot1 <- VariableFeaturePlot(seurat_obj_filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

#linear dimensional reduction
seurat_obj_filt <- RunPCA(seurat_obj_filt, verbose = FALSE)

print(seurat_obj_filt[["pca"]], dims = 1:5, nfeatures = 5)

png(file = file.path(glue("{out_path}/Plots/Filtereddata_Genes_in_PC1_PC2.png")), units = "in", width = 12, height = 12, res = 800) 
VizDimLoadings(seurat_obj_filt, dims = 1:2, reduction = "pca")
dev.off()

png(file = file.path(glue("{out_path}/Plots/Filtereddata_PCAplot_PC1_PC2.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, reduction = "pca") + NoLegend()
dev.off()

png(file = file.path(glue("{out_path}/Plots/Filtereddata_DimHeatmap_PC1_PC2.png")), units = "in", width = 12, height = 12, res = 800) 
DimHeatmap(seurat_obj_filt, dims = 1, cells = 500, balanced = TRUE)
dev.off()

png(file = file.path(glue("{out_path}/Plots/Filtereddata_DimHeatmap_PC1_PC15.png")), units = "in", width = 12, height = 12, res = 800) 
DimHeatmap(seurat_obj_filt, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

png(file = file.path(glue("{out_path}/Plots/Filtereddata_ElbowPlot.png")), units = "in", width = 12, height = 12, res = 800) 
ElbowPlot(seurat_obj_filt,ndims = 50)
dev.off()


#non-linear dimensional reduction (UMAP/tSNE)
seurat_obj_filt <- RunUMAP(seurat_obj_filt, dims = 1:30, verbose = FALSE)
seurat_obj_filt <- FindNeighbors(seurat_obj_filt, dims = 1:30, verbose = FALSE)
seurat_obj_filt <- FindClusters(seurat_obj_filt, verbose = FALSE)
png(file = file.path(glue("{out_path}/Plots/Filtereddata_UMAP.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, label = TRUE,reduction = "umap")
dev.off()


##SNE
seurat_obj_filt <-RunTSNE(seurat_obj_filt, dims = 1:30, verbose = FALSE)
png(file = file.path(glue("{out_path}/Plots/Filtereddata_TNSE.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, label = TRUE,reduction = "tsne")
dev.off()

#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurat_obj_filt, ident.1 = 2)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(seurat_obj_filt, ident.1 = 0, ident.2 = c(1, 19))
head(cluster5.markers, n = 5)



###-----------
monaco.ref <- celldex::MonacoImmuneData()
sce <- as.SingleCellExperiment(DietSeurat(seurat_obj_filt))
sce

monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)

table(monaco.main$pruned.labels)

table(monaco.fine$pruned.labels)

seurat_obj_filt@meta.data$monaco.main <- monaco.main$pruned.labels
seurat_obj_filt@meta.data$monaco.fine <- monaco.fine$pruned.labels

seurat_obj_filt <- SetIdent(seurat_obj_filt, value = "monaco.fine")

png(file = file.path(glue("{out_path}/Plots/Filtereddata_TNSE_with_genemarkers.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, label = T , repel = T, label.size = 3,reduction = "tsne") + NoLegend()
dev.off()

png(file = file.path(glue("{out_path}/Plots/Filtereddata_UMAP_with_genemarkers.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, label = T , repel = T, label.size = 3,reduction = "umap") + NoLegend()
dev.off()

seurat_obj_filt <- FindNeighbors(seurat_obj_filt, dims = 1:10)
FindClusters(seurat_obj_filt)
head(Idents(seurat_obj_filt), 5)




png(file = file.path(glue("{out_path}/Plots/Filtereddata_DimHeatmap_tsne_PC1_PC2.png")), units = "in", width = 12, height = 12, res = 800) 
DimHeatmap(seurat_obj_filt, dims = 1, cells = 500, balanced = TRUE,reduction = 'umap')
dev.off()


seurat_obj_filt <- RunTSNE(seurat_obj_filt, dims = 1:50, verbose = FALSE)

png(file = file.path(glue("{out_path}/Plots/Filtereddata_Dimplot_50tsne_PC1_PC2.png")), units = "in", width = 12, height = 12, res = 800) 
DimHeatmap(seurat_obj_filt, reduction = "tsne")
dev.off()




# These are now standard steps in the Seurat workflow for visualization and clustering
seurat_obj_filt <- RunPCA(seurat_obj_filt, verbose = FALSE)
seurat_obj_filt <- RunUMAP(seurat_obj_filt, dims = 1:30, verbose = FALSE)

seurat_obj_filt <- FindNeighbors(seurat_obj_filt, dims = 1:30, verbose = FALSE)
seurat_obj_filt <- FindClusters(seurat_obj_filt, verbose = FALSE)
png(file = file.path(glue("{out_path}/Plots/Filtereddata_Dimplot_50tsne_PC1_PC2.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, label = TRUE)
dev.off()



# SCTranform

seurat_obj_filt <- RunPCA(seurat_obj_filt, verbose = FALSE)
seurat_obj_filt <- RunUMAP(seurat_obj_filt, dims = 1:30, verbose = FALSE)

seurat_obj_filt <- FindNeighbors(seurat_obj_filt, dims = 1:30, verbose = FALSE)
seurat_obj_filt <- FindClusters(object = seurat_obj_filt,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4,1.8))

Idents(object = seurat_obj_filt) <- "integrated_snn_res.1.8"
seurat_integrated <- RunUMAP(seurat_obj_filt, 
                 reduction = "pca", 
                 dims = 1:40)


png(file = file.path(glue("{out_path}/Plots/Filtereddata_Dimplot_50tsne_PC1_PC2_.1.8.png")), units = "in", width = 12, height = 12, res = 800) 
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()

png(file = file.path(glue("{out_path}/Plots/Filtereddata_Dimplot_50tsne_PC1_PC2_.1.8_2.png")), units = "in", width = 12, height = 12, res = 800) 
DimPlot(seurat_obj_filt, label = T , repel = T, label.size = 3,reduction = "umap") + NoLegend()
dev.off()