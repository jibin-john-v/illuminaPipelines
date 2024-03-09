# Load required libraries
library(RColorBrewer)
library(edgeR)
library(dplyr)
library(rtracklayer)
library(tibble)
library(randomcoloR)
library(ggplot2)
library(glue)
library(limma)
library(Glimma)
library(gplots)
library(data.table)
suppressMessages(library(argparse))
library(DESeq2)
suppressMessages(library(pheatmap))
library(ComplexHeatmap)
library(EnhancedVolcano)

# Function to check if a file exists
file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("Error: File not found -", file_path))
  }
}

# Function to check if required arguments are provided
check_argument <- function(arg_value, arg_name) {
  if (is.null(arg_value)) {
    stop(paste("Error: Argument '--", arg_name, "' is missing"))
  }
}

# Command line argument parsing
parser <- ArgumentParser()

parser$add_argument("--ReadCount", help="Provide read count file, it should only contain Geneid column and gene count of each samples")
parser$add_argument("--SampleInfo", help="Provide file with sample info; headers [id,condition,covariates1,2 etc]; Should be tab separated")
parser$add_argument("--contrast", help="Provide file with contrast info; headers [Reference,Contrast]; Should be tab separated")
parser$add_argument("--OutputName", help="Provide OutputFile/Folder Prefix")

args <- parser$parse_args()

# Check if required arguments are provided
check_argument(args$ReadCount, "ReadCount")
check_argument(args$SampleInfo, "SampleInfo")
check_argument(args$contrast, "contrast")
check_argument(args$OutpuName, "OutputName")

# Check if input files exist
file_exists(args$ReadCount)
file_exists(args$SampleInfo)
file_exists(args$contrast)
#Loading all files
countdata <-read.delim(args$ReadCount, header = T, row.names = NULL, check.names = F)
metadata <- read.delim(args$SampleInfo, header = T,sep = ",", row.names = NULL)
contrast <- read_delim(args$contrast,delim = ",",show_col_types = FALSE)
metadata2 <- read_delim(args$SampleInfo,delim = ",",show_col_types = FALSE)


countdata <-fread('/Users/JJOHN41/Desktop/projects/bincy_company/mirna_data/mirgenedb_miRNA_unstrand.tsv')
metadata <-fread("/Users/JJOHN41/Desktop/projects/bincy_company/SampleInfo.txt")
contrast_df <-fread("/Users/JJOHN41/Desktop/projects/bincy_company/Contraxt.txt")


##Convert all data table data frame 
countdata<-as.data.frame(countdata)
row.names(countdata)<-countdata$Geneid
countdata <- subset(countdata, select = -Geneid)
metadata<-as.data.frame(metadata)
contrast_df<-as.data.frame(contrast_df)


check_inputdata <- function(countdata, metadata, contrast) {
        # Checking the order of the sample names
        if (all(sort(colnames(countdata)) == sort(metadata$id))) {
            # Get the order of columns based on the 'id' column in metadata
            column_order <- metadata$id
            
            # Reorder the columns of countdata
            countdata <- countdata[, column_order]
            
            # Now countdata has columns reordered based on the 'id' column in metadata
        } else {
            stop("Error: Column names do not match between countdata and metadata.")
        }

        # Checking if sample names in the ReadCount file and SampleInfo files are the same or not; if not, it will stop
        if (any(colnames(countdata) != metadata$id)) { 
            stop("\n\nSample names in the ReadCount file and SampleInfo file are not the same or not in the same order; please correct..\n\n")
        }

        # Checking for the presence of the column named 'condition' in the metadata file
        if ("condition" %in% colnames(metadata)) {
            message("Column named 'condition' is present in the provided SampleInfo file.")
        } else {
            stop("Error: Column named 'condition' is not present in the SampleInfo file provided; please check column names.")
        }

        # Introducing metadata and contrast checking
        if (all(c(as.character(contrast$Reference), as.character(contrast$Contrast)) %in% metadata$condition)) { 
            message("Reference and contrast are present in the 'condition' column of the SampleInfo file provided.")
        } else {
            stop("Error: Contrast file names are not present in the 'condition' column of the SampleInfo file provided; please check contrast file and SampleInfo.")
        }
        }


## Create Design matrix
create_model_matrix <- function(metadata) {
  # Get covariates excluding "id" and "condition"
  covariates <- subset(colnames(metadata), colnames(metadata) != "id" & colnames(metadata) != "condition")

  # Create the formula string
  formula_string <- paste("~ 0 +", "condition", " + ", paste(covariates, collapse = " + "))

  # Convert the formula string to a formula object
  formula <- as.formula(formula_string)

  # Use the formula in model.matrix
  model_matrix_result <- model.matrix(formula, data = metadata)

  return(model_matrix_result)
}


#################-------------------------------------set different color for each group--------------------------------------------
create_pca_mds_plots <- function(y, metadata, design, plot_type = c("MDS", "PCA")) {
  if ("MDS" %in% plot_type) {
    mds <- plotMDS(y)
    covariates <- subset(colnames(metadata), colnames(metadata) != "id" & colnames(metadata) != "condition")

    # Check if covariates is not empty
    if (length(covariates) > 0) {
      pdf("Rawdata_MDS.pdf", height = 10, width = 12)

      for (covariate in covariates) {
        toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, colour = metadata$condition, shape = metadata[, covariate])

        ggplot(toplot, aes(Dim1, Dim2, colour = colour, shape = shape)) +
          geom_point() +
          labs(x = glue("MDS1 ({round(mds$var.explained[1], 2) * 100}%)"),
               y = glue("MDS2 ({round(mds$var.explained[2], 2) * 100}%)"),
               color = "condition",
               shape = covariate) +
          theme_minimal()  # Optional: Adjust the theme as per your preference

        # Save the current plot
        ggsave(glue("Rawdata_MDS_{covariate}.pdf"), plot = last_plot(), device = "pdf")

        # Print a message or do additional processing if needed
        print(glue("MDS Plot with covariate '{covariate}' saved."))
      }

      dev.off()

    } else {
      # If covariates is empty, create a plot without shape
      toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, colour = metadata$condition)

      ggplot(toplot, aes(Dim1, Dim2, colour = colour)) +
        geom_point() +
        labs(x = glue("MDS1 ({round(mds$var.explained[1], 2) * 100}%)"),
             y = glue("MDS2 ({round(mds$var.explained[2], 2) * 100}%)"),
             color = "condition") +
        theme_minimal()  # Optional: Adjust the theme as per your preference

      # Save the plot
      ggsave("Rawdata_MDS_no_covariates.pdf", plot = last_plot(), device = "pdf")
      print("MDS Plot without covariates saved.")
    }
  }

  if ("PCA" %in% plot_type) {
    edgeR.DDS <- DESeqDataSetFromMatrix(countData = round(y$counts), colData = y$samples, design = design)
    transform.edgeR.DDS <- rlog(edgeR.DDS, blind = TRUE)
    pcaData <- plotPCA(transform.edgeR.DDS, intgroup = c("group"), returnData = TRUE, ntop = 1000)  
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    covariates <- subset(colnames(metadata), colnames(metadata) != "id" & colnames(metadata) != "condition")

    # Check if covariates is not empty
    if (length(covariates) > 0) {
      for (covariate in covariates) {
        toplot <- data.frame(PC1 = pcaData$PC1, PC2 = pcaData$PC2, colour = metadata$condition, shape = metadata[, covariate])
        ggplot(toplot, aes(PC1, PC2, colour = colour, shape = shape)) +
          geom_point() +
          labs(x = glue("PC1 ({percentVar[1]}%)"),
               y = glue("PC2 ({percentVar[2]}%)"),
               color = "condition",
               shape = covariate) +
          theme_minimal()  # Optional: Adjust the theme as per your preference

        # Save the current plot
        ggsave(glue("Rawdata_PCA_{covariate}.pdf"), plot = last_plot(), device = "pdf")

        # Print a message or do additional processing if needed
        print(glue("PCA Plot with covariate '{covariate}' saved."))
      }

      dev.off()

    } else {
      # If covariates is empty, create a plot without shape
      toplot <- data.frame(PC1 = pcaData$PC1, PC2 = pcaData$PC2, colour = metadata$condition)
      ggplot(toplot, aes(PC1, PC2, colour = colour)) +
        geom_point() +
        labs(x = glue("PC1 ({percentVar[1]}%)"),
             y = glue("PC2 ({percentVar[2]}%)"),
             color = "condition") +
        theme_minimal()  # Optional: Adjust the theme as per your preference

      # Save the plot
      ggsave("Rawdata_PCA.pdf", plot = last_plot(), device = "pdf")
      print("PCA Plot without covariates saved.")
    }
  }
}


#Heatmap function
create_heatmap <- function(n_genes, input_matrix, metadata) {
  # Calculate variance for each gene
  var_genes <- apply(input_matrix, 1, var)

  # Get the gene names for the top N most variable genes
  select_var <- names(sort(var_genes, decreasing = TRUE))[1:n_genes]

  # Subset logcounts matrix
  highly_variable_lcpm <- input_matrix[select_var, ]

  # Order metadata by condition
  metadata2 <- metadata[order(metadata$condition), ]

  # Subset highly_variable_lcpm by ordered metadata
  highly_variable_lcpm <- highly_variable_lcpm[, metadata2$id]
  
  if (any(colnames(highly_variable_lcpm) != metadata2$id)) 
        { 
            stop("\n\nSample names in the ReadCount file and SampleInfo file are not the same or not in the same order; please correct..\n\n")
        }else{

  # Check if 'id' column exists in highly_variable_lcpm
    ann <- data.frame(subset(metadata2, select = -id))
    colnames(ann) <- c('condition', 'sample_series_id')
    colAnn <- HeatmapAnnotation(df = ann)

    # Set up pdf file
    pdf_file <- glue("Top_{n_genes}_Variable_Genes_heatmap_raw_scaled.pdf")
    pdf(pdf_file, height = 10, width = 12)

    # Generate heatmap
    print(Heatmap(
      t(scale(t(highly_variable_lcpm))),
      name = "expression",
      top_annotation = colAnn,
      cluster_columns = FALSE,
      show_row_names = n_genes <= 50,
      show_column_names = nrow(metadata2) <= 50
    ))

    # Close pdf file
    dev.off()
}}



check_inputdata(countdata, metadata, contrast_df)
design <- create_model_matrix(metadata)

## createe DGEList object
y <- DGEList(counts=countdata, group=metadata$condition,remove.zeros = TRUE)

##Remove genes with low expression
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]  #keep.lib.sizes=FALSE causes the library sizes to be recomputed after the filtering. 

create_pca_mds_plots(y, metadata, design, plot_type="MDS")
create_pca_mds_plots(y, metadata, design, plot_type="PCA")
create_heatmap(n_genes = 50, cpm(y,log=TRUE), metadata)
create_heatmap(n_genes = 25, cpm(y,log=TRUE), metadata)
create_heatmap(n_genes = 100, cpm(y,log=TRUE), metadata)
create_heatmap(n_genes = 500, cpm(y,log=TRUE), metadata)





#Normalization for composition bias
y <- calcNormFactors(y)

#Dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)



pdf("biological_oefficient_of_variation.pdf",height=6,width =8)
plotBCV(y)
dev.off()


##Fit model
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

#Plot the quasi-likelihood dispersion
pdf("quasi-likelihood_dispersion.pdf",height=4,width =6)
plotQLDisp(fit)
dev.off()


fdr_cutof=0.05
log_fc_cutof=1
##Make contrast
# Use paste to dynamically create the contrast formula
contrast_formula <- paste0('condition',contrast_df[1,'Contrast'], " - ", 'condition',contrast_df[1,'Reference'])
# Create the contrast
contrast_defined <- makeContrasts(contrast_formula, levels = design)
## get the normalised reads
log2_cpm_counts <- cpm(y,log=TRUE)

## Extract the normalised mean count in each group
mean_group_cpm<-as.data.frame(cpmByGroup(y,normalized.lib.sizes = TRUE,log = FALSE, prior.count = 1))
mean_group_cpm$gene_id<-rownames(mean_group_cpm)

#Differential expression analysis
out_prefix<-deparse(substitute(contrast_defined))


res <- glmQLFTest(fit, contrast=contrast_defined)
is_res <- decideTestsDGE(res,adjust.method="BH",lfc=log_fc_cutof,p.value=fdr_cutof)
summary(is_res)
res_df <-as.data.frame(topTags(res,adjust.method = "BH",n="inf"))
res_df$gene_id <- rownames(res_df)
rownames(res_df) <- NULL
res_df <- res_df %>% mutate(category = ifelse(FDR < fdr_cutof & logFC < -log_fc_cutof, "Down regulated",
                    ifelse(FDR <fdr_cutof & logFC >log_fc_cutof, "Up regulated",
                    ifelse(FDR <fdr_cutof & logFC > -log_fc_cutof & logFC < log_fc_cutof, "No change", "Not significant"))))


contrast_name<-paste0(contrast_df[1,'Contrast'], "vs", 'condition',contrast_df[1,'Reference'])

res_df<- merge(res_df, mean_group_cpm, by = "gene_id", all.x = TRUE)
write.csv(res_df,glue("{out_prefix}_{contrast_name}_Gene.csv"),row.names=FALSE)


## MD plots
pdf(glue("{out_prefix}_{contrast_name}_md_plots.pdf"),height=8,width =10)
plotMD(res,status=is_res)
dev.off()

## volcano_plots
pdf(glue("{out_prefix}_{contrast_name}_EnhancedVolcano_volcano_plots.pdf"),height=8,width =10)
EnhancedVolcano(res_df,lab = res_df$gene_name,
            x = 'logFC',y = 'FDR',pointSize = 1,labSize =2.5,pCutoff = 1e-06,colAlpha=0.8,labCol="blue")
dev.off()

pdf(glue("{out_prefix}_volcano_plots.pdf"),height=8,width =10)
selected_genes <- res_df %>% arrange(FDR) %>% head(10) %>% select(gene_id)
res_df <- res_df %>% mutate(Label = ifelse(gene_id %in% selected_genes$gene_id, gene_name, ""))
ggplot(res_df,aes(x = logFC, y = -log10(FDR), label=Label,color=category)) + geom_point(alpha=0.8) + geom_text(col="black",size=2,fontface='bold')
dev.off()


filtered_df <- res_df[res_df$category %in% c("Up regulated", "Down regulated"), ]
filtered_df <- head(filtered_df[order(filtered_df$FDR), ],n=50)
selected_rows <- log2_cpm_counts[rownames(log2_cpm_counts) %in% filtered_df$gene_id, ]
selected_rows<-selected_rows[,metadata$id]
identical(colnames(x), sample$SRRID)


pdf(glue("{out_prefix}_Top_Heatmap_plots.pdf"),height=8,width =10)
# Plot the heatmap
heatmap.2(selected_rows,Colv=FALSE,
      main="Top 500 most variable genes across samples",ColSideColors=gse_tcga_meta2$color,scale="row",
      dendrogram="row",trace="none", margin=c(8,9))
legend("topright", legend = unique(gse_tcga_meta2$group), fill = unique(gse_tcga_meta2$color), title = "Legend Title")
dev.off()




