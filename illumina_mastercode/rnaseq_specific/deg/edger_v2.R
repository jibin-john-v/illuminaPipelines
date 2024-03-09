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
library(plyr)
library(data.table)

# Function to check if a file exists
check_input <- function(arg_value = NULL, arg_name = NULL) {
    # Check if required argument is provided
    if (!is.null(arg_value) && is.null(arg_name)) {
        stop("Error: Argument name is missing.")}
    
    if (is.null(arg_value) && !is.null(arg_name)) {
        stop(paste("Error: Argument '--", arg_name, "' is missing"))}

    # Check if file exists
    if (!is.null(arg_value)) {
        if (!file.exists(arg_value)) {
          stop(paste("Error: File not found -", arg_value))
        }}
    }


# Command line argument parsing
parser <- ArgumentParser()

parser$add_argument("--readCount", help="Provide read count file, first column should be Geneid/Gene name and rest of the column should be read count of each samples")
parser$add_argument("--sampleInfo", help="Provide file with sample info; headers [id,condition,covariates1,2 etc]; Should be tab separated")
parser$add_argument("--contrast", help="Provide file with contrast info; headers [Reference,Contrast]; Should be tab separated")
parser$add_argument("--outputName", help="Provide OutputFile/Folder Prefix")
parser$add_argument("--gtffile", help="Provide the gtf file, that used in the gene quantification")
parser$add_argument("--featuretype", help="Feature type in the gtf file used to quantify")
parser$add_argument("--gtf_atribute", help="Name of the atribute from the gtf file used to quantify")

args <- parser$parse_args()

# Check if required arguments are provided and the files exist or not
check_input(args$readCount, "readCount")
check_input(args$sampleInfo, "sampleInfo")
check_input(args$contrast, "contrast")
check_input(args$outputName, "outputName")
check_input(args$gtffile, "gtffile")
check_input(args$featuretype, "featuretype")


#Loading all files
countdata <-fread(args$ReadCount, header = T, row.names = NULL, check.names = F)
metadata <- fread(args$SampleInfo, header = T,sep = ",", row.names = NULL)
contrast <- fread(args$contrast,delim = ",",show_col_types = FALSE)
gtffile<-rtracklayer::import(args$gtffile) ## rtracklayer using to import
feature_type<-args$featuretype
outputName<-args$outputName
outputfolder<-outputName
gtf_atribute<-args$gtf_atribute

countdata <-fread('/Users/JJOHN41/Desktop/projects/bincy_company/mirna_data/gencode_miRNA_unstrand.tsv')
metadata <-fread("/Users/JJOHN41/Desktop/projects/bincy_company/SampleInfo.txt")
contrast_df <-fread("/Users/JJOHN41/Desktop/projects/bincy_company/Contraxt.txt")
gtffile<-rtracklayer::import("/Users/JJOHN41/Desktop/projects/bincy_company/mirna_data/gtf_files/miRNA_gencode.v45.chr_patch_hapl_scaff.annotation.gff3")
feature_type<-"transcript"
outputName<-"gencode"
outputfolder<-outputName
gtf_atribute<-"ID"



##Convert all data table data frame 
countdata<-as.data.frame(countdata)
geneid_column <- colnames(countdata)[1] 
row.names(countdata) <- countdata[, geneid_column]
countdata <- subset(countdata, select = -get(geneid_column))
metadata<-as.data.frame(metadata)
contrast_df<-as.data.frame(contrast_df)


##create gene details with dataframe
gtf_to_df <- function(gtf,feature_type){ 
        gtf_df<-as.data.frame(gtf)
        gtf_gene_df<-gtf_df[gtf_df$type==feature_type,]

        commands <- list(
        expression(gtf_df <- gtf_df[, c("ID","gene_id", "seqnames", "start", "end", "width", "strand", "gene_type", "gene_name")]),
        expression(gtf_df <- gtf_df[, c("gene_id", "seqnames", "start", "end", "width", "strand", "gene_type", "gene_name")]),
        expression(gtf_df <- gtf_df[, c("ID", "seqnames", "start", "end", "width", "strand", "Name")]),
        expression(gtf_df <- gtf_df[, c("ID", "seqnames", "start", "end", "width", "strand")]))

        # Attempt each command sequentially
        for (cmd in commands) {
        tryCatch({
            # Execute the command
            eval(cmd)
            # If execution succeeds, break out of the loop
            break
        }, error = function(e) {
            # If an error occurs, do nothing and continue to the next command
        })}
        return (gtf_df)
        }

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
        return(list(countdata = countdata, metadata = metadata))
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
create_pca_mds_plots <- function(y, metadata, design, plot_type = c("MDS", "PCA"),outputfolder=outputName) {
    dir.create(file.path(outputfolder, "summary"), recursive = TRUE, showWarnings = FALSE)
    out_path<-glue("{outputfolder}/summary/")

    if ("MDS" %in% plot_type) {
      mds <- plotMDS(y)
      covariates <- subset(colnames(metadata), colnames(metadata) != "id" & colnames(metadata) != "condition")

      # Check if covariates is not empty
      if (length(covariates) > 0) {
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
          ggsave(glue("{out_path}MDS_{covariate}.pdf"), plot = last_plot(), device = "pdf")

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
        ggsave("{out_path}MDS_no_covariates.pdf", plot = last_plot(), device = "pdf")
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
          ggsave(glue("{out_path}PCA_{covariate}.pdf"), plot = last_plot(), device = "pdf")

          # Print a message or do additional processing if needed
          print(glue("PCA Plot with covariate '{covariate}' saved."))

        }


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
        ggsave("{out_path}PCA.pdf", plot = last_plot(), device = "pdf")
        print("PCA Plot without covariates saved.")
      }
    }
}


#Heatmap function
create_heatmap <- function(n_genes, input_matrix, metadata,outputfolder=outputName) {
    dir.create(file.path(outputfolder, "summary"), recursive = TRUE, showWarnings = FALSE)
    out_path<-glue("{outputfolder}/summary/")

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
      colnames(ann) <- colnames(metadata2)[-1]
      colAnn <- HeatmapAnnotation(df = ann)
      # Set up pdf file
      pdf_file <- glue("{out_path}Top_{n_genes}_Variable_Genes_heatmap_raw_scaled.pdf")
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



mds_volcano_plots <- function(res_df2, outputfolder, out_prefix, contrast_name) {
    # MD plots
    pdf(glue("{outputfolder}/DES/plots/{out_prefix}_{contrast_name}_md_plots.pdf"), height = 8, width = 10)
    plotMD(res, status = is_res)
    dev.off()
    
    # Check if 'gene_name' column is not present in the dataframe
    if (!"gene_name" %in% colnames(res_df2)) {
        res_df2$gene_name <- unique(res_df2$ID)
    }
    
    # volcano_plots
    first_10_rows <- res_df2 %>%
        filter(category %in% c("Up regulated", "Down regulated")) %>%
        arrange(FDR) %>%
        head(10)
    pvalue_cutof <- first_10_rows$FDR[min(10, length(first_10_rows$FDR))]
    
    pdf(glue("{outputfolder}/DES/plots/{out_prefix}_{contrast_name}_EnhancedVolcano_volcano_plots.pdf"), height = 8, width = 10)
    print(EnhancedVolcano(res_df2, lab = res_df2$gene_name,
                    x = 'logFC', y = 'FDR',
                    pointSize = 3, labSize = 4,
                    pCutoff = pvalue_cutof, FCcutoff = 1,
                    colAlpha = 0.8,
                    labFace = 'bold',
                    boxedLabels = FALSE,
                    labCol = "red"))
    dev.off()
    
    pdf(glue("{outputfolder}/DES/plots/{out_prefix}_{contrast_name}_volcano_plots.pdf"), height = 8, width = 10)
    selected_genes <- first_10_rows[c(1, 1)]
    colnames(selected_genes) <- c("ID", "Label")
    res_df3 <- merge(res_df2, selected_genes, by = "ID", all.x = TRUE)
    
    print(ggplot(res_df3, aes(x = logFC, y = -log10(FDR), label = Label, color = category)) +
        geom_point(alpha = 0.8) +
        geom_text(col = "black", size = 2, fontface = 'bold'))
    dev.off()
    }


generate_heatmap <-function(res_df,n_genes,log2_cpm_counts,metadata,outputfolder,out_prefix,contrast_name){
    filtered_df <- res_df[res_df$category %in% c("Up regulated", "Down regulated"), ]
    filtered_df <- head(filtered_df[order(filtered_df$FDR), ],n=n_genes)
    selected_rows <- log2_cpm_counts[rownames(log2_cpm_counts) %in% filtered_df$gene_id, ]
    selected_rows<-selected_rows[,metadata$id]

    # Order metadata by condition
    metadata2 <- metadata[order(metadata$condition), ]
    # Subset highly_variable_lcpm by ordered metadata
    selected_rows <- selected_rows[, metadata2$id]

        # Check if 'id' column exists in highly_variable_lcpm
        ann <- data.frame(subset(metadata2, select = -id))
        colnames(ann) <- colnames(metadata2)[-1]
        colAnn <- HeatmapAnnotation(df = ann)
        # Set up pdf file
        pdf_file <- glue("{outputfolder}/DES/plots/{out_prefix}_{contrast_name}_Top50_Genes_Heatmap_plott.pdf")
        pdf(pdf_file, height = 10, width = 12)

        # Generate heatmap
        print(Heatmap(
            t(scale(t(selected_rows))),
            name = "expression",
            top_annotation = colAnn,
            cluster_columns = FALSE,
            show_row_names = n_genes <= 50,
            show_column_names = nrow(metadata2) <= 50
        ))
    dev.off()

    }
    




result <- check_inputdata(countdata, metadata, contrast_df)
countdata <- result$countdata
metadata <- result$metadata
gtf_df<-gtf_to_df(gtffile,feature_type)
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




dir.create(file.path(outputfolder, "DES","data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outputfolder, "DES","plots"), recursive = TRUE, showWarnings = FALSE)

#Normalization for composition bias
y <- calcNormFactors(y)

#Dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)

pdf(glue("{outputfolder}/DES/plots/biological_coefficient_of_variation.pdf"),height=6,width =8)
plotBCV(y)
dev.off()

##Fit model
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

#Plot the quasi-likelihood dispersion
pdf(glue("{outputfolder}/DES/plots/quasi-likelihood_dispersion.pdf"),height=4,width =6)
plotQLDisp(fit)
dev.off()


fdr_cutof=0.05
log_fc_cutof=1
n_genes=50

for (c_index in 1:nrow(contrast_df)){
      ##Make contrast
      # Use paste to dynamically create the contrast formula
      contrast_formula <- paste0('condition',contrast_df[c_index,'Contrast'], " - ", 'condition',contrast_df[c_index,'Reference'])
      contrast_name<-paste0(contrast_df[c_index,'Contrast'], "vs", 'condition',contrast_df[c_index,'Reference'])
      #Differential expression analysis
      out_prefix<-deparse(substitute(contrast_defined))

      # Create the contrast
      contrast_defined <- makeContrasts(contrast_formula, levels = design)
      ## get the normalised reads
      log2_cpm_counts <- cpm(y,log=TRUE)

      ## Extract the normalised mean count in each group
      mean_group_cpm<-as.data.frame(cpmByGroup(y,normalized.lib.sizes = TRUE,log = FALSE, prior.count = 1))
      mean_group_cpm$gene_id<-rownames(mean_group_cpm)

      ##-----This part need to incorporate SD and SE by group wise for future relase
      cpm_df<-as.data.frame(cpm(y,normalized.lib.sizes = TRUE,log = FALSE, prior.count = 1))

      ##-----This part need to incorporate SD and SE by group wise for future relase

      res <- glmQLFTest(fit, contrast=contrast_defined)
      is_res <- decideTestsDGE(res,adjust.method="BH",lfc=log_fc_cutof,p.value=fdr_cutof)
      #summary(is_res)
      res_df <-as.data.frame(topTags(res,adjust.method = "BH",n="inf"))
      res_df$gene_id <- rownames(res_df)
      rownames(res_df) <- NULL
      res_df <- res_df %>% mutate(category = ifelse(FDR < fdr_cutof & logFC < -log_fc_cutof, "Down regulated",
                          ifelse(FDR <fdr_cutof & logFC >log_fc_cutof, "Up regulated",
                          ifelse(FDR <fdr_cutof & logFC > -log_fc_cutof & logFC < log_fc_cutof, "No change", "Not significant"))))


      res_df<- merge(res_df, mean_group_cpm, by = "gene_id", all.x = TRUE)
      #write.csv(res_df,glue("{outputfolder}/DES/data/{out_prefix}_{contrast_name}_Gene.csv"),row.names=FALSE)
      metadata2<-t(metadata)
      colnames(metadata2) <- as.character(unlist(metadata2[1, ]))
      cpm_df$gene_id<-rownames(cpm_df)
      cpm_df2<-rbind.fill(as.data.frame(metadata2),as.data.frame(cpm_df))
      cpm_df2 <- cpm_df2[, c(length(colnames(cpm_df2)), 1:(length(colnames(cpm_df2)) - 1))]
      write.csv(cpm_df2,glue("{outputfolder}/DES/data/cpm_normalised_readcount.csv"),row.names=FALSE)

          
      # Group by the "ID" column and concatenate values of other columns ; some time in GTF file same gene may located in multiple chromosome an d or start and end position
      grouped_gtf_df <- gtf_df %>%group_by(ID) %>%summarise_all(~ paste(., collapse = ";"))
      res_df2<-merge(grouped_gtf_df,res_df,by.x=gtf_atribute,by.y="gene_id",all.y=TRUE)
      write.csv(res_df2,glue("{outputfolder}/DES/data/{out_prefix}_{contrast_name}_DES.csv"),row.names=FALSE)

      mds_volcano_plots(res_df2, outputfolder, out_prefix, contrast_name)
      generate_heatmap(res_df,n_genes,log2_cpm_counts,metadata,outputfolder,out_prefix,contrast_name)

}