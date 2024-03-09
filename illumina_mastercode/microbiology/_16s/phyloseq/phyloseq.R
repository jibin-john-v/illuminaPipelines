library("phyloseq")
library(glue)
library(qiime2R)
library(ggplot2)
library(cluster)
library("plyr")

##Import data ; reference https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121 # nolint
base_folder="/Users/JJOHN41/Desktop/projects/microbiome/"

features=glue("{base_folder}/Results/denoise/filtered/filtered-table-qc-passed.qza") # nolint
tree=glue("{base_folder}/Results/tree/rooted-tree.qza")
tax=glue("{base_folder}/Results/taxonomy/taxonomy.qza")
metadata=glue("{base_folder}/modified_meta_datafile.tsv")

physeq<-qza_to_phyloseq(
    features=features,
    tree=tree,tax,
    metadata=metadata) # nolint


myTaxa = names(sort(taxa_sums(physeq), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, physeq)
plot(phy_tree(ex1), show.node.label = TRUE)

svg('phylogenetic_tree.png')
plot_tree(ex1, color = "depth", label.tips = "Phylum", ladderize = "left", justify = "left" , size = "Abundance")
dev.off()


##############--------------------------------------gap statistics-----------------------------------------------------------
exord = ordinate(physeq, method="MDS", distance="jsd")

#Compute Gap Statistic
pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
x = phyloseq:::scores.pcoa(exord, display="sites")
gskmn = clusGap(x[, 1:2], FUN=pam1, K.max = 6, B = 50)

gap_statistic_ordination = function(ord, FUNcluster, type="sites", K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
    require("cluster")
    #   If "pam1" was chosen, use this internally defined call to pam
    if(FUNcluster == "pam1"){
        FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
    }
    # Use the scores function to get the ordination coordinates
    x = phyloseq:::scores.pcoa(ord, display=type)
    #   If axes not explicitly defined (NULL), then use all of them
    if(is.null(axes)){axes = 1:ncol(x)}
    #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
    clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}

#Plot Results
plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
    require("ggplot2")
    gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
    p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
    p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
    p = p + ggtitle(title)
    return(p)
}

gs = gap_statistic_ordination(exord, "pam1", B=50, verbose=FALSE)
print(gs, method="Tibs2001SEmax")

plot_clusgap(gs)

#Base graphics plotting, for comparison.
plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("Looks like 4 clusters is best, with 3 and 5 close runners up.")  


##############--------------------------------------gap statistics-----------------------------------------------------------
enterotype <- subset_taxa(physeq, Genus != "-1")

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
# This is the user-defined method:
dist_methods["designdist"]


# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]


plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
    # Calculate distance matrix
    iDist <- distance(enterotype, method=i)
    # Calculate ordination
    iMDS  <- ordinate(enterotype, "MDS", distance=iDist)
    ## Make plot
    # Don't carry over previous plot (if error, p will be blank)
    p <- NULL
    # Create plot, store as temp variable, p
    p <- plot_ordination(enterotype, iMDS, shape="depth",color="depth") # shape="Enterotype",color="depth"
    # Add title to each plot
    p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
    # Save the graphic to file.
    plist[[i]] = p
}