#### References

1. Sanger Institute  		: https://www.singlecellcourse.org/index.html
2. holab                  		:https://holab-hku.github.io/Fundamental-scRNA/
3. hbctraining         		: https://github.com/hbctraining/scRNA-seq_online/tree/master
4. Cambridge         		: https://bioinformatics-core-shared-training.github.io/SingleCell_RNASeq_Sept23/
5. babraham           		: https://www.bioinformatics.babraham.ac.uk/training.html#10xrna
6. bioconductor			:https://bioconductor.org/books/3.12/OSCA/
7. SingleCellMultiOmics 	:https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/
8. sib-swiss				:https://sib-swiss.github.io/single-cell-training-archived/2023.3/
9. erilu					:https://github.com/erilu/single-cell-rnaseq-analysis/tree/master

#### [METHODS](https://bioinformatics-core-shared-training.github.io/SingleCell_RNASeq_Sept23/UnivCambridge_ScRnaSeqIntro_Base/Slides/01_Introduction.pdf "METHODS")

* Cells in wells, traps and valves (nanowell, Flow sorting, CellenOne, FluidigmC1, SmartSeq, plexWell/seqWell)

  * Screen for and retrieve single cells of interest-Enrich for rare cells with decided properties
  * Control the cellular microenvironment-Monitor and control cell-cell interactions
  * Precise/extensive manipulation of single cells
* Droplets (Drop-seq, 10x Genomics)

  * Introduce distinct ‘packets’ of reagents to single cell (e.g. barcodes)
  * Perform amplification on individual cells
  * Sort large population of single cells
* Combinatorial indexing (SCI-seq, SPLiT-seq)

  * Economic use of reagents for cell separation
  * Efficiency of handling larger population than Drop-seq-Maintain
  * complexities of population without bias from droplet or well
* ‘One tube’ solutions (SciPioBioscience, Fluent  Bioscience, CS Genetics)
* [for details](https://www.annualreviews.org/doi/10.1146/annurev-biodatasci-080917-013452#)

#### A classical scRNA-seq workflow contains four main steps:

1. ###### Mapping the cDNA fragments to a reference;

   [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) is the default tool for processing 10x Genomics Chromium scRNAseq data. It uses STAR aligner, which performs splicing-aware alignment of reads to the genome. After this, it uses the transcript annotation GTF to bucket the reads into exonic, intronic, and intergenic, and by whether the reads align (confidently) to the genome.

   [Primary genome assembly](https://www.ncbi.nlm.nih.gov/grc/help/definitions) versions (i.e. without the ALT loci) are used for the alignment in all versions. The following biotypes are retained: protein coding, long noncoding RNA, antisense, and all biotypes belonging to BCR/TCR (i.e. V/D/J) genes.All pseudogenes and small noncoding RNAs are removed.
2. ##### Assigning reads to genes;
3. Assigning reads to cells (cell barcode demultiplexing);
4. Counting the number of unique RNA molecules (UMI deduplication).

##### QC Filtering

1. **Aim of QC is …**
   1. To remove undetected genes
   2. To remove empty droplets
   3. To remove droplets with dead cells
   4. To remove Doublet/multiplet
   5. Ultimately To filter the data to only include true cells that are of high quality

* The number of unique genes detected in each cell.
  * Low-quality cells or empty droplets will often have very few genes
  * Cell doublets or multiplets may exhibit an aberrantly high gene count
* Total number of molecules detected within a cell (correlates strongly with unique genes)
* The percentage of reads that map to the mitochondrial genome
  * Low-quality / dying cells often exhibit extensive mitochondrial contamination
  * We use the set of all genes starting with `MT-` as a set of mitochondrial genes
