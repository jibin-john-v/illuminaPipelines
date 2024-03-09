
suppressMessages(library(argparse))
suppressMessages(library(ChIPQC))

parser <- ArgumentParser()
parser$add_argument("--metadata",
    help='csv file with meta data')
#it should contain following columns, "SampleID","Factor","Replicate","bamReads","ControlID","bamControl","Peaks","PeakCaller","Tissue","Condition"]'
args <- parser$parse_args()

samples<-read.csv(args$metadata)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
chipObj <- ChIPQC(samples, annotation="hg38") 
ChIPQCreport(chipObj, reportName="ChIP_QC_report", reportFolder="Results/ChIPQCreport")

