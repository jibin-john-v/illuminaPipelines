##This is for macbook ; mamba activate methylseq


import sys,os,glob,subprocess
from pathlib import Path
import pandas as pd
import concurrent.futures


sys.path.insert(0, '../../master_codes/') 
from fastq_qc.fastq_qc import falco_run,fastp_run



##https://github.com/yupenghe/methylpy/tree/methylpy
##https://bioinformaticsdotca.github.io/EPI_2021_Module4_lab
##https://github.com/NIGMS/DNA-Methylation-Sequencing-Analysis-with-WGBS/blob/main/tutorial_2-metilene.ipynb
##https://github.com/sagc-bioinformatics/SAGC_5mCMethylation_tutorial/blob/main/tutorial/04_methylation_calling.md
##https://www.atsjournals.org/doi/10.1165/rcmb.2019-0150TR
##https://qcb.ucla.edu/collaboratory/resources/
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6775954/
##https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq_analysis.html
##https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/WGBS/WGBS
##https://smithlabresearch.org/software/methpipe/
##https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-018-0194-0
##https://as-botanicalstudies.springeropen.com/articles/10.1186/s40529-022-00366-5
##http://statgen.cnag.cat/GEMBS/UserGuide/_build/html/index.html
##https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html
