
import os 
#rMATS ; https://github.com/Xinglab/rmats-turbo ;https://rnaseq-mats.sourceforge.io/


def rmats_run(g1_bamfiles,g2_bamfiles,gtf_file,librarytype,readlength,threads):
    rmat_dir=f"Results/splice_variant/rmats/"
    os.system(f"mkdir -p {rmat_dir}")
    os.system(f'''python /mnt/disks/sdc/Software/rmats_turbo_v4_2_0/rmats.py \
                        --b1 {g1_bamfiles} \
                        --b2 {g2_bamfiles} \
                        --gtf {gtf_file} \
                        -t {librarytype} \
                        --readLength {readlength} \
                        --variable-read-length \
                        --allow-clipping \
                        --nthread {threads} \
                        --od {rmat_dir} \
                        --tmp {rmat_dir} \
                        --novelSS ''')


#MAJIQ: Modeling Alternative Junction Inclusion Quantification ; https://majiq.biociphers.org/
#Bisbee    :https://www.nature.com/articles/s41598-021-89938-2 
#SpliceVault : https://www.nature.com/articles/s41588-022-01293-8 


#SpliceTools ; https://github.com/flemingtonlab/SpliceTools

#AS-Quant ; https://github.com/compbiolabucf/AS-Quant/tree/main

#SUPPA2 : https://github.com/comprna/SUPPA

#betAS ;https://github.com/DiseaseTranscriptomicsLab/betAS/tree/dev?tab=readme-ov-file

#fraser  ;https://github.com/gagneurlab/fraser

#LeafCutter https://davidaknowles.github.io/leafcutter/index.html
    
#SpliceWiz ; https://github.com/alexchwong/SpliceWiz?tab=readme-ov-file

#ASTool ;http://zzdlab.com/ASTool/index.php

#Whippet ; https://github.com/timbitz/Whippet.jl

#dasper https://dzhang32.github.io/dasper/articles/dasper.html

#MISO ; https://miso.readthedocs.io/en/fastmiso/index.html#how-miso-works

#VAST-TOOLS ; https://github.com/vastgroup/vast-tools