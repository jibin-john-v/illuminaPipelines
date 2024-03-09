

##gff to bed 
gffread \
    --bed gencode.v45.chr_patch_hapl_scaff.annotation.gff3 \
     > gencode.v45.chr_patch_hapl_scaff.annotation.bed


#https://rnabio.org/module-02-alignment/0002/06/01/Alignment_QC/
gff2bed < ref_ribosome.gtf > ref_ribosome.bed


# Create a genePred file for our reference transcriptome
gtfToGenePred \
    -genePredExt chr22_with_ERCC92.gtf chr22_with_ERCC92.ref_flat.txt


genePredToBed \
    chr22_with_ERCC92.genePred \
    chr22_with_ERCC92.bed12



##convert exome bed file for qualimap

awk -F"\t" '{print $1,$2,$3,"."}' S04380219_Covered_AgilentV5_UTR_Exomecoverage.bed |tr " " "\t" \
                > S04380219_Covered_AgilentV5_UTR_Exomecoverage_qualimap.bed

grep -w "chr1" S04380219_Covered_AgilentV5_UTR_Exomecoverage_qualimap.bed > S04380219_Covered_AgilentV5_UTR_Exomecoverage_qualimap_chr1.bed
