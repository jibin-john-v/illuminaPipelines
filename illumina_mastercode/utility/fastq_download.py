

#1) nf-core/fetchngs


#2)pysradb ; https://github.com/saketkc/pysradb
os.system(f'''pysradb metadata SRP265425 --detailed > SRP265425_detailed_metadata.csv''')
os.system(f'''pysradb download -t 8 --use_ascp -p 	SRP082156''')
os.system(f'''mkdir -p seurat_objofastq && mkdir -p tmpdir''')
os.system(f'''parallel-fastq-dump \
                        --threads 4 \
                        --outdir seurat_objofastq/ \
                        --split-files \
                        --tmpdir tmpdir \
                        --gzip -s pysradb_downloads/SRP063852/SRX1254413/SRR2433794.sra''')
          


for file in pysradb_downloads/SRP082156/* ; do echo $file ; done
    
parallel-fastq-dump --sra-id pysradb_downloads/SRP082156/SRX2023278/SRR4031769.sra --threads 4 --outdir out/ --split-files --gzip


#3) fetchfastq ; https://github.com/pachterlab/ffq

os.system(f'''ffq --ftp SRP082156 | jq -r '.[] | .url' | xargs curl -O''')

#4) fasterq-dump ; https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump : SRA-Toolkit

#5)sra-explorer : https://sra-explorer.info/ 

#6) geofetch : https://github.com/pepkit/geofetch
