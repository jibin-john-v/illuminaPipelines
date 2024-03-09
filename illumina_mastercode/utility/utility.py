import sys,os,glob
from pathlib import Path
import pandas as pd
import concurrent.futures

##check a file is gziped or not
def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

#uncompress_if_gzipped
def uncompress_if_gzipped(file_path):
    def is_gzipped(file_path):
        with open(file_path, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'

    if is_gzipped(file_path):
        output_file_path = file_path.rstrip('.gz')
        os.system(f"gunzip {file_path}")
        print(f'File uncompressed to {output_file_path}')
        return output_file_path
    else:
        print(f'{file_path} is not gzipped.')
        return file_path


##Extract mirna fasta file of specific species
def extract_mirna_fasta(inputfasta,species,output_fasta):
    os.system(f'perl /mnt/disks/sdc/Software/mirdeep2_patch/extract_miRNAs.pl \
                    {inputfasta} {species} > {output_fasta}')


##create input csv file with sample name and fastq location
def fastq_files_location_df(fastq_path):
    # Supported fastq extensions
    supported_extensions_paired = {'_R1.fastq.gz': '_R2.fastq.gz', '_1.fastq.gz': '_2.fastq.gz',
                                   '_R1.fastq': '_R2.fastq', '_1.fastq': '_2.fastq',
                                   '_R1.fq.gz': '_R2.fq.gz', '_1.fq.gz': '_2.fq.gz',
                                   '_R1.fq': '_R2.fq', '_1.fq': '_2.fq'}
    supported_extensions_single = ['.fastq.gz','_fastq.gz']
    
    paired_end_fasq_df = pd.DataFrame()
    single_end_fasq_df = pd.DataFrame()
    
    for extension_1, extension_2 in supported_extensions_paired.items():
        fastq_files = glob.glob(f"{fastq_path}*{extension_1}")
        if fastq_files:
            fastq_df = pd.DataFrame(fastq_files, columns=['File'])
            fastq_df["SampleName"] = fastq_df["File"].str.replace(f"{extension_1}", "").str.split("/", expand=True).iloc[:, -1]
            fastq_df["R1"] = fastq_df["File"].str.replace(f"{extension_1}", "", regex=False) + extension_1
            fastq_df["R2"] = fastq_df["File"].str.replace(f"{extension_1}", "", regex=False) + extension_2
            fastq_df=fastq_df.drop("File",axis=1)
            paired_end_fasq_df = pd.concat([paired_end_fasq_df, fastq_df])
    
    for extension in supported_extensions_single:
        fastq_files = glob.glob(f"{fastq_path}*{extension}")
        if fastq_files:
            fastq_df = pd.DataFrame(fastq_files, columns=['File'])
            fastq_df["SampleName"] = fastq_df["File"].str.replace(f"{extension}", "", regex=False).str.split("/", expand=True).iloc[:, -1]
            fastq_df["R1"] = fastq_df["File"].str.replace(f"{extension}", "", regex=False) + extension
            fastq_df=fastq_df.drop("File",axis=1)
            single_end_fasq_df = pd.concat([single_end_fasq_df, fastq_df])
    
    if not paired_end_fasq_df.empty:
        for fastq_file in paired_end_fasq_df['R1'].unique():
            if not Path(fastq_file).is_file():
                print(f"{fastq_file} does not exist")
                sys.exit(1)
        for fastq_file in paired_end_fasq_df['R2'].unique():
            if not Path(fastq_file).is_file():
                print(f"{fastq_file} does not exist")
                sys.exit(1)
        paired_end_fasq_df.to_csv("paired_end_fastq_files_locations.csv", index=None)
    
    if not single_end_fasq_df.empty:
        for fastq_file in single_end_fasq_df['R1'].unique():
            if not Path(fastq_file).is_file():
                print(f"{fastq_file} does not exist")
                sys.exit(1)
        single_end_fasq_df.to_csv("single_end_fastq_files_locations.csv", index=None)
    
    if paired_end_fasq_df.empty and single_end_fasq_df.empty:
        sys.exit(1)
    return paired_end_fasq_df, single_end_fasq_df



#file gziping much faster https://www.zlib.net/pigz/

##effective genome size calculation; https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/# ;https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html