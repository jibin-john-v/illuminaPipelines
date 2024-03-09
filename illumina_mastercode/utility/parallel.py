import sys,os,glob
from pathlib import Path
import pandas as pd
import concurrent.futures

sys.path.append('../fastq_qc')
from fastq_qc.fastq_qc import falco_run,fastp_run

#Raw data fastqc
def run_fastqc_parallel(method,fastq_df, n_parallel=5,n_threads=3,out_name="output",qc_cutoff=20):
    with concurrent.futures.ThreadPoolExecutor(n_parallel) as executor:
        futures = []
        
        for name in fastq_df["SampleName"].unique():
            if 'R2' in fastq_df.columns:
                r1 = fastq_df[fastq_df["SampleName"] == name]['R1'].values[0]
                r2 = fastq_df[fastq_df["SampleName"] == name]['R2'].values[0]
                reads=[r1, r2]
            else:
                r1 = fastq_df[fastq_df["SampleName"] == name]['R1'].values[0]
                reads=[r1]
            if method=="falco":
                future = executor.submit(falco_run, name, reads,out_name)
                futures.append(future)
            elif method=="fastqc":
                future = executor.submit(fastqc_run, name, reads,n_threads)
                futures.append(future)
            elif method=="fastp":
                future = executor.submit(fastp_run, name, reads,n_threads,qc_cutoff)
                futures.append(future)
        # Wait for each job to complete and start the next one
        for future in concurrent.futures.as_completed(futures):
            # This will block until a job is completed, allowing the next one to start
            pass
