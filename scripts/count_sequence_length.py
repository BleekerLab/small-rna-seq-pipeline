#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from collections import Counter
import os
from pathlib import Path
import gzip
from functools import reduce


# Helper functions

def collect_list_of_fastq_files_from_directory(path_to_directory_with_fastq_files):
    """
    Given a directory, returns a list of fastq files with their path.
    Fastq files have to end with '.fastq', '.fastq.gz', '.fq', '.fq.gz'
    """
    # validates path to directory
    if Path(path_to_directory_with_fastq_files).is_dir():
        pass
    else:
        print("please provide a valid directory")
    
    # authorised file extensions
    extensions_ok = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    
    fastq_files = [path_to_directory_with_fastq_files + f for f in os.listdir(path_to_directory_with_fastq_files) if f.endswith(extensions_ok)]
    
    assert len(fastq_files) > 0, "Empty list of fastq files. Does your directory contains fastq files?"
    
    return(fastq_files)


def get_count_df_of_seq_length_from_fastq(fastq_file):
    """
    From a fastq file, generates a Pandas dataframe with the seq lengths + number of occurences.
    One sample = one column
    N seq lengths = N rows
    """
    if fastq_file.endswith(".fastq") or fastq_file.endswith(".fq"):
        with open(fastq_file, "r") as filin:
            seqs = [str(rec.seq) for rec in SeqIO.parse(filin,"fastq")]
            sequence_lengths = [len(seq) for seq in seqs]
            sequence_lengths_dict = Counter(sequence_lengths)        
    if fastq_file.endswith(".fastq.gz") or fastq_file.endswith(".fq.gz"):
        with gzip.open(fastq_file, "rt") as filin:
            seqs = [str(rec.seq) for rec in SeqIO.parse(filin,"fastq")]
            sequence_lengths = [len(seq) for seq in seqs]
            sequence_lengths_dict = Counter(sequence_lengths)        
        
    # sample name for column in Pandas dataframe
    if fastq_file.endswith(".fastq"):
        col_name = os.path.basename(fastq_file).split(".fastq")[0]
    if fastq_file.endswith(".fq"):
        col_name = os.path.basename(fastq_file).split(".fq")[0]
    if fastq_file.endswith(".fastq.gz"):
        col_name = os.path.basename(fastq_file).split(".fastq.gz")[0]
    if fastq_file.endswith(".fq.gz"):
        col_name = os.path.basename(fastq_file).split(".fq.gz")[0]
            
    # Pandas dataframe (one column, N rows per seq length)
    df = pd.DataFrame.from_dict(data = sequence_lengths_dict, orient="index", columns=[col_name])
    df.index.name = 'length'
            
    return(df)
        

### Main function
def create_df_of_seq_length_distributions(path_to_fastq_files, outfile = "./sequence_length_distribution.tsv"):
    """
    Given a path to a directory with fastq files, creates a Pandas dataframe of sequence length distribution
    Fastq files have to end with '.fastq', '.fastq.gz', '.fq', '.fq.gz'
    """
    list_of_fastq_files = collect_list_of_fastq_files_from_directory(path_to_fastq_files)
    print("Working on files:", list_of_fastq_files)
    print("Results are saved in:", outfile)
    
    # Creates one df per count sequences
    dfs_of_seq_len_distri = [get_count_df_of_seq_length_from_fastq(f) for f in list_of_fastq_files]
    
    # Merge into one df
    merged_df = reduce(lambda x, y: pd.merge(x, y, how = "outer", left_on = 'length', right_on = 'length'), 
                       dfs_of_seq_len_distri)
    
    # Sort by increasing length
    merged_df_sorted = merged_df.sort_index()
    
    # write to tabulated separated values file
    merged_df_sorted.to_csv(path_or_buf = outfile, sep = "\t", index = True)

if __name__ == "__main__":
   create_df_of_seq_length_distributions(path_to_fastq_files, outfile = "./sequence_length_distribution.tsv")

