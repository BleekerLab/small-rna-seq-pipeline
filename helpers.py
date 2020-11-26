#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from collections import Counter
import os
from pathlib import Path
import gzip
from functools import reduce
import re
import sys

def concatenate_shorstacks_and_assign_unique_cluster_ids(list_of_shortstack_dfs, 
                                                         outfile = "concatenated_shortstacks.tsv"):
    """
    Given a list of Shortstack individual dataframes, returns a single dataframe with uniquely numbered clusters. 
    The final resulting dataframe is written to a tab-separated file
    """
    dfs = [pd.read_csv(f, sep = "\t") for f in list_of_shortstack_dfs]
    df = pd.concat(dfs) # row-bind the individual dataframes
    
    # We want every cluster to have a unique id number in the XXXX format where the number of X is the min required numeral number
    total_nb_of_clusters = df.shape[0]                         # this returns 30610 for instance
    required_number_of_digits = len(str(total_nb_of_clusters)) # this returns 5 if total_nb_of_clusters = 30610 for instance

    cluster_unique_ids = ["cluster_" + str(i+1).zfill(required_number_of_digits) for i in range(0,df.shape[0],1)]  
    cluster_unique_ids = pd.Series(cluster_unique_ids)    
    
    # Add cluster_unique_id at the beginning of the dataframe
    df.insert(0, "cluster_unique_id", cluster_unique_ids)
    
    df.to_csv(outfile, sep="\t", index = False, header = True, na_rep = "NaN")


def check_samples_tsv_file(sample_tsv_file = "config/samples.tsv"):
    """
    A function to check the validity of the input sample file provided. 
    
    Checks implemented:
      1) Checks whether the column are properly named. If not, will rename columns. 
          column1 => sample
          column2 => fastq
          column3 => genome
      2) Checks whether a dot (.) is present in the 'sample' column of the provided sample file. 
         If that is the case, stop the pipeline execution and returns an explicit error message. 
         This is to provide compatibility with multiQC.
    """
    
    # check naming of columns
    df = pd.read_csv(sample_tsv_file, sep="\t")
    colnames = df.columns.to_list()
    assert colnames[0] == "sample", "Your first column in your samples.tsv file should be named 'sample' "
    assert colnames[1] == "fastq", "Your first column in your samples.tsv file should be named 'fastq' "
    assert colnames[2] == "genome", "Your first column in your samples.tsv file should be named 'genome' "

    # check if sample names have a dot inside their name
    pattern_to_find = re.compile("\.")
    df = df.set_index("sample")
    for sample_name in list(df.index):
        if bool(pattern_to_find.search(sample_name)) == True:
            sys.exit("Please replace '.' (dots) in your sample names with another character '_' (underscore")
        else:
            return(df)
        




def create_microrna_dataframe_from_shortstack_results(list_of_sample_names,list_of_shortstack_result_files,how="outer",outfile="micrornas.merged.tsv"):
    """
    1. Takes a list of shortstack result files and create a list of Pandas dataframes
    2. For each dataframe, keep only the microRNAs (MIRNAs == Y)
    3. Select the "MajorRNA" and "MajorRNAReads" columns.
    4. Rename the MajorRNAReads column using the sample name
    5. Perform a recursive outer merge and output a single Pandas dataframe

    The option "how" can take either "outer" or "inner" as values 

    Example of output file:
    microrna    sample1    sample2
    UUUTC       10         20
    UUTTC       30         50
    ...

    """
    assert isinstance(list_of_sample_names,list), "You must provide a list containing the sample names"
    assert isinstance(list_of_shortstack_result_files,list), "You must provide a list containing paths to ShortStack results files"

    # reads each ShortStack result file
    dfs = [pd.read_csv(df,sep="\t") for df in list_of_shortstack_result_files]
    # filter the dataframes (steps 2 and 3)
    dfs_filtered = [df.query("MIRNA == 'Y'") for df in dfs]
    dfs_filtered = [df.loc[:,["MajorRNA","MajorRNAReads"]] for df in dfs_filtered]
    # rename the MajorRNAReads column (step 4)
    dfs_filtered = [df.rename(columns={"MajorRNAReads":sample}) for df,sample in zip(dfs_filtered,list_of_sample_names)]
    # perform the recursive merge
    if how == "outer":
        df_merged = reduce(lambda left,right: pd.merge(
            left.drop_duplicates("MajorRNA"),
            right.drop_duplicates("MajorRNA"),
            on=["MajorRNA"],how="outer"),
        dfs_filtered
        )
    elif how == "inner":
        df_merged = reduce(lambda left,right: pd.merge(
            left.drop_duplicates("MajorRNA"),
            right.drop_duplicates("MajorRNA"),
            on=["MajorRNA"],how="inner"),
        dfs_filtered
        )
    else:
        print("The argument 'how' only accepts 'outer' or 'inner' to merge the Shortstack results")
    
    return df_merged



def collect_clusterfiles_path(path_to_mirna_folder):
    """
    Takes the MIRNAs/ folder and returns a list of cluster files from that folder
    """
    list_of_cluster_files = [os.path.join(path_to_mirna_folder,f) for f in os.listdir(path_to_mirna_folder)]
    return list_of_cluster_files

# defines a function to read the hairpin sequence from one MIRNA cluster file
def extract_hairpin_name_and_sequence(file,sampleName):
    """
    Reads one MIRNA cluster file
    It returns the corresponding cluster name and hairpin sequence in a Python dictionary
    Dictionary keys: cluster names
    Dictionary values: cluster sequences
    """
    with open(file,"r") as filin:
        lines = filin.readlines()
        clusterName = lines[0].split(" ")[0].strip()
        hairpinSequence = lines[2].strip()
        d = {clusterName:hairpinSequence}
        return d

# takes the dictionary of sequences generated by extract_hairpin_name_and_sequence
# converts it to a fasta file
def converts_list_of_sequence_dictionary_to_fasta(list_of_sequence_dictionaries,outfastafile):
    """
    Takes the dictionary of sequences generated by extract_hairpin_name_and_sequence
    and converts it to a fasta file
    """
    with open(outfastafile,"w") as fileout:
        for sequence_dictionary in list_of_sequence_dictionaries:
            for clusterName, hairpinSequence in sequence_dictionary.items():
                fileout.write(">" + clusterName + "\n" + hairpinSequence + "\n")

# add blast header to blast result files
def add_blast_header_to_file(blast_file_without_header,blast_file_with_header):
    """takes a blast result file (outformat 6 not customized) and add a header
    """
    if os.path.getsize(blast_file_without_header) == 0:
        with open(blast_file_with_header, "w") as fileout:
            fileout.write("No blast hit. Check the miRBase databases you have used (correct species?).\n")
            fileout.write("If you run the pipeline on subset test files, you can ignore this message.") 
    
    else: 
        blast_header = ["qseqid",
                "subject_id",
                "pct_identity",
                "aln_length",
                "n_of_mismatches",
                "gap_openings",
                "q_start",
                "q_end",
                "s_start",
                "s_end",
                "e_value",
                "bit_score"]
        df = pd.read_csv(blast_file_without_header,sep="\t",header=None)
        df.columns = blast_header
        df.to_csv(blast_file_with_header,sep="\t",header=True,index=False)



def add_sample_name_to_shortstack_results(
    path_to_shortstack_results,
    sample_name):
    """
    Takes a "Results.txt" dataframe as produced by Shortstack.
    Add the sample name as a new column (sample name is repeated N times (number of rows))
    Returns a pandas dataframe as specified by outfile
    """
    df = pd.read_csv(path_to_shortstack_results,sep="\t")
    df["sample"] = sample_name
    return df


def make_dataframe_from_hairpin_fasta_file(hairpin_fasta_file):
    """
    This function will one fasta file containing all hairpins from one sample. 
    It will make a Pandas dataframe with 2 columns:
      * column 1: the hairpin sequence.
      * column 2: the sample name from which it has been extracted from.
    """
    with open(hairpin_fasta_file,"r") as filin:
        hairpins = [hairpin for hairpin in SeqIO.parse(filin,"fasta")]
    
    hairpin_identifiers = [str(hairpin.id) for hairpin in hairpins]    
    hairpin_sequences = [str(hairpin.seq) for hairpin in hairpins]
    
    df = pd.DataFrame(list(zip(hairpin_sequences,hairpin_identifiers)), columns=["hairpin","Name"])

    return df

def add_hairpin_sequence_to_shortstack_results(
    shortstack_results_dataframe,
    hairpin_fasta_file):
    """
    Add the hairpin sequence (taken from the sample hairpin fasta file) to the Shorstack dataframe.
    It uses the "Name" column as the common key for the left merge.
    """
    hairpin_df = make_dataframe_from_hairpin_fasta_file(hairpin_fasta_file)

    df = pd.merge(shortstack_results_dataframe,
                 hairpin_df,
                how="left",
                on="Name") # cluster name

    return df


def add_sample_name_and_hairpin_seq_to_shortstack(path_to_shortstack_results, sample_name, hairpin_fasta_file, outfile):
    """
    Takes a Shortstack Results.txt file, a sample name and a hairpin fasta file. 
    It returns a Shortstack result file with two additional columns (sample, hairpin) that contains the 
    sample name and hairpin sequences. 

    """
    # add sample name
    df_with_name = pd.read_csv(path_to_shortstack_results,sep="\t")
    df_with_name["sample"] = sample_name
    
    # add hairpin sequence
    df_with_name_and_hairpin = add_hairpin_sequence_to_shortstack_results(df_with_name, hairpin_fasta_file)

    # write to file
    df_with_name_and_hairpin.to_csv(outfile, sep="\t", index=False, header=True, na_rep = "NaN")


def concatenate_shortstacks_and_assign_unique_cluster_ids(list_of_shortstack_files,
                                                         outfile = "shortstack_concatenated.tsv"):
    list_of_shortstack_dfs = [pd.read_csv(f,sep="\t") for f in list_of_shortstack_files]
    df = pd.concat(list_of_shortstack_dfs)
    df["cluster_unique_id"] = ["cluster_" + str(i+1).zfill(10) for i in range(0,df.shape[0],1)]  
    df.to_csv(outfile, sep="\t", index=False, header=True)

def extract_hairpins_from_concatenated_shortstack_file(concatenated_shortstack_file, outfile):
    """
    Extract the hairpin sequences from the concatenated shortstack dataframe.
    The hairpin identifier is a unique cluster identifier.
    """
    concatenated_shortstack_df = pd.read_csv(concatenated_shortstack_file, sep="\t")
    df_with_only_mirnas = concatenated_shortstack_df.query(" MIRNA == 'Y' ") # only true MIRNAs have a hairpin sequence

    cluster_ids = df_with_only_mirnas["cluster_unique_id"]
    hairpin_seqs = df_with_only_mirnas["hairpin"]

    with open(outfile, "w") as fileout:
        for cluster_id, hairpin_seq in zip(cluster_ids, hairpin_seqs):
                fileout.write(">" + cluster_id + "\n" + hairpin_seq + "\n")


def extract_mature_micrornas_from_concatenated_shortstack_file(concatenated_shortstack_file, outfile):
    """
    Extract the mature miRNA sequences from the concatenated shortstack dataframe.
    The mature miRNA identifier is a unique cluster identifier.
    """
    concatenated_shortstack_df = pd.read_csv(concatenated_shortstack_file, sep="\t")
    df_with_only_mirnas = concatenated_shortstack_df.query(" MIRNA == 'Y' ") # only true MIRNAs are kept

    cluster_ids = df_with_only_mirnas["cluster_unique_id"]
    major_rna_seqs = df_with_only_mirnas["MajorRNA"]

    with open(outfile, "w") as fileout:
        for cluster_id, major_rna_seq in zip(cluster_ids, major_rna_seqs):
                fileout.write(">" + cluster_id + "\n" + major_rna_seq + "\n")


def extract_mature_mirna_fasta_file_from_shortstack_file(path_to_shortstack_file, out_fasta_file):
    """
    Takes a Shortstack Results.txt file and
    output a fasta file with all predicted mature miRNAs
    """
    # converts all mature miRNAs into a dictionary 
    df = pd.read_csv(path_to_shortstack_file, sep="\t", index_col=1)
    df_mirnas = df[df["MIRNA"] == "Y"]
    mirnas_dict = df_mirnas["MajorRNA"].to_dict()

    # convert to fasta format
    with open(out_fasta_file, "w") as fileout:
        for name, sequence in mirnas_dict.items():
            fileout.write(">" + name + "\n" + sequence + "\n")


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
    merged_df_sorted.to_csv(path_or_buf = outfile, sep = "\t", index = True, na_rep = "NA")

if __name__ == "__main__":
    print("main script")


