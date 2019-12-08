# helper functions
import os
import pandas as pd
from functools import reduce
from Bio import SeqIO

# 

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
    df_with_name_and_hairpin.to_csv(outfile, index=False, header=True, na_rep = "NaN")