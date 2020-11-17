#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from Bio import SeqIO
import subprocess
from collections import Counter
import os


# # Helper functions

# In[10]:


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
            
            # sample name for column in Pandas dataframe
            if fastq_file.endswith(".fastq"):
                col_name = os.path.basename(fastq_file).split(".fastq")[0]
            if fastq_file.endswith(".fq"):
                col_name = os.path.basename(fastq_file).split(".fq")[0]
            
        # Pandas dataframe (one column, N rows per seq length)
    df = pd.DataFrame.from_dict(data = sequence_lengths_dict, orient="index", columns=[col_name])
            
    return(df)
        


# In[11]:


get_count_df_of_seq_length_from_fastq("test/test.fastq")


# In[ ]:




