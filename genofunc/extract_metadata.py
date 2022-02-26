#!/usr/bin/env python3

"""
Name: extract_metadata.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Extract relevent metadata based on index fields from metadata file using fasta files. If field is empty, it will not be extracted and flagged in log file.

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

import sys
import csv
import pandas as pd
from Bio import SeqIO
from genofunc.utils import *

"""
Extract relevent metadata based on index fields from metadata file using fasta files.

:param in_fasta: Input fasta file (Required)
:param in_metadata: Input metadata file (Required)
:param column_name: column name for extraction (Required)
:param column_rename: column name for re-naming 
:param id_column: column name containing sequence ID (Required)
:param out_fasta: Output fasta file (Default: filtered.fasta)
:param out_metadata: Output metadata file (Default: extracted_metadata.csv)
:return:
"""

def extract_metadata(in_fasta,in_metadata,column_name,column_rename,id_column,out_fasta,out_metadata,log_file):

    output_fasta = open(out_fasta,"w")
    log_handle = get_log_handle(log_file)
    output_metadata = open(out_metadata,"w")
    seq_list = []

    if file_check(in_fasta):
        for record in SeqIO.parse(in_fasta, "fasta"):
            seq_list.append(record.id)
    else:
        print("Fasta file does not exist. Please enter a new file path. Program Exiting")
        sys.exit()

    if file_check(in_metadata):
        seperator = ','
        if in_metadata.endswith('tsv'):
            seperator = '\t'
        metadata = pd.read_csv(in_metadata,sep=seperator)
        metadata.columns = map(str.lower, metadata.columns)
    else:
        print("Metadata file does not exist. Please enter a new file path. Program Exiting")
        sys.exit()

    for column in column_name:
        if column not in metadata.columns:
            print(column + " is not a column name in metadata. Please re-check column names. Program Exiting.")
            sys.exit()

    if id_column not in column_name:
        column_name.append(id_column)

    filtered_df = metadata.loc[:, metadata.columns.isin(column_name)]
    if column_rename != []:
        filtered_df.columns = column_rename

    extracted_df = filtered_df[filtered_df[id_column].isin(seq_list)]
    out_df = extracted_df.dropna(axis=0, how='any', thresh=None, subset=None)
    out_df = out_df.drop_duplicates(subset=id_column, keep='first')

    out_list = out_df[id_column].tolist()

    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id in out_list:
            SeqIO.write(record, output_fasta, "fasta-2line")
        else:
            log_handle.write(record.id + " has been filtered due to unmatched metadata")
    
    out_df.to_csv(output_metadata, index = False)

    close_handle(output_fasta)
    close_handle(output_metadata)
    close_handle(log_handle)
