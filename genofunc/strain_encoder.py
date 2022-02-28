#!/usr/bin/env python3

"""
Name: strain_encoder.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Encoded strain id into non-defining ids. 

Options:
    :param in_fasta: Input folder containing alignments in fasta format (Required)
    :param in_metadata: Input metadata for cohort information (Required)
    :param encoding_column: Column for the base for encoding information (Required)
    :param out_dir: Output folder including encoded files (Default: ./)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import glob
from genofunc.utils import *

def strain_encoder(in_fasta,in_metadata,encoding_column,out_dir,log_file):

    fasta_files = glob.glob(in_fasta + "*.fasta")

    if fasta_files == []:
        print("No fasta files in input directory. Please re-enter new one. Program Exiting.")
        sys.exit()

    if not file_check(in_metadata):
        print("Metadata file does not exist. Please enter a new file path. Program Exiting")
        sys.exit()


    output_metadata = open(out_dir+"encoded_metadata.csv","w")
    log_handle = get_log_handle(log_file)
    column_dic = {}
    encoding_dic = {}

    with open(in_metadata,"r") as f:
        seperator = ','
        if in_metadata.endswith('tsv'):
            seperator = '\t'
        for line in f:
            rows = line.rstrip().split(seperator)
            if rows[encoding_column] not in column_dic.keys():
                column_dic[rows[encoding_column]] = [rows[0]]
            else:
                column_dic[rows[encoding_column]].append(rows[0])

    for k,v in column_dic.items():
        for i in range(len(v)):
            temp_encode = "-".join([k,str(i+1)])
            encoding_dic[v[i]] = temp_encode
            log_handle.write("\t".join([v[i],temp_encode])+"\n")

    for files in fasta_files:
        file_name = files[files.rfind("/")+1:]
        output_sequence = open(out_dir+file_name[:-6]+"_encoded.fasta","w")
        for record in SeqIO.parse(files, "fasta"):
            new_id = encoding_dic[record.id]
            new_record = SeqRecord(record.seq,new_id,description="")
            SeqIO.write(new_record, output_sequence, "fasta-2line")
        output_sequence.close() 

    with open(in_metadata,"r") as f:
        seperator = ','
        if in_metadata.endswith('tsv'):
            seperator = '\t'
        for line in f:
            rows = line.rstrip().split(seperator)
            if rows[0] == "strain":
                output_metadata.write(line)
                continue
            rows[0] = encoding_dic[rows[0]]
            output_metadata.write(",".join(rows)+"\n")

    close_handle(output_metadata)
    close_handle(log_handle)
