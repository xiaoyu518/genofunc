#!/usr/bin/env python3

"""
Name: filter_fasta.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Filter fasta file based on minimum length threshold. 

Options:
    :param in_dir: Input directory (Required)
    :param in_metadata: Input metadata for filtering alongside sequence file (Required)
    :param gene_list: Which genes fasta files are needed for filtering (Required)
    :param min_length: Minimum percentage for sequences to be filtered if under this threshold (Required)
    :param symmetric: Require all gene regions to be available for the same sequence (Default: False)
    :param out_dir: Output directory after filtering (Default: ./)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob
import sys
import os
import math
import csv
from genofunc.utils import *
from functools import reduce

def filter_fasta(in_dir,in_metadata,gene_list,min_length,symmetric,out_dir,log_file):
    log_handle = get_log_handle(log_file)
    fasta_files = glob.glob(in_dir + "*.fasta")
    gene_dic = {}
    all_list = []

    if fasta_files == []:
        print("This folder is empty. Please re-enter a new input folder for processing. Program Exiting.")
        sys.exit()

    for files in fasta_files:
        seq_length_list = []
        output_file = out_dir + files[files.rfind("/")+1:-6]+"_filtered.fasta"
        out_f = open(output_file,"w")
        for record in SeqIO.parse(files, "fasta"):
            seq_length_list.append(len(record.seq))
        average_len = sum(seq_length_list)/len(seq_length_list)
        min_threshold = math.trunc(min_length*average_len)
        flag = False
        for gene in gene_list:
            if files.find(gene) > -1:
                gene_dic[gene] = []
                file_gene = gene
                flag = True
        if flag == False:
            log_handle.write("Folder does not contain " + gene + " file.\n")                       
        for record in SeqIO.parse(files, "fasta"):
            if len(record.seq) >= min_threshold:
                if record.id.find("|") > -1:
                    record.id = record.id[:record.id.find("|")]
                gene_dic[file_gene].append(record.id)
                records = SeqRecord(record.seq, record.id, description= '')
                SeqIO.write(records, out_f, "fasta-2line")
            else:
                log_handle.write("Sequence " + record.id + " contains less than " + str(min_length) + " base positions in file " + output_file + "\n")
        close_handle(out_f)
        all_list.append(gene_dic[file_gene]) 

    filtered_files = glob.glob(out_dir + "*_filtered.fasta")

    with open(in_metadata,"r") as f:
        seperator = ','
        if in_metadata.endswith('tsv'):
            seperator = '\t'
        reader = csv.reader(f,delimiter=seperator)
        metadata = list(reader)
    
    os.remove(in_metadata)

    if symmetric is False:
        for files in filtered_files:
            os.rename(files,out_dir + files[files.rfind("/")+1:-15] + "_filter.fasta")

        metadata_dic = {}

        for gene in gene_dic.keys():
            metadata_dic[gene] = [metadata[0]]

            for rows in metadata:
                if rows[0] in gene_dic[gene]:
                    metadata_dic[gene].append(rows)
                else:
                    log_handle.write("Sequence " + rows[0] + " has been filtered out from " + gene + " metatdata due to sequence not existing in fasta file.\n")

            with open(in_metadata[:in_metadata.rfind("/")+1] + gene + "_metadata.csv","w") as f:
                writer = csv.writer(f)
                writer.writerows(metadata_dic[gene])

    else: 
        sequence_list = list(reduce(lambda i, j: i & j, (set(x) for x in all_list)))

        for files in filtered_files:
            output_file = open(out_dir + files[files.rfind("/")+1:-15] + "_filter.fasta","w")
            for record in SeqIO.parse(files, "fasta"):
                if record.id in sequence_list:
                    SeqIO.write(record, output_file, "fasta-2line")
                else:
                    log_handle.write("Sequence " + record.id + " does not have a sequence for all three genome regions and is removed from sequence file.\n")
            close_handle(output_file)
            os.remove(files)
        
        new_metadata = [metadata[0]]

        for rows in metadata:
            if rows[0] in sequence_list:
                new_metadata.append(rows)
            else:
                log_handle.write("Sequence " + rows[0] + " has been filtered out due to sequence not existing in fasta file.\n")

        with open(in_metadata,"w") as f:
            writer = csv.writer(f)
            writer.writerows(new_metadata)

    close_handle(log_handle)
