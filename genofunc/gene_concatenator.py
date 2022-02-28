#!/usr/bin/env python3

"""
Name: gene_concatenator.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Concatenate fasta file based on similar sequence names of multiple genomic regions. 

Options:
    :param in_fasta: Multiple fasta files for concatenation (Required)
    :param out_fasta: Output fasta file (Default: concatenate.fasta)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
from genofunc.utils import *

def gene_concatenator(in_fasta,out_fasta,log_file):
    output_fasta = open(out_fasta,"w")
    log_handle = get_log_handle(log_file)
    seq_dic = {}

    for files in in_fasta:
        if not file_check(files):
            print("Fasta file " + files + " does not exist. Please enter new one. Program Exiting.")
            sys.exit()

    for files in in_fasta:
        for record in SeqIO.parse(files, "fasta"):
            if record.id not in seq_dic.keys():
                seq_dic[record.id] = str(record.seq)
                log_handle.write(record.id + " has been added to fasta file. \n")
            else:
                seq_dic[record.id] = seq_dic[record.id] + str(record.seq)
                log_handle.write(record.id + " has been concatenated \n")

    for key, value in seq_dic.items():
        records = SeqRecord(Seq(value), key, description= '')
        SeqIO.write(records, output_fasta, "fasta-2line")
            
    close_handle(output_fasta)
    close_handle(log_handle)
