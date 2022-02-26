#!/usr/bin/env python3

"""
Name: concatenate_fasta.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Concatenate all fasta file with the same gene regions into one sequence fasta file. 

Options:
    :param in_prefix: Input prefix within a directory for specific fasta files to be read (Required)
    :param gene_region: Gene regions to concatenate by which fasta file names are included in (Required)
    :param out_dir: Output directory after concatenating all fasta files with the same gene name (Default: ./)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from Bio import SeqIO
import glob
from genofunc.utils import *

def concatenate_fasta(in_prefix,gene_region,out_dir,log_file):
    log_handle = get_log_handle(log_file)

    for genes in gene_region:
        fasta_files = glob.glob(in_prefix + "*" + genes + '*.fasta')
        if fasta_files == []:
            log_handle.write("Prefix does not contain " + genes + " files to concatenate.\n")
        output_file = out_dir + genes + ".fasta"
        out_fasta = open(output_file,"w")
        for files in fasta_files:
            for record in SeqIO.parse(files, "fasta"):
                SeqIO.write(record, out_fasta, "fasta-2line")
            log_handle.write(files + " has concatenated into " + output_file + "\n")    
        close_handle(out_fasta)

    close_handle(log_handle)
