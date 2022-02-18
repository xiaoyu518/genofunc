#!/usr/bin/env python3

"""
Name: rename_fasta.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Renaming fasta sequence names based on character splits.

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from genofunc.utils import *

"""
Renaming fasta sequence names based on character splits.

:param in_fasta: Input fasta file (Required)
:param piping_character: Input character the fasta sequence name is split based on (Required)
:param out_fasta: Output fasta file (Default: cleaned_sequences.fasta)
:return:
"""

def rename_fasta(in_fasta,piping_character,out_fasta,log_file):
    output_fasta = open(out_fasta,"w")
    log_handle = get_log_handle(log_file)

    if file_check(in_fasta):
        for record in SeqIO.parse(in_fasta, "fasta"):
            if record.id.find(piping_character) > -1:
                new_id = record.id[:record.id.find(piping_character)]
                new_record = SeqRecord(record.seq,new_id,description="")
                SeqIO.write(new_record, output_fasta, "fasta-2line")
            else:
                log_handle.write(record.id + " does not contain the split character. Whole ID used.\n")
                SeqIO.write(record, output_fasta, "fasta-2line")

    else:
        print("Fasta file does not exist. Please enter a new file path. Program Exiting")
        sys.exit()

    close_handle(output_fasta)
    close_handle(log_handle)
