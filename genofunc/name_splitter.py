#!/usr/bin/env python3

"""
Name: name_splitter.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Split the sequence name into metadata based on piping character.

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from genofunc.utils import *

"""
Split the sequence name into metadata based on piping character.

:param in_fasta: Input fasta file (Required)
:param piping_character: Input character the fasta sequence name is split based on (Required)
:param header: Header for the output metadata table (Default: "")
:param out_metadata: Output metadata file (Default: metadata.csv)
:return:
"""

def name_splitter(in_fasta,piping_character,header,out_metadata,log_file):
    log_handle = get_log_handle(log_file)

    if header == "":
        metadata = []
    else:
        metadata = [header]

    if file_check(in_fasta):
        for record in SeqIO.parse(in_fasta, "fasta"):
            if record.id.find(piping_character) == -1:
                log_handle.write(record.id + " does not contain the split character. No metadata created\n")
                sys.exit()
            else:
                row = record.id.split(piping_character)
                metadata.append(row) 
    else:
        print("Fasta file does not exist. Please enter a new file path. Program Exiting")
        sys.exit()

    with open(out_metadata,"w") as f:
        writer = csv.writer(f)
        writer.writerows(metadata)

    close_handle(log_handle)
