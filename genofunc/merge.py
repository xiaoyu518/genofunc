"""
Name: merge.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Merges two or more fasta files avoiding duplicates based on matches to metadata.

Takes the first appearance according to the sequence of files within the input.
At least two fasta files must be within the input command and only those sequences matching metadata
will be processed into output fasta file.

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import sys
from genofunc.utils import *

"""
Merges two or more fasta files avoiding duplicates based on matches to metadata

:param in_fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. Only
fasta files are taken as input. (Required)
:param in_metadata: list of matching metadata file with same naming convention as fasta file (index-column). Those
that does not match or have duplicates will be flagged within the log file for post-processing. Metadata file must
be in csv format (Required)
:param index_column: The column ID with matching sequence IDs with fasta file (Required)
:param out_dir: Output metadata file with merged columns from multiple inputs (Default: ./).
:return:
"""

def merge(in_fasta,in_metadata,index_field,out_dir,log_file):
    for files in in_fasta:
        if not file_check(files):
            print("Fasta file " + files + " does not exist. Please enter new one. Program Exiting.")
            sys.exit()

    for metadata in in_metadata:
        if not file_check(metadata):
            print("Metadata file " + metadata + " does not exist. Please enter new one. Program Exiting.")
            sys.exit()

    out_fasta = open(out_dir + "merged.fasta","w")
    out_metadata = open(out_dir + "merged_metadata.csv","w",newline='')
    log_handle = get_log_handle(log_file)

    sequence_dictionary = {}
    metadata_dictionary = {}

    for metadata_file in in_metadata:
        with open(metadata_file,"r") as f:
            seperator = ','
            if metadata_file.endswith('tsv'):
                seperator = '\t'
            reader = csv.DictReader(f,delimiter=seperator)
            reader.fieldnames = [name.lower() for name in reader.fieldnames]
            column_names = reader.fieldnames
            metadata = [r for r in reader]
        for sequence in metadata:
            if index_field not in column_names:
                print("Index column not in metadata. Please re-enter a new one. Program exiting.")
                sys.exit()
            else:
                taxon_name = sequence[index_field]
            if taxon_name not in metadata_dictionary.keys():
                metadata_dictionary[taxon_name] = sequence
            else:
                log_handle.write("Sequence " + taxon_name + " had a duplicate in metadata and the new metadata is removed\n")

    sequence_list = list(metadata_dictionary.keys())
    out_list = list(metadata_dictionary.values())

    f = csv.DictWriter(out_metadata,fieldnames=column_names)
    f.writeheader()
    f.writerows(out_list)
    close_handle(out_metadata)

    for fasta_file in in_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in sequence_list and record.id not in sequence_dictionary.keys():
                sequence_dictionary[record.id] = record.seq
            elif record.id not in sequence_list:
                log_handle.write(record.id + " sequence is not in metadata file or the name is wrong (in file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                log_handle.write(record.id + " is a duplicate (in file " + fasta_file + ")\n")

    if len(sequence_dictionary.keys()) == 0:
        print("There is no matching sequences to metadata. Program exiting")
        sys.exit()

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, out_fasta, "fasta-2line")

    close_handle(out_fasta)
    close_handle(log_handle)
