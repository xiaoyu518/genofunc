#!/usr/bin/env python3

"""
Name: feature_extractor.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Extract gene(s) sequences from annotated json file based on user input gene regions.

Options:
    :param in_annotation: Annotated json file containing all sequences (Required)
    :param gene_region: Gene regions to be extracted (Required)
    :param strip_gap: Strip gap bases within gene regions (Default: False)
    :param filter_span: Minimum gene sequence length to be filtered (Default: 0)
    :param output_prefix: Output prefix for output sequences (Default: extracted_)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import datetime as dt
from genofunc.utils import *

def feature_extractor(in_annotation,gene_region,strip_gap,filter_span,output_prefix,log_file):
    time_start = dt.datetime.now()
    output_dic = {}
    
    if file_check(in_annotation):
        with open(in_annotation,"r") as f:
            input_data = json.load(f)
    else:
        print("Input file does not exist. Please check the file path.")
        sys.exit()

    log_handle = get_log_handle(log_file)

    for gene in gene_region:
        output_dic[gene] = {}

    for sequences in input_data.keys():
        strain_id = sequences
        for gene in gene_region:
            coordinates = input_data[sequences][gene].split("|")
            temp_seq = ""
            span = 0.0
            if len(coordinates) == 2:
                begin = int(coordinates[0])
                end = int(coordinates[1])
                temp_seq = input_data[sequences]["sequence"][begin-1:end-1]
                length = end-begin
            if len(coordinates) == 4:
                begin = [int(coordinates[0]),int(coordinates[2])]
                end = [int(coordinates[1]),int(coordinates[3])]
                temp_seq = input_data[sequences]["sequence"][begin[0]-1:end[0]-1] + input_data[sequences]["sequence"][begin[1]-1:end[1]-1]
                length = (end[0] + end[1]) - (begin[0] + begin[1])
            if strip_gap:
                temp_seq = temp_seq.replace("-","")
            span = len(temp_seq)/length
            if span < filter_span:
                log_handle.write(strain_id + " " + gene + " gene region sequence length is shorter than the minimum required span length "
                + str(filter_span) + " and therefore filtered out.\n")
                continue
            else:
                output_dic[gene][strain_id] = temp_seq

    for k in output_dic.keys():
        outfile = open(output_prefix + k + "_extracted.fasta","w")
        for id, seq in output_dic[k].items():
            new_record = SeqRecord(Seq(seq),id,description="")
            SeqIO.write(new_record, outfile, "fasta-2line")
        close_handle(outfile)

    close_handle(log_handle)

    time_ran = dt.datetime.now() - time_start
    print("Time Lapse:", time_ran.total_seconds() / 60, "Minutes")
