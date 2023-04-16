#!/usr/bin/env python3

"""
Name: feature_extractor.py
Author: Xiaoyu Yu
Date: 10 April 2023
Description: Extract gene(s) sequences from annotated json file based on user input gene regions.

Options:
    :param in_annotation: Annotated json file containing all sequences (Required)
    :param gene_region: Gene regions to be extracted (Required)
    :param strip_gap: Strip gap bases within gene regions (Default: False)
    :param filter_coverage: Minimum base ratio of non N in gene sequence length to be filtered (Default: 0)
    :param filter_span: Minimum sequence length to be filtered (Default: 0)
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

def feature_extractor(in_annotation,gene_region,strip_gap,filter_coverage,filter_span,output_prefix,log_file):
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
            coordinates = input_data[strain_id][gene].split("|")
            temp_seq = ""
            coverage = 0.0
            length = 0
            if len(coordinates) > 1:
                for i in range(0,len(coordinates),2):
                    begin = int(coordinates[i])
                    end = int(coordinates[i+1])
                    if begin < 0:
                        temp_seq += str(Seq(input_data[strain_id]["sequence"][abs(begin) - 1:abs(end)]).reverse_complement())
                    else:
                        temp_seq += input_data[strain_id]["sequence"][begin - 1:end]
                    length += abs(end) - abs(begin) + 1           
            else:
                print("Sequence " + strain_id + " with gene " + gene + " does not contain valid gene coordinate for extraction. Program Exiting.")
                sys.exit()
                
            temp_seq = temp_seq.replace("N","-")
            if strip_gap:
                temp_seq = temp_seq.replace("-","")            
            coverage = round(float(len(temp_seq))/length,2)
            if len(temp_seq) < filter_span:
                log_handle.write(strain_id + " " + gene + " gene region sequence length is shorter than the minimum required span length "
                + str(filter_span) + " and therefore filtered out.\n")
                continue
            else:
                output_dic[gene][strain_id] = temp_seq                
            if coverage < filter_coverage:
                log_handle.write(strain_id + " " + gene + " gene region sequence base ration less than the minimum required coverage "
                + str(filter_coverage) + " and therefore filtered out.\n")
                continue
            else:
                output_dic[gene][strain_id] = temp_seq

    for k in output_dic.keys():
        outfile = open(output_prefix + k + ".fasta","w")
        for id, seq in output_dic[k].items():
            new_record = SeqRecord(Seq(seq),id,description="")
            SeqIO.write(new_record, outfile, "fasta-2line")
        close_handle(outfile)

    time_ran = dt.datetime.now() - time_start
    print("Time Lapse:", time_ran.total_seconds() / 60, "Minutes")

    close_handle(log_handle)
