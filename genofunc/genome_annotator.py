#!/usr/bin/env python3

"""
Name: genome_annotator.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Annotate genomes based on closest reference sequence annotation. 

Options:
    :param raw_fasta: Raw sequences with reference in name tag in fasta format (Required)
    :param reference_sequence: Annotated reference sequences in json format (Required)
    :param annotated_json: Output list of sequences annotated in json format (Default: referenced.json)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2021 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

import json
import sys
from Bio import SeqIO
import datetime as dt
from genofunc.utils import *

def genome_annotator(raw_fasta,reference_sequence,annotated_json,log_file):
    time_start = dt.datetime.now()
    location_dic = {}
    reference_dic = {}
    features = []

    if not file_check(raw_fasta):
        print("Input file does not exist. Please check the file path.")
        sys.exit()
    if not file_check(reference_sequence):
        print("Reference file does not exist. Please check the file path.")
        sys.exit()    

    log_handle = get_log_handle(log_file)

    with open(reference_sequence,"r") as f:
        reference_data = json.load(f)

    for k,v in reference_data.items():
        for i in v:
            location_dic[i["accession"]] = {}
            location_dic[i["accession"]]["sequence"] = i["sequence"]
            for genes in i["locations"]["location"]:
                temp_list = []
                if genes["featureId"] not in features:
                    features.append(genes["featureId"])
                if isinstance(genes["fragment"], list):
                    for j in genes["fragment"]:
                        temp_list.append(j["from"])
                        temp_list.append(j["to"])
                else:
                    for k,v in genes["fragment"].items():
                        temp_list.append(v)
                location_dic[i["accession"]][genes["featureId"]] = "|".join(temp_list)


    for record in SeqIO.parse(raw_fasta, "fasta"):
        id = record.id[:record.id.find("|")]
        reference_strain = record.id.split(".")[-1]
        reference_dic[id] = {}
        temp_string = str(record.seq).replace("?","-")
        reference_dic[id]["sequence"] = temp_string.upper()
        for genes in features:
            if genes not in location_dic[reference_strain].keys():
                log_handle.write(record.id + " does not contain " + genes + " region.\n")
                continue
            else:           
                reference_dic[id][genes] = location_dic[reference_strain][genes]

    with open(annotated_json, 'w') as f:
        json_dumps_str = json.dumps(reference_dic, indent=4)
        print(json_dumps_str, file=f)

    time_ran = dt.datetime.now() - time_start
    print("Time Lapse:", time_ran.total_seconds() / 60, "Minutes")

    close_handle(log_handle)
