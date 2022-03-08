#!/usr/bin/env python3

"""
Name: reference_matcher.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Map sequence to the closest reference sequence list based on mini-map2. 

Options:
    :param in_fasta: Raw sequences needed to be referenced to reference list in fasta format (Required)
    :param reference_sequence: Reference list in fasta format (Required)
    :param out_fasta: Output list of sequences referenced (Default: referenced.fasta)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

import mappy as mp
from Bio import SeqIO
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import datetime as dt
from genofunc.utils import *

def reference_matcher(in_fasta,reference_sequence,out_fasta,log_file):
    time_start = dt.datetime.now()

    if not file_check(in_fasta):
        print("Input file does not exist. Please check the file path.")
        sys.exit()
    if not file_check(reference_sequence):
        print("Reference file does not exist. Please check the file path.")
        sys.exit()    

    outfile = open(out_fasta,"w")
    log_handle = get_log_handle(log_file)
    log_handle.write("Strain\tReference_Sequence\tMatchingBase\tNumberofMismatches\tCigarStart\tCigarArray\n")

    reference_list = mp.Aligner(reference_sequence)

    for record in SeqIO.parse(in_fasta, "fasta"):
        matching_length = 0
        mismatches = 8000
        for hit in reference_list.map(record.seq):
            if hit.mlen >= matching_length and hit.NM <= mismatches:
                matching_length = hit.mlen
                mismatches = hit.NM
                reference_id = hit.ctg
                start_ref = hit.r_st
                end_ref = hit.r_en
                cigar_array = hit.cigar
            new_id = "|".join([record.id,str(matching_length),reference_id])
            new_record = SeqRecord(record.seq,new_id,description="")
        SeqIO.write(new_record, outfile, "fasta-2line")
        log_handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.id,hit.ctg,matching_length,mismatches,start_ref,end_ref,cigar_array))
    
    time_ran = dt.datetime.now() - time_start
    print("Time Lapse:", time_ran.total_seconds() / 60, "Minutes")

    close_handle(log_handle)
