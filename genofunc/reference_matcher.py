#!/usr/bin/env python3

"""
Name: reference_matcher.py
Author: Xiaoyu Yu
Date: 16 April 2023
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
        mismatches = round(float(percentage_match)*len(str(record.seq)),0)
        if virus_type.upper() == "DNA":
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
                log_handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.id,hit.ctg,matching_length,mismatches,start_ref,end_ref,cigar_array))
                if matching_length < mismatches:
                    temp_seq = record.seq.reverse_complement()
                    for hit in reference_list.map(temp_seq):
                        if hit.mlen >= matching_length and hit.NM <= mismatches:
                            reverse_matching_length = hit.mlen
                            reverse_mismatches = hit.NM
                            reverse_reference_id = hit.ctg
                            reverse_start_ref = hit.r_st
                            reverse_end_ref = hit.r_en
                            reverse_cigar_array = hit.cigar
                    if reverse_matching_length > matching_length:
                        new_id = "|".join([record.id,str(reverse_matching_length),reverse_reference_id])
                        new_record = SeqRecord(temp_seq,new_id,description="")
                        log_handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.id,hit.ctg,reverse_matching_length,reverse_mismatches,reverse_start_ref,
                                                                            reverse_end_ref,reverse_cigar_array))                       
            if matching_length != 0 or reverse_matching_length != 0:
                SeqIO.write(new_record, outfile, "fasta-2line")
            
        elif virus_type.upper() == "RNA":
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
            if matching_length != 0:
                SeqIO.write(new_record, outfile, "fasta-2line")
            log_handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.id,hit.ctg,matching_length,mismatches,start_ref,end_ref,cigar_array))
        else:
            print("Virus type not valid, please re-enter correct virus type option (RNA/DNA). Program Exiting")
            sys.exit()
    
    time_ran = dt.datetime.now() - time_start
    print("Time Lapse:", time_ran.total_seconds() / 60, "Minutes")

    close_handle(log_handle)
