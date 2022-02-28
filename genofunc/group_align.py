#!/usr/bin/env python3

"""
Name: group_align.py
Author: Xiaoyu Yu
Date: 01 April 2022
Description: Split the fasta file into groups of sequences set by a user threshold and align them in groups against reference. post group aligned sequences will be concatenated into a single alignment. 

Options:
    :param in_dir: Input directory (Required)
    :param group_size: Group size for fasta file to split by (Required)
    :param reference_dir: Reference sequence directory. Reference sequence used based on matching sequence and reference file names (Required)
    :param out_dir: Output folder directory (Default: ./)

This file is part of PANGEA HIV project (www.pangea-hiv.org).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from Bio import SeqIO
import datetime as dt
import glob
import os
import sys
import subprocess as sp
from genofunc.utils import *

def group_align(in_dir,group_size,reference_dir,out_dir,log_file):
    log_handle = get_log_handle(log_file)
    time_start = dt.datetime.now()
    fasta_files = glob.glob(in_dir + "*.fasta")
    ref_files = glob.glob(reference_dir + "*.fasta")
    ref_table = {"gag":"gag_ref.fasta",
                 "pol":"pol_ref.fasta",
                 "vif":"vif_ref.fasta",
                 "vpr":"vpr_ref.fasta",
                 "tat":"tat_ref.fasta",
                 "vpu":"vpu_ref.fasta",
                 "env":"env_ref.fasta",
                 "nef":"nef_ref.fasta"
                 }

    if fasta_files == []:
        print("No fasta files in input directory. Please re-enter new one. Program Exiting.")
        sys.exit()
    if ref_files == []:
        print("No reference files in input directory. Please re-enter new one. Program Exiting.")
        sys.exit()

    for files in fasta_files:
        file_counter = 1
        seq_counter = 0
        out_file = out_dir + "out_fasta_" + str(file_counter) + ".fasta"
        output_fasta = open(out_file,"w")
        for gene in ref_table.keys():
            if files.find(gene) > -1:
                ref_seq = reference_dir + ref_table[gene]
                file_name = gene + "_alignment.fasta"
                break
        for record in SeqIO.parse(files, "fasta"):    
            if seq_counter == group_size:
                log_handle.write("File " + out_file + " has been created with " + str(seq_counter) + " sequences.\n")
                seq_counter = 0
                file_counter += 1
                output_fasta.close()            
                align_file = out_file[:-6] + "_aligned.fasta"
                align_args = ["augur","align","--sequences",out_file,"--reference-sequence",ref_seq,"--remove-reference","--output",align_file]
                log_handle.write("File " + align_file + " has been created with reference sequence " + ref_seq + ".\n")
                alignment = sp.Popen(align_args)
                alignment.wait()            
                os.remove(out_file)
                out_file = out_dir + "out_fasta_" + str(file_counter) + ".fasta"
                output_fasta = open(out_file,"w")
            seq_counter += 1
            SeqIO.write(record, output_fasta, "fasta-2line")

        #For the last file with sequences not reaching group size
        log_handle.write("File " + out_file + " has been created with " + str(seq_counter) + " sequences.\n")
        output_fasta.close()
        align_file = out_file[:-6] + "_aligned.fasta"
        align_args = ["augur","align","--sequences",out_file,"--reference-sequence",ref_seq,"--remove-reference","--output",align_file]
        log_handle.write("File " + align_file + " has been created with reference sequence " + ref_seq + ".\n")
        alignment = sp.Popen(align_args)
        alignment.wait()
        os.remove(out_file)  

        cat_file = out_dir + os.path.basename(files)[:-6] + "_alignment.fasta"
        cat_command = "cat " + out_dir + "*_aligned.fasta > " + cat_file
        os.system(cat_command)
        log_handle.write("File " + cat_file + " has been created with all groups from " + files +".\n")

        remove_file = glob.glob(out_dir + "out_*") 
        for removal in remove_file:
            os.remove(removal)
            log_handle.write("File " + removal + " has been removed.\n")

        remove_file = glob.glob(out_dir + "*_aligned.fasta") 
        for removal in remove_file:
            os.remove(removal)
            log_handle.write("File " + removal + " has been removed.\n")

        log_handle.write("File " + file_name + " has been created without any reference sequences.\n")

    time_ran = dt.datetime.now() - time_start
    print("Time Lapse:", time_ran.total_seconds() / 60, "Minutes")

    close_handle(log_handle)
