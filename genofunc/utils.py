"""
Name: utils.py
Author: Rachel Colquhoun & Xiaoyu Yu
Date: 01 April 2022
Description: Utility functions used within genofunc

This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import sys
import os

def get_log_handle(log_file):
    if log_file:
        log_handle = open(log_file,"w")
    else:
        log_handle = sys.stdout
    return log_handle

def get_out_handle(out_fasta):
    if out_fasta == '' or out_fasta == [] or out_fasta == ['']:
        out_fasta = None
    if not out_fasta:
        out_handle = sys.stdout
    else:
        out_handle = open(out_fasta,"w")
    return out_handle

def get_in_handle(in_fasta):
    if not in_fasta:
        in_handle = sys.stdin
    else:
        in_handle = open(in_fasta)
    return in_handle

def close_handle(handle):
    if handle:
        handle.close()

def file_check(file_name):
    if os.path.exists(file_name):
        return True
    else:
        return False