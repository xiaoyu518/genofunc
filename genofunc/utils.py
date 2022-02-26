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

def close_handle(handle):
    if handle:
        handle.close()

def file_check(file_name):
    if os.path.exists(file_name):
        return True
    else:
        return False
