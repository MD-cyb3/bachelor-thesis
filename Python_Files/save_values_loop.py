# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 17:43:42 2018

@author: mdyck
"""
import csv
import os

def append_values(base_folder, folder_name, file_name, values, header):

    fullfile = os.path.join(base_folder, folder_name, file_name)
    directory = os.path.dirname(fullfile)
    if not os.path.exists(directory):
        os.makedirs(directory)    

    if os.path.exists(fullfile):
        append_write = 'a' # append if already exists
    else:
        append_write = 'w' # make a new file if not

    with open(fullfile, append_write) as f:
        writer = csv.writer(f)
        writer.writerow([header, values])

