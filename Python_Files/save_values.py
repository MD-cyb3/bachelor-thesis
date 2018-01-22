# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 17:43:42 2018

@author: mdyck
"""
import numpy as np
import csv

def save_values_to_file(file_name, x, u, m1, abs_m1, theta, n, p):
    name_of_file = 'csv_files_simulation_study/' + file_name + '.csv'
    headers = [['x'], ['u'], ['m1'], ['abs_m1'], ['theta'], ['n'], ['p']]

    values = [headers[0], x, headers[1], u, headers[2], m1, headers[3], abs_m1,
              headers[4], theta, headers[5], n, headers[6], p]
    
    with open(name_of_file, 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(values)