# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 17:43:42 2018

@author: mdyck
"""
import numpy as np
import csv
import os

def save_time_dependent_values(folder_name ,file_name, u, m1, abs_m1):
    name_of_folder = 'csv_files_simulation_SECONDstudy/' + folder_name + '/'
    name_of_path = name_of_folder + file_name + '.csv'

    directory = os.path.dirname(name_of_path)
    if not os.path.exists(directory):
        os.makedirs(directory)    
    
    headers = [['u'], ['m1'], ['abs_m1']]

    values = zip(u, m1, abs_m1)
    
    with open(name_of_path, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(values)
        
def save_distribution_values(folder_name ,file_name, x, n_begin, n_end, p_begin, p_end):
    name_of_folder = 'csv_files_simulation_SECONDstudy/' + folder_name + '/'
    name_of_path = name_of_folder + file_name + '.csv'

    directory = os.path.dirname(name_of_path)
    if not os.path.exists(directory):
        os.makedirs(directory)    
        
    headers = [['x'], ['n(x,0)'], ['n(x,t_end)'], ['p(x,0)'], ['p(x,t_end)']]

    values = zip(x, n_begin, n_end, p_begin, p_end)
    
    with open(name_of_path, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(values)
        
def save_phase_dependent_values(folder_name, file_name, theta_begin, theta_end):
    name_of_folder = 'csv_files_simulation_SECONDstudy/' + folder_name + '/'
    name_of_path = name_of_folder + file_name + '.csv'

    directory = os.path.dirname(name_of_path)
    if not os.path.exists(directory):
        os.makedirs(directory)    
        
    headers = [['theta(0)'], ['theta(t_end)']]

    values = zip(theta_begin, theta_end)
    
    with open(name_of_path, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(values)
        
def save_information_file(folder_name ,file_name, info_text):
    name_of_folder = 'csv_files_simulation_SECONDstudy/' + folder_name + '/'
    name_of_path = name_of_folder + file_name + '.txt'

    directory = os.path.dirname(name_of_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    with open(name_of_path, 'wb') as f:
        f.write(info_text)
