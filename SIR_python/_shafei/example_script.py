# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 20:02:19 2023

@author: Vincent Bazinet
"""

import os
import h5py
import simulated_atrophy as sim
import numpy as np
from scipy import stats
from statistics import mean

import time

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

def pearson_corr(array_2d, array_1d):
    # Standardize the arrays
    array_1d_mean = np.mean(array_1d)
    array_1d_std = np.std(array_1d)

    # Standardize 2D array
    array_2d_mean = np.mean(array_2d, axis=0)
    array_2d_std = np.std(array_2d, axis=0)

    # Compute the correlations
    correlations = np.dot((array_2d - array_2d_mean), (array_1d - array_1d_mean)) / (array_1d_std * array_2d_std)

    # Find the maximum correlation and its index
    max_corr = np.max(correlations)
    max_corr_index = np.argmax(correlations)

    print("Maximal correlation:", max_corr)
    print("Index of maximal correlation in 2D array:", max_corr_index)

    return max_corr


# os.chdir(os.path.expanduser("~") +
#          "/Dropbox/Sid/R_files/SV2A_pet/scripts/SIR_code/")
#
# import simulated_atrophy as sim
# import pickle
#
# # Load data for the example
# with open("data.pickle", 'rb') as handle:
#     data = pickle.load(handle)
#
# # Simulate atrophy:
# atrophy = sim.simulate_atrophy(data['SC_den'], data['SC_len'], 109, data['ROI_size'])

# Load workspace variables from HDF5 file
ws_file = f'../data/workspace_S132.h5'

with h5py.File(ws_file, 'r') as f:
    gene_expr = f['gene_expr'][:]
    ifod_len_35 = f['ifod_len_35'][:]
    ifod_den_35 = f['ifod_den_35'][:]
    roi_size = f['ROIsize'][:]
    bgs = f['bgs'][:]
    cobre = f['cobre'][:]
    hcpep = f['hcpep'][:]
    stages = f['stages'][:]

# reshape files
roi_size = roi_size.squeeze()
bgs = bgs.squeeze()
cobre = cobre.squeeze()
hcpep = hcpep.squeeze()
stages = stages.squeeze()

tic()
for i in range(20):
    gba = np.random.normal(0, 1, 66)
    snca = np.random.normal(0, 1, 66)
    atrophy = sim.simulate_atrophy(ifod_den_35, ifod_len_35, 41, np.transpose(roi_size), GBA=gba, SNCA = snca)
toc()

# corr_bgs = pearson_corr(atrophy, bgs)
# corr_cobre = pearson_corr(atrophy, cobre)
# corr_hcpep = pearson_corr(atrophy, hcpep)
# corr_stages = pearson_corr(atrophy, stages)
# corr_mean = mean([corr_bgs, corr_cobre, corr_hcpep])
#
# print(corr_mean, corr_stages, corr_bgs, corr_cobre, corr_hcpep)
#
