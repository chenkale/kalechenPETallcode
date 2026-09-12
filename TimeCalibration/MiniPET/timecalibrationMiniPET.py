# Time calibration script

'''
Notes for use:
- Paths should be changed as needed, geo_offsets_dir should be on the Box.
- Make sure row length matches the setting used in compressor.py
- Script works as of 01/20/2026
'''


# File containing calibration data, make sure to have run through compressor.py
calibration_data_dir = '/home/kale-chen/Documents/PET/MiniPET/GAGG/LYSO60min.bin'

# Directory to output calibration tsvs to
output_dir = '/home/kale-chen/Documents/PET/TimeCalibration/MiniPET'

# Start from preexisting iteration, leave as '' if start from iteration 0
#start_from_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\tsvs\time_offset_calibrationit10.tsv'
#start_from_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\tsvs\time_offset_calibrationit5.tsv'
start_from_dir = ''
start_from = 0

chunksize = 128 # Chunk size to calibrate in, 128 is stable, adjust according to RAM, must be a factor of 128
numiterations = 20 # Number of calibration iterations
row_length = 3 # Number of items per row


'''BEGIN SCRIPT'''

import numpy as np
import os
from tqdm import tqdm
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# Reading from photopeak cut data
num_rows = os.path.getsize(calibration_data_dir) // (row_length * 2)      # 6 = 3 (items per row) * 2 (size of int16)
data = np.memmap(calibration_data_dir, dtype = np.int16, mode = 'r', shape=(num_rows, row_length))

# Event counts
eventcounts = np.zeros((128 * 128), dtype = np.int32)
for i in tqdm(range(0, num_rows, 1000000), desc='Event counts'):
    chunk = data[i:i+1000000]
    ids = np.asarray(chunk[:, [0, 1]], dtype = np.int32)
    indices = ids[:, 0] * 128 + ids[:, 1]
    np.add.at(eventcounts, indices, 1)
eventcounts = eventcounts.reshape((128, 128))

# For splitting into double peak
half1=np.asarray(list(range(0, 27)) + list(range(28, 33)) + list(range(64, 73)) + list(range(76, 97)) + [74, 102])
def asichalf1(events):
    return events[np.isin(events[:, 0] % 128, half1)]
def asichalf2(events):
    return events[~np.isin(events[:, 0] % 128, half1)]

# For fitting peaks, gaussian + background
def gaussianbg(x, a, mu, c, bg):
    return a * np.exp(-(x - mu)**2 / (2 * c**2)) + bg

# For histogramming
bins = np.linspace(-2500, 2500, 200) # Consistent histogram bins
bincenters = [(bins[n] + bins[n + 1]) / 2 for n in range(len(bins) - 1)]

# Create lookup tables, either from scratch or from a previously completed iteration
if start_from_dir == '':
    lookupl = np.zeros(128, dtype = np.float32)
    lookupr = np.zeros(128, dtype = np.float32)
else:
    lookup = np.genfromtxt(start_from_dir, delimiter = '\t')[:,4]
    lookupl = np.asarray(lookup[128:256], dtype = np.float32)
    lookupr = np.asarray(lookup[0:128], dtype = np.float32)


'''MAIN LOOP'''
for iteration in range(numiterations):

    # Loop
    for chunknum in tqdm(range(128 // chunksize), desc = f'IDL'):
        chunkidls = range(chunknum * chunksize, (chunknum + 1) * chunksize)
        backadj = chunknum * chunksize

        chunkchannelevents = [np.zeros((np.sum(eventcounts[idl]), 2), dtype = np.float32) for idl in chunkidls]
        writeheads = np.zeros(chunksize, dtype=np.int32)

        # Collect events for idls in chunkidls
        for i in tqdm(range(0, num_rows, 1000000), desc = 'Read chunk', leave = False):
            for idl in chunkidls:
                chunk = data[i:i+1000000] 
                chunk = chunk[chunk[:,0] == idl]
                chunk = chunk[:, [1, 2]]
                chunkchannelevents[idl - backadj][writeheads[idl - backadj]:writeheads[idl - backadj] + len(chunk)] = chunk
                writeheads[idl - backadj] += len(chunk)

        # Apply geometric offset and lookup values
        for idl in tqdm(chunkidls, desc = 'Geo offsets', leave = False):
            idrs = np.asarray(chunkchannelevents[idl - backadj][:,0], dtype = np.int16)
            timediffs = np.asarray(chunkchannelevents[idl - backadj][:,1], dtype = np.float32)
            timediffs = timediffs - lookupl[int(idl)] + lookupr[idrs]
            chunkchannelevents[idl - backadj][:,1] = timediffs

        # Histogram to find time offsets per channel
        for idl in tqdm(chunkidls, desc = 'Hist', leave = False):
            try:
                hist1 = np.histogram(asichalf1(chunkchannelevents[idl - backadj])[:,1], bins = bins)[0]
                popt1, pcov1 = curve_fit(gaussianbg, bincenters, hist1, p0 = [max(hist1), bincenters[list(hist1).index(max(hist1))], 500, 10])
                hist2 = np.histogram(asichalf2(chunkchannelevents[idl - backadj])[:,1], bins = bins)[0]
                popt2, pcov2 = curve_fit(gaussianbg, bincenters, hist2, p0 = [max(hist2), bincenters[list(hist2).index(max(hist2))], 500, 10])
                lookupl[idl] = lookupl[idl] + (popt1[1] + popt2[1]) / 2 / 2 # One of these /2 is to average popt1 and popt2, the other is to split the calibration effect between L and R
            except:
                print(f'Fit failed {idl}')
        #break 

    for chunknum in tqdm(range(128 // chunksize), desc = f'IDR Chunk'):
        chunkidrs = range(chunknum * chunksize, (chunknum + 1) * chunksize)
        backadj = chunknum * chunksize

        chunkchannelevents = [np.zeros((np.sum(eventcounts[:, idr]), 2), dtype = np.float32) for idr in chunkidrs]
        writeheads = np.zeros(chunksize, dtype=np.int32)

        # Collect events for idls in chunkidls
        for i in tqdm(range(0, num_rows, 1000000), desc = 'Read chunk', leave = False):
            for idr in chunkidrs:
                chunk = data[i:i+1000000] 
                chunk = chunk[chunk[:,1] == idr] 
                chunk = chunk[:, [0, 2]]
                chunkchannelevents[idr - backadj][writeheads[idr - backadj]:writeheads[idr - backadj] + len(chunk)] = chunk
                writeheads[idr - backadj] += len(chunk)
    
        # Apply geometric offset
        for idr in tqdm(chunkidrs, desc = 'Geo offsets', leave = False):
            idls = np.asarray(chunkchannelevents[idr - backadj][:,0], dtype = np.int16)
            timediffs = np.asarray(chunkchannelevents[idr - backadj][:,1], dtype = np.float32)
            timediffs = timediffs - lookupl[idls] + lookupr[int(idr)]
            chunkchannelevents[idr - backadj][:,1] = timediffs

        # Histogram to find
        for idr in tqdm(chunkidrs, desc = 'Hist', leave = False):
            try:
                hist1 = np.histogram(asichalf1(chunkchannelevents[idr - backadj])[:,1], bins = bins)[0]
                popt1, pcov1 = curve_fit(gaussianbg, bincenters, hist1, p0 = [max(hist1), bincenters[list(hist1).index(max(hist1))], 500, 10])
                hist2 = np.histogram(asichalf2(chunkchannelevents[idr - backadj])[:,1], bins = bins)[0]
                popt2, pcov2 = curve_fit(gaussianbg, bincenters, hist2, p0 = [max(hist2), bincenters[list(hist2).index(max(hist2))], 500, 10])
                lookupr[idr] = lookupr[idr] - (popt1[1] + popt2[1]) / 2 / 2 # One of these /2 is to average popt1 and popt2, the other is to split the calibration effect between L and R
            except:
                print(f'Fit failed {idr}')

    # Prepare table to write to output tsv
    output = pd.DataFrame()

    for r in range(128):
        output = pd.concat([output, pd.DataFrame([[0, r // 1024, (r % 1024) // 64, r % 64, lookupr[r]]])], ignore_index=True)
    for l in range(128):
        output = pd.concat([output, pd.DataFrame([[1, l // 1024, (l % 1024) // 64, l % 64, lookupl[l]]])], ignore_index=True)

    # Write to tsv
    output.to_csv(os.path.join(output_dir, 'tsvs',f'time_offset_calibrationit{str(iteration + 1 + start_from)}.tsv'), sep='\t', header=None, float_format = '%.7e', index=False)
    print(f'TSV written iteration {iteration}')
    #break

