# Time calibration script

''''''
# File containing calibration data, make sure to have run through compressor.py
calibration_data_dir = '/home/kale-chen/Documents/PET/TimeCalibration/LineSourceStudy/remapping/center.bin'

# File containing geometric offsets table, can generate from map with geotabler.py
geo_offsets_dir = '/home/kale-chen/Documents/PET/TimeCalibration/ToUpload/geometrictimeoffsets2.bin'

# Directory to output calibration tsvs to
output_dir = f'/home/kale-chen/Documents/PET/TimeCalibration/LineSourceStudy/remapping/'

chunksize = 128 # Chunk size to calibrate in, 128 is stable, adjust according to RAM, must be a factor of 3072

''''''

import numpy as np
import os
from tqdm import tqdm
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#import pandas as pd


# Reading from photopeak cut data
num_rows = os.path.getsize(calibration_data_dir) // 6
data = np.memmap(calibration_data_dir, dtype = np.int16, mode = 'r', shape=(num_rows, 3))

# Event counts
eventcounts = np.zeros((3072 * 3072), dtype = np.int32) 
for i in tqdm(range(0, num_rows, 1000000), desc='Event counts'):
    chunk = data[i:i+1000000]
    ids = np.asarray(chunk[:, [0, 1]], dtype = np.int32)
    indices = ids[:, 0] * 3072 + ids[:, 1]
    np.add.at(eventcounts, indices, 1)
eventcounts = eventcounts.reshape((3072, 3072))

# For applying geometric offset
geooffsets = np.fromfile(geo_offsets_dir, dtype = np.float32).reshape((3072, 3072))

# Lookup table
lookup = np.zeros((3072, 3072), dtype = np.int32)

for chunknum in tqdm(range(3072 // chunksize), desc = f'Big loop'):
    chunkidls = range(chunknum * chunksize, (chunknum + 1) * chunksize)
    backadj = chunknum * chunksize

    chunkchannelevents = [np.zeros((np.sum(eventcounts[idl]), 2), dtype = np.int16) for idl in chunkidls]
    writeheads = np.zeros(chunksize, dtype=np.int32)

    # Collect events for idls in chunkidls
    for i in tqdm(range(0, num_rows, 1000000), desc = 'Read chunk', leave = False):
        for idl in chunkidls:
            chunk = data[i:i+1000000] 
            chunk = chunk[chunk[:,0] == idl]
            chunk = chunk[:, [1, 2]]
            chunkchannelevents[idl - backadj][writeheads[idl - backadj]:writeheads[idl - backadj] + len(chunk)] = chunk
            writeheads[idl - backadj] += len(chunk)

    # Collect events for idrs
    for idl in tqdm(chunkidls, desc = 'IDR loop', leave = False):
        channelevents = chunkchannelevents[idl - backadj]
        for idr in range(3072):
            timediffs = channelevents[channelevents[:,0] == idr][:,1]
            try:
                lookup[idl][idr] = -np.mean(timediffs) + geooffsets[idl][idr]
            except:
                pass

lookup.tofile(output_dir + 'hailmaryfullgeo.bin')


