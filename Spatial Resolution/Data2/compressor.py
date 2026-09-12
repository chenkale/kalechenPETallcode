# Data file compressor for time calibration, also cuts on photopeak


'''
Notes for use:
- Paths should be changed as needed, photopeak bounds should be in the same directory as the script.
- Place coinc.dat file in the same directory as the script, specify name in data_name variable.
- Script works as of 01/20/2026
'''
data_name = '2532.dat'
cut_on_photopeak = False # If false, row_length = 6
make_histogram = False
keep_timeL = False # If set to True, row_length = 4, otherwise row_length = 3
only_times = False # If True, only save TimeL in int64
start_time = 0

photopeak_dir = data_name.replace('.dat', '_ppbounds.txt')

'''BEGIN SCRIPT'''


import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

current_dir = os.path.dirname(os.path.abspath(__file__))

calibration_data_file = os.path.join(current_dir, data_name)
if cut_on_photopeak:
    out_data_file = os.path.join(current_dir, data_name.replace('.dat', '_pcut.bin'))
elif only_times:
    out_data_file = os.path.join(current_dir, data_name.replace('.dat', '_times.bin'))
else:
    out_data_file = os.path.join(current_dir, data_name.replace('.dat', '.bin'))
if start_time != 0:
    out_data_file = os.path.join(current_dir, out_data_file.replace('.bin', f'_start{start_time}.bin'))
# If cutting on photopeak, need to specify:
photopeak_dir = os.path.join(current_dir, photopeak_dir)


if make_histogram:
    lor_pops = [0 for lr in range(3072 * 3072)]
if cut_on_photopeak:
    photopeakL = np.genfromtxt(photopeak_dir, dtype = np.float32)[3072:] * 100
    photopeakR = np.genfromtxt(photopeak_dir, dtype = np.float32)[:3072] * 100

if os.path.exists(out_data_file):
    os.remove(out_data_file)

# To binary with dtype int16, keep only needed data: EnergyL, IDL, EnergyR, IDR, TimeDiff
# Performs time difference calculation, IDL and IDR mod conversions
# int16 is exact for IDL, IDR, and TimeDiff;
# EnergyL and EnergyR are multiplied by 100 to keep (essentially) 2 decimal precision, exact precision is not needed for photopeak cuts
chunk_count = 0
event_count = 0
time0 = -2
for chunk in tqdm(pd.read_csv(calibration_data_file, sep = '\t', chunksize = 10000000, usecols = [2, 3, 4, 7, 8, 9], header = None)):
    chunk.columns = ['TimeL', 'EnergyL', 'IDL', 'TimeR', 'EnergyR', 'IDR']
    chunk['TimeDiff'] = chunk['TimeL'] - chunk['TimeR']     # Calculate time difference, smaller number than individual times
    chunk = chunk[['TimeL','EnergyL', 'IDL', 'EnergyR', 'IDR', 'TimeDiff']]
    if time0 == -1:
        time0 = chunk.iloc[0,0]
    if 'BtF' in calibration_data_file:
        chunk['TimeL'] = 1800 - ((chunk['TimeL'] - time0) // 1000000000000)
    else:
        if not only_times:
            chunk['TimeL'] = (chunk['TimeL'] - time0) // 1000000000000

    chunk = chunk[chunk['TimeL'] > start_time]
    
    chunk['IDL'] = (chunk['IDL'] - 131072) % 3072           # Convert IDs to 0-3071 range for ease of handling, convert back when writing tsv
    chunk['IDR'] = chunk['IDR'] % 3072                      
    chunk['EnergyL'] = chunk['EnergyL'] * 100               # Multiply Energy by 100 to allow for storage as int16 with 2 decimal precision
    chunk['EnergyR'] = chunk['EnergyR'] * 100

    if cut_on_photopeak:
        energyL = chunk["EnergyL"].to_numpy(dtype=float)
        energyR = chunk["EnergyR"].to_numpy(dtype=float)
        lowerboundLs = photopeakL[chunk["IDL"], 0]
        upperboundLs = photopeakL[chunk["IDL"], 1]
        lowerboundRs = photopeakR[chunk["IDR"], 0]
        upperboundRs = photopeakR[chunk["IDR"], 1]
        #mask = ((energyL > lowerboundLs) & (energyL < upperboundLs)
        #    & (energyR > lowerboundRs) & (energyR < upperboundRs))
        mask = ((energyL > lowerboundLs)
            & (energyR > lowerboundRs))
        chunk = chunk[mask]
        if keep_timeL:
            chunk = chunk[['IDL', 'IDR', 'TimeDiff', 'TimeL']]
        else:
            chunk = chunk[['IDL', 'IDR', 'TimeDiff']]
    else:
        chunk = chunk[['EnergyL', 'IDL', 'EnergyR' ,'IDR', 'TimeDiff', 'TimeL']]

    if only_times:
        chunk = chunk[['TimeL']]
        chunk = np.asarray(chunk, dtype = np.int64)
    else:
        chunk = np.asarray(chunk, dtype = np.int16)
    
    with open(out_data_file, 'ab') as f:
        chunk.tofile(f)
    if make_histogram:
        for row in chunk:
            lor_pops[int(row[1] * 3072 + row[3])] += 1
    chunk_count += 1
    event_count += len(chunk)

print('Total events: ' + str(event_count))
# Histogram of LOR populations
if make_histogram:
    print('Total events: ' + str(sum(lor_pops)))
    plt.hist(lor_pops, bins = range(0, 250, 1), color = 'blue', alpha = 0.4)
    plt.yscale('log')
    plt.title('LOR Populations')
    plt.xlabel('Pop')
    plt.ylabel('Counts')
    plt.grid()
    plt.show()
