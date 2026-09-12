'''
Different approach for organizing analysis. Here a class is created for a "chain" centered around a data file.
All steps are organized as methods in the class. 
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from tqdm import tqdm
import time

class AnalysisChain:

    dat_dir = None
    bin_dir = None
    binpcut_dir = None
    photopeak_dir = None

    data = None
    num_cols = 0
    num_rows = 0

    data_pcut = None
    num_cols_pcut = 0
    num_rows_pcut = 0

    def __init__(self, dat_dir, photopeak_dir):
        self.dat_dir = dat_dir
        self.photopeak_dir = photopeak_dir
        self.bin_dir = self.dat_dir.replace('.dat', '1.bin')
        self.binpcut_dir = self.dat_dir.replace('.dat', '_pcut.bin')

    # Takes 6144 x n array and writes to txt file
    def write_txt(self, array, file_name):
        with open(file_name, 'w') as f:
            for i in range(array.shape[0]):
                for j in range(array.shape[1]):
                    f.write(f'{round(array[i, j], 8)}\t')
                f.write('\n')

    def compress(self, cut_on_photopeak = False, keep_timeL = True):
        chunk_count = 0
        event_count = 0
        time0 = -1
        # Delete previous file if exists
        outbin = self.bin_dir if not cut_on_photopeak else self.binpcut_dir
        if os.path.exists(outbin):
            os.remove(outbin)

        # Setup photopeak boundaries
        if cut_on_photopeak:
            photopeakL = np.genfromtxt(self.photopeak_dir, dtype=np.float32)[3072:] * 100  # 3072 x 2
            photopeakR = np.genfromtxt(self.photopeak_dir, dtype=np.float32)[:3072] * 100  # 3072 x 2

        for chunk in tqdm(pd.read_csv(self.dat_dir, sep = '\t', chunksize = 10000000, usecols = [2, 3, 4, 7, 8, 9], header = None)):
            chunk.columns = ['TimeL', 'EnergyL', 'IDL', 'TimeR', 'EnergyR', 'IDR']
            chunk['TimeDiff'] = chunk['TimeL'] - chunk['TimeR']
            chunk = chunk[['TimeL', 'EnergyL', 'IDL', 'TimeR', 'EnergyR', 'IDR', 'TimeDiff']]

            if time0 == -1:
                time0 = chunk.iloc[0,0]

            # Original code for aligning time calibration data            
            if 'BtF' in self.dat_dir:
                chunk['TimeL'] = 1800 - ((chunk['TimeL'] - time0) // 1000000000000)
            else:
                chunk['TimeL'] = (chunk['TimeL'] - time0) // 1000000000000

            # Convert IDs to 0-3071 range for ease of handling
            chunk['IDL'] = (chunk['IDL'] - 131072) % 3072
            chunk['IDR'] = chunk['IDR'] % 3072 

            # Multiply Energy by 100 to allow for storage as int16 with 2 decimal precision
            chunk['EnergyL'] = chunk['EnergyL'] * 100
            chunk['EnergyR'] = chunk['EnergyR'] * 100

            if cut_on_photopeak:
                energyL = chunk["EnergyL"].to_numpy(dtype=float)
                energyR = chunk["EnergyR"].to_numpy(dtype=float)
                idls = chunk["IDL"].to_numpy(dtype=int)
                idrs = chunk["IDR"].to_numpy(dtype=int)
                lowerboundLs = photopeakL[idls, 0]
                lowerboundRs = photopeakR[idrs, 0]
                mask = ((energyL > lowerboundLs) & (energyR > lowerboundRs))
                chunk = chunk[mask]

                if keep_timeL:
                    chunk = chunk[['IDL', 'IDR', 'TimeDiff', 'TimeL']]
                else:
                    chunk = chunk[['IDL', 'IDR', 'TimeDiff']]
            elif keep_timeL:
                chunk = chunk[['EnergyL', 'IDL', 'EnergyR', 'IDR', 'TimeDiff', 'TimeL']]
            else:
                chunk = chunk[['EnergyL', 'IDL', 'EnergyR', 'IDR', 'TimeDiff']]

            # Write as int16 array
            chunk_np = np.asarray(chunk, dtype=np.int16)
            with open(outbin, 'ab') as f:
                chunk_np.tofile(f)
            
            chunk_count += 1
            event_count += len(chunk_np)

        print(f'Written to {outbin}')
        print(f'Event count: {event_count}')

    def event_counts(self, mode = 'perchannel', num_time_bins = None):
        
        if mode == 'perlor':
            event_counts = np.zeros((3072 * 3072), dtype = np.int32)
            for i in tqdm(range(0, num_rows, 1000000), desc='Event counts'):
                chunk = data[i:i+1000000]
                chunk.columns = ['IDL', 'IDR']
                chunk['IDL'] = chunk['IDL'] + 3072
                np.add.at(event_counts, chunk['IDL'], 1)
                np.add.at(event_counts, chunk['IDR'], 1)
            return event_counts
        elif mode == 'perchannel_timebinned' and num_time_bins is not None:
            event_counts = np.zeros((6144, num_time_bins), dtype = np.int32)
            for i in tqdm(range(0, num_rows, 1000000), desc='Event counts'):
                chunk = data[i:i+1000000]
                chunk.columns = ['IDL', 'IDR']
                chunk['IDL'] = chunk['IDL'] + 3072
                np.add.at(event_counts, chunk['IDL'], 1)
                np.add.at(event_counts, chunk['IDR'], 1)
            return event_counts

def main():
    dat_dir = '/home/kale-chen/Documents/PET/TimeCalibration/LineSourceStudy/MovingSourceRun_vertical_OptimizedV2_OV_T1_2600_TrigOn_coinc.dat'
    ac = AnalysisChain(dat_dir, None)
    ac.compress(keep_timeL = True)

if __name__ == '__main__':
    main()