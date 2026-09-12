# Photopeak cut finder and cutter
import numpy as np
import os
from tqdm import tqdm
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches
from scipy.optimize import curve_fit
#import pandas as pd

names = ['hcenter', 'mcenter', 'hcirc', 'hvert', 'pointsr', 'pointsR', 'htof']

data_dir = '/home/kale-chen/Documents/PET/TPPT2026/data/'
chunksize = 64
cut_range = 2.5

use_previous_photopeaks = False

def countidl(idl):
    return np.sum(eventcounts[idl])
def countidr(idr):
    return np.sum(eventcounts[:,idr])

bins = np.linspace(1500, 3500, 80)
bin_centers = [(bins[q] + bins[q+1])/2 for q in range(0,len(bins)-1)]
def gaussian(x, a, mu, c):
    return a * np.exp(-(x - mu)**2 / (2 * c**2))
# Function for finding photopeak, modified from ChannelPairBuilder.py
def findPhotopeak(channelevents, id):
    values = np.histogram(channelevents, bins = bins)[0] # values are to be fitted
    # Find the photopeak within 3 bins starting from the RHS of the graph
    photopeak = 0
    for i in range(1, len(values)-4):
        if values[-1*i] > values[-1*(i+1)] and values[-1*i] > values[-1*(i+2)] and values[-1*i] > values[-1*(i+3)] and values[-1*(i+1)] != 0 and values[-1*i] >= 0.50 * np.max(values):
            photopeak = len(values)-1*i
            break
    fit_p = [0,1,1]
    fit_x = list(bin_centers[photopeak-10:photopeak+20])
    fit_y = list(values[photopeak-10:photopeak+20])

    if len(fit_x) != 0 and len(fit_y) != 0:
        try:
            fit_p, fit_co = curve_fit(gaussian, fit_x, fit_y, p0=[max(fit_y), fit_x[fit_y.index(max(fit_y))], 100])
        except:
            print(f'Photopeak fit failed {id}')
    return fit_p

for name in names:
    print(f'Processing {name}')
    txt_out_dir = os.path.join(data_dir, name + '_ppbounds.txt')
    num_rows = os.path.getsize(os.path.join(data_dir, name + '.bin')) // 12
    data = np.memmap(os.path.join(data_dir, name + '.bin'), dtype = np.int16, mode = 'r', shape=(num_rows, 6))
    print(num_rows)
    if not use_previous_photopeaks:
        # Event counts
        eventcounts = np.zeros((3072 * 3072), dtype = np.int32)
        for i in tqdm(range(0, num_rows, 1000000), desc='Event counts'):
            chunk = data[i:i+1000000]
            ids = np.asarray(chunk[:, [1, 3]], dtype = np.int32)
            indices = ids[:, 0] * 3072 + ids[:, 1]
            np.add.at(eventcounts, indices, 1)
        eventcounts = eventcounts.reshape((3072, 3072))

        photopeakboundsL = np.zeros((3072, 2), dtype = np.float32)
        photopeakboundsR = np.zeros((3072, 2), dtype = np.float32)

        for chunknum in tqdm(range(3072 // chunksize), desc = 'ID Chunk'):
            chunkidls = range(chunknum * chunksize, (chunknum + 1) * chunksize)
            backadj = chunknum * chunksize

            chunkchannelevents = [np.zeros(countidl(chunkidls[idl - backadj]), dtype = np.int16) for idl in chunkidls]
            writeheads = np.zeros(chunksize, dtype=np.int32)

            # Collect events for idls in chunkidls
            for i in tqdm(range(0, num_rows, 1000000), desc = 'Read chunk', leave = False):
                for idl in chunkidls:
                    chunk = data[i:i+1000000] 
                    chunk = chunk[chunk[:,1] == idl] # Get chunk with only idl
                    chunk = chunk[:, 0] # EnergyL
                    chunkchannelevents[idl - backadj][writeheads[idl - backadj]:writeheads[idl - backadj] + len(chunk)] = chunk
                    writeheads[idl - backadj] += len(chunk)
            
            # Find photopeaks for idls in chunkidls 
            for idl in chunkidls:
                fit = findPhotopeak(chunkchannelevents[idl - backadj], idl)
                if fit[0] != 0:
                    photopeakboundsL[idl][0] = fit[1] - cut_range * fit[2]
                    photopeakboundsL[idl][1] = fit[1] + cut_range * fit[2]
                if idl % 100 == 0:
                    #print(f'\nIDL: {idl}, photopeak bounds: {photopeakboundsL[idl]}')
                    pass

        for chunknum in tqdm(range(3072 // chunksize), desc = 'ID Chunk'):
            chunkidrs = range(chunknum * chunksize, (chunknum + 1) * chunksize)
            backadj = chunknum * chunksize

            chunkchannelevents = [np.zeros(countidr(chunkidrs[idr - backadj]), dtype = np.int16) for idr in chunkidrs]
            writeheads = np.zeros(chunksize, dtype=np.int32)

            # Collect events for idrs in chunkidrs
            for i in tqdm(range(0, num_rows, 1000000), desc = 'Read chunk', leave = False):
                for idr in chunkidrs:
                    chunk = data[i:i+1000000] 
                    chunk = chunk[chunk[:,3] == idr] # Get chunk with only idr
                    chunk = chunk[:, 2] # EnergyR
                    chunkchannelevents[idr - backadj][writeheads[idr - backadj]:writeheads[idr - backadj] + len(chunk)] = chunk
                    writeheads[idr - backadj] += len(chunk)
            
            # Find photopeaks for idrs in chunkidrs 
            for idr in chunkidrs:
                fit = findPhotopeak(chunkchannelevents[idr - backadj], idr)
                if fit[0] != 0: 
                    photopeakboundsR[idr][0] = fit[1] - cut_range * fit[2]
                    photopeakboundsR[idr][1] = fit[1] + cut_range * fit[2]
                if idr % 100 == 0:
                    #print(f'\nIDR: {idr}, photopeak bounds: {photopeakboundsR[idr]}')
                    pass

        with open(txt_out_dir, 'w') as f:
            for row in photopeakboundsR:
                f.write(f'{row[0] / 100}\t{row[1] / 100}\n')
            for row in photopeakboundsL:
                f.write(f'{row[0] / 100}\t{row[1] / 100}\n')

    if not use_previous_photopeaks:
        photopeaksL = photopeakboundsL[:, 0]
        photopeaksR = photopeakboundsR[:, 0]
    else:
        photopeaksL = np.loadtxt(os.path.join(data_dir, name + '_ppbounds.txt'))[3072:, 0] * 100
        photopeaksR = np.loadtxt(os.path.join(data_dir, name + '_ppbounds.txt'))[:3072, 0] * 100
    print(photopeaksL[0:5])
    print(photopeaksR[0:5])
    # Cut data
    bin_out_dir = os.path.join(data_dir, name + 'p.bin')
    mask_out_dir = os.path.join(data_dir, name + '_pcut_mask.bin') # binary mask of the cut
    if os.path.exists(bin_out_dir):
        os.remove(bin_out_dir)
    if os.path.exists(mask_out_dir):
        os.remove(mask_out_dir)
    keptevents = 0
    for i in tqdm(range(0, num_rows, 1000000), desc = 'Cut chunk'):
        chunk = data[i:i+1000000]
        mask = np.zeros(len(chunk), dtype = np.bool)
        mask = ((chunk[:,0] > photopeaksL[chunk[:,1]]) & (chunk[:,2] > photopeaksR[chunk[:,3]]))
        chunk = chunk[mask]
        with open(bin_out_dir, 'ab') as f:
            chunk.tofile(f)
        with open(mask_out_dir, 'ab') as f:
            mask.tofile(f)
        keptevents += len(chunk)
    print(f'Data cut and saved to {bin_out_dir} and {mask_out_dir}')
    print(f'Kept {keptevents} events')