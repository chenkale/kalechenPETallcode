# Calibration data finalizer
# Selects sample of remaining LORs that could not be filled

data_storage_dir = 'C:/Users/burri/Documents/PET/TimeCalibration/CalibrationDataPreparer/test/'
num_lor_chunks = 8
# Create histogram showing the sizes (how many coincidences) of incomplete LORs
incomplete_lors_sizes_histogram = True
delete_temporary_lor_files = False

import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
import statistics as stat
import random
import time
import matplotlib.pyplot as plt
from tqdm import tqdm

start = time.time() # For keeping track of time

# First, combine all filled LORs into 1 file
for lorchunk in tqdm(range(num_lor_chunks)):
    # Read in incomplete LORs, skip already filled LORs to save computation and storage space
    filled_lors_file = data_storage_dir + 'filled_lors_file_' + str(lorchunk).zfill(2) + '.bin'
    num_rows = os.path.getsize(filled_lors_file) // 4
    if num_rows != 0:
        filled_lors = np.memmap(filled_lors_file, dtype=np.float32, mode='r', shape=(num_rows,))
        for i in range(0, num_rows, 1000000):
            chunk = np.asarray(filled_lors[i:i+1000000]) # Read as chunks
            with open(data_storage_dir + r'\final_lors.bin', 'ab') as f:
                f.write(chunk.tobytes())
        del filled_lors, chunk # Delete memmap and dependent objects after use
        print(f'{round(time.time() - start, 3)} appended filled lors file {lorchunk}')
    else:
        print(f'{round(time.time() - start, 3)} filled lors file {lorchunk} is empty, could not be appended')

# Second, go through remaining incomplete LORs, take sample of LORs that could not be fitted, otherwise attempt fit on lors with insufficient statistics

# Simple gaussian
def gauss(x,A,mu,sigma):
    return A*np.exp(-(x - mu)**2 / (2 * sigma**2))

count_filled = 0
lorchunksize = int(3072 * 3072 / num_lor_chunks)
filled_lor_counts = []
for lorchunk in tqdm(range(num_lor_chunks)):
    # Read incomplete LORs
    incomplete_lors_file = data_storage_dir + 'incomplete_lors_file_' + str(lorchunk).zfill(2) + '.bin'
    coinc_list = [[] for lr in range(lorchunksize)] # idL = lr // 3072, idR = lr % 3072
    if os.path.exists(incomplete_lors_file): 
        count = 0
        num_rows = os.path.getsize(incomplete_lors_file) // 12 # 12 = 4 (size of float32) * 3 (entries per row: timediff, energyL, energyR)
        incomplete_lors = np.memmap(incomplete_lors_file, dtype=np.float32, mode='r', shape=(num_rows,3))
        temp = [] # For storing each LORs events before stored to coinc_list
        lr = 0 # For matching file read to proper LOR
        for i in range(0, num_rows, 1000000):
            chunk = np.asarray(incomplete_lors[i:i+1000000]) # Read as chunks
            for row in chunk:
                if row[0] > 0 and row[0] < 1: # Indicator for new LOR is a non integer
                    if len(temp) != 0:
                        coinc_list[lr] = temp # Write temp to coinc_list[lr]
                        temp = [] # Reset temp
                    lr = int(row[1] * 3072 + row[2]) # Match proper LOR
                else:
                    temp.append(list(row))
                    count += 1
        del incomplete_lors, chunk, row, temp # Delete memmap and dependent objects after use
        print(f'{round(time.time() - start, 3)} read incomplete lors {count}')
    else: 
        print(f'{round(time.time() - start, 3)} read incomplete lors 0')
    
    # Perform fits without consideration for sufficient statistics
    filled_lors = []
    for lr in range(lorchunksize):
        if len(coinc_list[lr]) > 0:
            temp = np.asarray(coinc_list[lr])
            # Photopeak
            valuesL, binsL = np.histogram(temp[:,1], bins = np.linspace(0,45,45)) # valuesL are to be fitted
            valuesR, binsR = np.histogram(temp[:,2], bins = np.linspace(0,45,45)) # valuesR are to be fitted
            # Binning will be same for both so just define bin centers once w.r.t. binsR
            bin_centers = [(binsR[q] + binsR[q+1])/2 for q in range(0,len(binsR)-1)]
            
            # Find the photopeak within 3 bins starting from the RHS of the graph
            for i in range(1, len(valuesR)-4):
                PhotopeakR = 0
                if valuesR[-1*i] > valuesR[-1*(i+1)] and valuesR[-1*i] > valuesR[-1*(i+2)] and valuesR[-1*i] > valuesR[-1*(i+3)] and valuesR[-1*(i+1)] != 0 and valuesR[-1*i] >= 20:
                    PhotopeakR = len(valuesR)-1*i
                    break
            for i in range(1, len(valuesL)-4):
                PhotopeakL = 0
                if valuesL[-1*i] > valuesL[-1*(i+1)] and valuesL[-1*i] > valuesL[-1*(i+2)] and valuesL[-1*i] > valuesL[-1*(i+3)] and valuesL[-1*(i+1)] != 0 and valuesL[-1*i] >= 20:
                    PhotopeakL = len(valuesL)-1*i
                    break
            
            fitR_p = [0,1,1] # Initialized to amplitude=0, mean=1, sigma=1
            fitL_p = [0,1,1]
            # Set the range of x and y values to fit: from 7 bins to the left of the photopeak, up to the end of the spectra
            fitR_x = list(bin_centers[PhotopeakR-4:])
            fitL_x = list(bin_centers[PhotopeakL-4:])
            fitR_y = list(valuesR[PhotopeakR-4:])
            fitL_y = list(valuesL[PhotopeakL-4:])

            # Make sure fit ranges are non-empty, perform the fits, keep parameters
            if len(fitR_y) != 0 and len(fitR_x) != 0:
                try:
                    # curve_fit(function, x_data, y_data, p0=param_guesses, bounds=lower and upper param bounds)
                    fitR_p, fitR_co = curve_fit(gauss, fitR_x, fitR_y, p0=[max(fitR_y), fitR_x[fitR_y.index(max(fitR_y))], 0.5], bounds=[[max(fitR_y)-50, fitR_x[fitR_y.index(max(fitR_y))]-2, 0],[max(fitR_y)+50, fitR_x[fitR_y.index(max(fitR_y))]+2, 1.5*abs(stat.pstdev(fitR_x))]])
                except:
                    pass
            # Repeat the fit for the left charge spectrum
            if len(fitL_x) != 0 and len(fitL_y) != 0:
                try:
                    fitL_p, fitL_co = curve_fit(gauss, fitL_x, fitL_y, p0=[max(fitL_y), fitL_x[fitL_y.index(max(fitL_y))], 0.5], bounds=[[max(fitL_y)-50, fitL_x[fitL_y.index(max(fitL_y))]-2, 0],[max(fitL_y)+50, fitL_x[fitL_y.index(max(fitL_y))]+2, 1.5*abs(stat.pstdev(fitL_x))]])
                except:
                    pass
            
            if fitR_p[0] != 0:
                # Use the fitted mean and sigma to define an energy cut based on the photopeaks
                # Cut value is define to be 2.5 sigma below (to the left of) the mean
                R_cut = fitR_p[1] - 2.5*fitR_p[2]
                # Only keep events with energies equal to or above the cut values
                temp = temp[temp[:,2] >= R_cut]

            if fitL_p[0] != 0:
                L_cut = fitL_p[1] - 2.5*fitL_p[2]
                temp = temp[temp[:,1] >= L_cut]


            filled_lor_counts.append(len(temp))
            if len(temp) > 200:
                temp = np.asarray(random.sample(list(temp), 200))
            coinc_list[lr] = []
            count_filled += 1
            # Add kept events to filled lors
            filled_lors.append(20000) # 20000 arbitrary number greater than 10000 used for indicator
            filled_lors.append(int(lr // 3072)) # IDL 
            filled_lors.append(int(lr % 3072)) # IDR
            for timediff in temp[:,0]:
                filled_lors.append(int(timediff))
            print(f'successfully filled L:{lr // 3072} R:{lr % 3072}')
    
    # Datafile to keep filled LORs
    filled_lors = np.asarray(filled_lors, dtype = np.int16)
    filled_lors_file = data_storage_dir + 'filled_lors_file_' + str(lorchunk).zfill(2) + '.bin'
    with open(data_storage_dir + r'\final_lors.bin', 'ab') as f:
        f.write(filled_lors.tobytes())
    print(f'{round(time.time() - start, 3)}, {count_filled} incomplete lors appended to final lors file')

# Create incomplete lor size histograms
if incomplete_lors_sizes_histogram:
    print('Number of zero counts: ' + str(len([x for x in filled_lor_counts if x == 0])))
    print('Number of nonzero counts: ' + str(len([x for x in filled_lor_counts if x > 0])))
    #np.savetxt(data_storage_dir + '/testfilledlors.txt', filled_lor_counts, delimiter='\t')
    plt.figure(1)
    plt.hist([x for x in filled_lor_counts if x > 0], bins = np.linspace(0, 500, 200), color = 'blue', alpha = 0.4)
    plt.grid()
    plt.title('Distribution of LOR populations')
    plt.xlabel('Events')
    plt.ylabel('Number of LORs')
    plt.show()


# Delete temporary LOR storage files
if delete_temporary_lor_files:
    for lorchunk in num_lor_chunks:
        os.remove(data_storage_dir + 'filled_indicators_file_' + str(lorchunk).zfill(2) + '.bin') 
        os.remove(data_storage_dir + 'incomplete_lors_file_' + str(lorchunk).zfill(2) + '.bin') 
        os.remove(data_storage_dir + 'filled_lors_file_' + str(lorchunk).zfill(2) + '.bin') 
