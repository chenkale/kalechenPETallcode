# Time calibration for MiniPET

# Center source alibration data and miniPET geometry directories
center_source_data_dir = 'C:/Users/burri/Downloads/MiniPET_Test_100mm_15min_coinc.dat'
miniPET_geo_map_dir = 'C:/Users/burri/Downloads/testMiniPETmap.csv' # Only used if generate_coincidence_heatmaps = True
output_dir = 'C:/Users/burri/Downloads/calibrationtest'

# Use gaussian fit to find mean for time offset, otherwise use simple weighted average around max
# Weighted average has minimal loss in accuracy (< 5%) but improves computation time drastically
use_guassian_fit = True
# Sample size threshold, will only perform fit or weighted average for each table entry if sample size is above this value
sample_size_threshold = 0
# Bin range above and below peak for weighted average
binrange = 5 # Only used if use_gaussian_fit = False

# Plots to be generated
# Generate before and after calibration time difference histograms with basic gaussian fit
generate_time_difference_histograms = True
# Generate before and after coincidence heatmaps
generate_coincidence_heatmaps = True

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Read calibration data into pd DataFrame
df = pd.DataFrame()
for chunk in pd.read_csv(center_source_data_dir, usecols=[2, 3, 4, 7, 8, 9], delimiter="\t", header=None, chunksize=10000):
    # Name columns
    chunk.columns = ['TimeL', 'EnergyL', 'IDL', 'TimeR', 'EnergyR', 'IDR']
    # Convert Left channel to mod 128
    chunk['IDL'] = chunk['IDL'] % 128
    # Time difference column
    chunk['TimeDiff'] = chunk['TimeR'] - chunk['TimeL']
    
    df = pd.concat([df, chunk], ignore_index=True)

# Convert to numpy array, faster
arr = df.to_numpy()

# Matrix for storing arrays of coincidences
# each row corresponds to a left channel, each col to a right channel, 
# so then each entry is a unique LOR  
coincidences_matrix = [[[] for i in range(128)] for j in range(128)]
# Creates coincidence matrix containing TimeDiffs for every LOR 
for x in range(len(arr)):
    coincidences_matrix[int(arr[x][2])][int(arr[x][5])].append(arr[x][6])

# Table that stores final time offsets
time_offset_table = np.asarray([np.asarray([0.0 for x in range(128)]) for y in range(128)])

# Gaussian
if use_guassian_fit:

    # Define gaussian fit
    def gaussian(x, a, mu, c):
        return a * np.exp(-(x - mu)**2 / (2 * c**2))

    # Loop through coincidences matrix and calculate a mean time offset for each list of coincidences
    for l in range(128):
        for r in range(128):
            # If sufficient sample size, make gaussian to find mean for adjustment
            if len(coincidences_matrix[l][r]) > sample_size_threshold:
                # Create histogram
                counts, bins = np.histogram(coincidences_matrix[l][r], bins = 100)
                # Initial guess for amplitude
                inita = 0.8 * np.max(counts)
                # Initial guess for mean
                initmu = bins[:-1][np.where(counts == np.max(counts))[0][0]]
                
                try:
                    # try fit, with bad stats the fit sometimes doesn't work
                    popt = curve_fit(gaussian, bins[:-1], counts, p0 = [inita, initmu, 250])[0]
                    # Check that mean is reasonable 
                    if abs(popt[1]) < 10000: # 10000
                        # write time offset to table as negative so that it can be ADDED when used
                        time_offset_table[l][r] = -1.0 * popt[1]
                except:
                    # otherwise leave as default 0
                    pass

# Weighted averages
else:
    for l in range(128):
        for r in range(128):
            if len(coincidences_matrix[l][r]) > sample_size_threshold:
                # Here no fit is performed to save time, minimal loss in accuracy
                counts, bins = np.histogram(coincidences_matrix[l][r], bins = 100)
                maxindex = np.where(counts == np.max(counts))[0][0]
                try:
                    time_offset_table[l][r] = -1.0 * np.average(bins[maxindex-binrange:maxindex+binrange+1], weights=counts[maxindex-binrange:maxindex+binrange+1])
                except:
                    pass

# Write table to csv
np.savetxt(output_dir + '/timeoffsettabletest.csv', time_offset_table, delimiter=',', fmt = '%.2f')

# Generate time difference histograms
if generate_time_difference_histograms:

    # Create master time difference array
    master_time_diffs = []
    for x in arr:
        master_time_diffs.append(x[6])

    # Create time difference array after calibration
    master_time_diffs_calibrated = []
    for l in range(128):
        for r in range(128):
            for timediff in coincidences_matrix[l][r]:
                master_time_diffs_calibrated.append(timediff + time_offset_table[l][r])
    
    def gaussian(x, a, mu, c):
        return a * np.exp(-(x - mu)**2 / (2 * c**2))

    def formatlabel(a, mu, c):
        return f'a={round(a)}, mu={round(mu,1)}, sig={round(abs(c),1)}'

    # Histogram before time calibration
    hist_before = plt.hist(master_time_diffs, bins = 1000, color = 'red', alpha = 0.4)
    bincenters0 = [(hist_before[1][x]+hist_before[1][x+1])/2 for x in range(len(hist_before[0]))]
    popt0 = curve_fit(gaussian, bincenters0, hist_before[0], p0 = [6000, 0, 1000])[0]
    #print(*popt0)
    plt.plot(np.linspace(-3000, 3000, 300), gaussian(np.linspace(-3000, 3000, 300), *popt0), color = 'red', label = formatlabel(*popt0))
    
    # Histogram after time calibration
    hist_after = plt.hist(master_time_diffs_calibrated, bins = 1000, color = 'blue', alpha = 0.4)
    bincenters1 = [(hist_after[1][x]+hist_after[1][x+1])/2 for x in range(len(hist_after[0]))]
    popt1 = curve_fit(gaussian, bincenters1, hist_after[0], p0 = [6000, 0, 1000])[0]
    #print(*popt1)
    plt.plot(np.linspace(-3000, 3000, 300), gaussian(np.linspace(-3000, 3000, 300), *popt1), color = 'blue', label = formatlabel(*popt1))

    plt.legend()
    plt.grid()
    plt.xlim(-3000, 3000)
    plt.title('Time differences calibrated')
    plt.xlabel('Time difference (ps)')
    plt.ylabel('Counts')
    plt.savefig(output_dir + '/time_diffence_histograms.png')
    plt.close()

# Generate heatmaps before and after for visual check
if generate_coincidence_heatmaps:
    
    # Read map, and convert to nparray to save time when referencing map
    mapdf = pd.read_csv(miniPET_geo_map_dir, header=None)
    map = []
    for row in mapdf.iloc:
        map.append(np.asarray(row[1:]))
    map = np.asarray(map)

    # Function for getting coincidence location
    def getcoinclocation(idl, idr, timediff):
        xL, yL, zL = map[idl]
        xR, yR, zR = map[idr + 128]
        #print(xL, yL, zL)
        #print(xR, yR, zR)
        # Parallel vector to LOR with length equal to light-picosecond in mm
        length = ((xR-xL)**2 + (yR-yL)**2 + (zR-zL)**2)**0.5
        xv = 0.3 * (xR-xL) / length
        yv = 0.3 * (yR-yL) / length
        zv = 0.3 * (zR-zL) / length
        return (xL+xR)/2 + xv*timediff, (yL+yR)/2 + yv*timediff, (zL+zR)/2 + zv*timediff

    coincxs, coincys, coinczs = [], [], []
    coincxs1, coincys1, coinczs1 = [], [], []
    for l in range(128):
        for r in range(128):
            for timediff in coincidences_matrix[l][r]:
                temp = getcoinclocation(l, r, timediff)
                coincxs.append(temp[0])
                coincys.append(temp[1])
                coinczs.append(temp[2])
                temp = getcoinclocation(l, r, timediff + time_offset_table[l][r])
                coincxs1.append(temp[0])
                coincys1.append(temp[1])
                coinczs1.append(temp[2])

    
    plt.hist2d(coincxs, coincys, range=[[-210, 210], [-30, 30]], bins=400, cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(output_dir + '/xybefore.png')
    plt.close()
    plt.hist2d(coincxs, coinczs, range=[[-210, 210], [-15, 15]], bins=400, cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.savefig(output_dir + '/xzbefore.png')
    plt.close()
    plt.hist2d(coincys, coinczs, range=[[-20, 20], [-20, 20]], bins=400, cmap='viridis')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.savefig(output_dir + '/yzbefore.png')
    plt.hist2d(coincxs1, coincys1, range=[[-210, 210], [-30, 30]], bins=400, cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(output_dir + '/xyafter.png')
    plt.close()
    plt.hist2d(coincxs1, coinczs1, range=[[-210, 210], [-15, 15]], bins=400, cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.savefig(output_dir + '/xzafter.png')
    plt.close()
    plt.hist2d(coincys1, coinczs1, range=[[-20, 20], [-20, 20]], bins=400, cmap='viridis')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.savefig(output_dir + '/yzafter.png')
    plt.close()
