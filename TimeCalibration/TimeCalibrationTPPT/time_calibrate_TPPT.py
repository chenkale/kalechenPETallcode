# Python script for generating time_offset_calibration.tsv for TPPT 

'''
ONLY NEED TO CHANGE: center_source_data_dir_BOTTOM, center_source_data_dir_TOP, TPPT_geo_map_dir, output_dir
'''

# Center source calibration data and TPPT geometry directories
# If only one file, use BOTTOM and leave TOP as ''
center_source_data_dir_BOTTOM = 'C:/Users/burri/Downloads/MovingSource_LowInt_60s_Bottom_HWTrigOff_coinc.dat'
center_source_data_dir_TOP = 'C:/Users/burri/Downloads/MovingSource_LowInt_60s_Top_HWTrigOff_coinc.dat'
TPPT_geo_map_dir = 'C:/Users/burri/Downloads/TPPT_Scanner_mapcopy.csv'
output_dir = 'C:/Users/burri/Downloads/TPPTtimecalibrationtest/'
filled_lors_file = ''

# Generate time difference histograms
generate_time_difference_histograms = True
# Generate coincidence heatmaps for visual check
generate_coincidence_heatmaps = False

num_lor_chunks = 8

# Fitted bins for lookup values
fitted_bins = 200
# Bin range for cuts
bin_range = 50

# Histogram limit for time difference histograms
histogramlim = 10000
# Bin width for histogram
binwidth = 50


import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import tqdm

# Read calibration data into pd DataFrame
df = pd.DataFrame()
for chunk in pd.read_csv(center_source_data_dir_BOTTOM, usecols=[2, 4, 7, 9], delimiter="\t", header=None, chunksize=100000):
    # Name columns
    chunk.columns = ['TimeL', 'IDL', 'TimeR', 'IDR']
    # Time difference column, left minus right lets positive time differences be on the right side
    chunk['TimeDiff'] = chunk['TimeL'] - chunk['TimeR']
    
    df = pd.concat([df, chunk], ignore_index=True)
if center_source_data_dir_TOP != '':
    for chunk in pd.read_csv(center_source_data_dir_TOP, usecols=[2, 4, 7, 9], delimiter="\t", header=None, chunksize=100000):
        # Name columns
        chunk.columns = ['TimeL', 'IDL', 'TimeR', 'IDR']
        # Time difference column, left minus right lets positive time differences be on the right side
        chunk['TimeDiff'] = chunk['TimeL'] - chunk['TimeR']
    
        df = pd.concat([df, chunk], ignore_index=True)
# Convert to numpy array, faster
arr = df.to_numpy()

num_rows = os.path.getsize(filled_lors_file) // 2
filled_lors = np.memmap(filled_lors_file, dtype=np.int16, mode='r', shape=(num_rows,))
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

# Geometric time offset matrix, stores time difference corresponding to the geometric center for each LOR
# each row corresponds to a left channel, each col to a right channel, such that each entry is a unique LOR
geometric_time_offsets = [[0 for i in range(3072)] for j in range(3072)]

# Read map
mapdf = pd.read_csv(TPPT_geo_map_dir, usecols=[0,1,2], header=None)
map = mapdf.to_numpy()

# Get geometric time offset for idl and idr channels to be added to time difference
def geometric_time_offset(idl, idr):
    xL, yL, zL = map[idl + 3072]
    xR, yR, zR = map[idr]
    # Get midpoint between two channels
    mX, mY, mZ = (xL+xR)/2, (yL+yR)/2, (zL+zR)/2

    # If midpoint lies on x = 0, return 0
    if mX == 0:
        return 0
    # else need to calculate distance between x = 0 and midpoint along LOR
    else:
        # Parameterize, below are slopes with x being the parameter
        ay, az = (yR-yL)/(xR-xL), (zR-zL)/(xR-xL)
        # Find y and z for positive distance xL away from xL (y0,z0 at x = 0)
        y0, z0 = yL + ay * -1 * xL, zL + az * -1 * xL
        # Distance, converted to light picoseconds
        d = ((mX)**2 + (y0-mY)**2 + (z0-mZ)**2)**0.5 / 0.299792458 #mm/ps
        # Find direction (midpoint left or right of x=0)
        if -1 * xL > xR:
            return -1 * d
        else:
            return d

for idl in tqdm(3072):



# Matrix for storing ARRAYS of coincidences, adjusted for geometric offset
# each row corresponds to a left channel, each col to a right channel, 
# so then each entry is a unique LOR 
coincidences_matrix = [[[] for i in range(3072)] for j in range(3072)]

# Creates coincidence matrix containing TimeDiffs (relative to x=0) for every LOR 
for row in arr:
    # Convert to 0-3072, mod3072 works because 0=0*1024+3072n, 4096=1*1024+3072n, 8196=2*1024+3072n
    idl = (int(row[1]) - 131072) % 3072
    idr = int(row[3]) % 3072 
    # Create coincidence matrix with geometric time offsets
    coincidences_matrix[idl][idr].append(row[4] + geometric_time_offsets[idl][idr])

# Array for storing left and right module channel time offsets, these will be used to calculate per channel offset values
coincidences_arr_l = [[] for i in range(3072)]
coincidences_arr_r = [[] for i in range(3072)]

# Fill arrays
for l in range(3072):
    for r in range(3072):
        coincidences_arr_l[l].extend(coincidences_matrix[l][r])
        coincidences_arr_r[r].extend(coincidences_matrix[l][r])

# Create lookup table for each channel in l and r
lookup_l = np.asarray([0.0 for i in range(3072)])
lookup_r = np.asarray([0.0 for i in range(3072)])

# Define gaussian fit
def gaussian(x, a, mu, c):
    return a * np.exp(-(x - mu)**2 / (2 * c**2))

# Calibration function
def calibrate():
    # Left side
    for l in range(3072):
        # Create histogram
        counts, bins = np.histogram(coincidences_arr_l[l], bins = fitted_bins)
        # Initial guess for amplitude
        inita = 0.8 * np.max(counts)
        # Initial guess for mean
        maxindex = np.where(counts == np.max(counts))[0][0]
        initmu = bins[:-1][maxindex]
        # Fit range
        bl, br = maxindex - bin_range, maxindex + bin_range
        try:
            # try fit, with bad stats the fit sometimes doesn't work
            popt = curve_fit(gaussian, bins[:-1][bl:br], counts[bl:br], p0 = [inita, initmu, 500])[0]
            # Check that mean is reasonable 
            if abs(popt[1]) < 10000: # 10000
                # Write time offset to lookup table, value is subtracted from event time before groups processing
                lookup_l[l] = popt[1] / 2 + lookup_l[l]
        except:
            pass

    # Right side
    for r in range(3072):
        # Create histogram
        counts, bins = np.histogram(coincidences_arr_r[r], bins = fitted_bins)
        # Initial guess for amplitude
        inita = 0.8 * np.max(counts)
        # Initial guess for mean
        maxindex = np.where(counts == np.max(counts))[0][0]
        initmu = bins[:-1][maxindex]
        # Fit range
        bl, br = maxindex - bin_range, maxindex + bin_range           
        try:
            # try fit, with bad stats the fit sometimes doesn't work
            popt = curve_fit(gaussian, bins[:-1][bl:br], counts[bl:br], p0 = [inita, initmu, 500])[0]
            # Check that mean is reasonable 
            if abs(popt[1]) < 10000: # 10000
                # Write time offset to lookup table, value is subtracted from event time before groups processing
                lookup_r[r] = -1 * popt[1] / 2 + lookup_r[r]
        except:
            pass

# Updates coincidence matrix and left/right arrays
def update_coincidences():
    # Update coincidence matrix
    for l in range(3072):
        for r in range(3072):
            temp = []
            for t in coincidences_matrix[l][r]:
                temp.append(t - lookup_l[l] + lookup_r[r])
            coincidences_matrix[l][r] = temp

    # Update left/right arrays
    global coincidences_arr_l
    coincidences_arr_l = [[] for i in range(3072)]
    global coincidences_arr_r
    coincidences_arr_r = [[] for i in range(3072)]

    # Fill arrays
    for l in range(3072):
        for r in range(3072):
            coincidences_arr_l[l].extend(coincidences_matrix[l][r])
            coincidences_arr_r[r].extend(coincidences_matrix[l][r])

# We can try calibrating multiple times to tighten distribution, twice works best so far
calibrate()
update_coincidences()
calibrate()
#update_coincidences()
#calibrate()

# Prepare table to write to output file
output = pd.DataFrame()

for r in range(3072):
    output = pd.concat([output, pd.DataFrame([[0, r // 1024, (r % 1024) // 64, r % 64, lookup_r[r]]])], ignore_index=True)
for l in range(3072):
    output = pd.concat([output, pd.DataFrame([[1, l // 1024, (l % 1024) // 64, l % 64, lookup_l[l]]])], ignore_index=True)

# Write to tsv
output.to_csv(output_dir + 'time_offset_calibration.tsv', sep='\t', header=None, float_format = '%.7e', index=False)

# Generate time difference histograms
if generate_time_difference_histograms:

    # Create master time difference array
    master_time_diffs = []
    for x in arr:
        master_time_diffs.append(x[4] + geometric_time_offsets[(int(x[1])-131072)%3072][int(x[3])%3072])
    
    # Calibrated time difference array
    master_time_diffs_calibrated = []
    for x in arr:
        master_time_diffs_calibrated.append(x[0] - lookup_l[(int(x[1])-131072)%3072] - x[2] + lookup_r[int(x[3])%3072] + geometric_time_offsets[(int(x[1])-131072)%3072][int(x[3])%3072])


    def formatlabel(a, mu, c):
        return f'a={round(a)}, mu={round(mu,1)}, sig={round(abs(c),1)}'

    # Define double gaussian
    def doublegaussian(x, a1, mu1, c1, a2, mu2, c2, noise):
        return a1 * np.exp(-(x - mu1)**2 / (2 * c1**2)) + a2 * np.exp(-(x - mu2)**2 / (2 * c2**2)) + noise

    
    # Histogram before time calibration
    hist_before = plt.hist(master_time_diffs, bins=range(int(min(master_time_diffs)), int(max(master_time_diffs)) + binwidth, binwidth), color = 'red', alpha = 0.4)
    bincenters0 = [(hist_before[1][x]+hist_before[1][x+1])/2 for x in range(len(hist_before[0]))]
    popt0 = curve_fit(doublegaussian, bincenters0, hist_before[0], p0 = [200000, 0, 1000, 10000, 0, 2500, 10000])[0]
    plt.plot(np.linspace(-1*histogramlim, histogramlim, 300), doublegaussian(np.linspace(-1*histogramlim, histogramlim, 300), *popt0), color = 'red', label = 'before')
    print(popt0)

    # Histogram after time calibration
    hist_after = plt.hist(master_time_diffs_calibrated, bins=range(int(min(master_time_diffs_calibrated)), int(max(master_time_diffs_calibrated)) + binwidth, binwidth), color = 'blue', alpha = 0.4)
    bincenters1 = [(hist_after[1][x]+hist_after[1][x+1])/2 for x in range(len(hist_after[0]))]
    popt1 = curve_fit(doublegaussian, bincenters1, hist_after[0], p0 = [250000, 0, 500, 50000, 0, 1000, 10000])[0]
    plt.plot(np.linspace(-1*histogramlim, histogramlim, 300), doublegaussian(np.linspace(-1*histogramlim, histogramlim, 300), *popt1), color = 'blue', label = 'after')
    print(popt1)

    plt.legend()
    plt.grid()
    plt.xlim(-1*histogramlim, histogramlim)
    plt.title('Time differences calibrated')
    plt.xlabel('Time difference (ps)')
    plt.ylabel('Counts')
    plt.savefig(output_dir + 'time_diffence_histograms.png')
    plt.close()

# Generate heatmaps before and after for visual check
if generate_coincidence_heatmaps:

    # Function for getting coincidence location
    def getcoinclocation(idl, idr, timediff):
        xL, yL, zL = map[idl]
        xR, yR, zR = map[idr + 3072]
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
    for l in range(3072):
        for r in range(3072):
            for timediff in coincidences_matrix[l][r]:
                temp = getcoinclocation(l, r, timediff)
                coincxs.append(temp[0])
                coincys.append(temp[1])
                coinczs.append(temp[2])
                temp = getcoinclocation(l, r, timediff - lookup_l[l] + lookup_r[r])
                coincxs1.append(temp[0])
                coincys1.append(temp[1])
                coinczs1.append(temp[2])

    if generate_coincidence_heatmaps:
        plt.hist2d(coincxs, coincys, range=[[-210, 210], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(output_dir + 'heatmap_before_xy.png')
        plt.close()
        plt.hist2d(coincxs, coinczs, range=[[-210, 210], [-15, 15]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.savefig(output_dir + 'heatmap_before_xz.png')
        plt.close()
        plt.hist2d(coincys, coinczs, range=[[-20, 20], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.savefig(output_dir + 'heatmap_before_yz.png')
        plt.hist2d(coincxs1, coincys1, range=[[-210, 210], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(output_dir + 'heatmap_after_xy.png')
        plt.close()
        plt.hist2d(coincxs1, coinczs1, range=[[-210, 210], [-15, 15]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.savefig(output_dir + 'heatmap_after_xz.png')
        plt.close()
        plt.hist2d(coincys1, coinczs1, range=[[-20, 20], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.savefig(output_dir + 'heatmap_after_yz.png')
        plt.close()