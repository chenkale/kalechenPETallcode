# Python script for generating time_offset_calibration.tsv for TPPT 

'''
ONLY NEED TO CHANGE: center_source_data_dir_BOTTOM, center_source_data_dir_TOP, TPPT_geo_map_dir, output_dir
'''

# Center source calibration data and TPPT geometry directories
# If only one file, use BOTTOM and leave TOP as ''
data_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\compressed.bin'
TPPT_geo_map_dir = 'C:/Users/burri/Downloads/TPPT_Scanner_mapcopy.csv'
output_dir = 'C:/Users/burri/Downloads/TPPTtimecalibrationtest/'

# Generate time difference histograms
generate_time_difference_histograms = False
# Generate coincidence heatmaps for visual check
generate_coincidence_heatmaps = False


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
import stat

# Memmap of data avoids reading entire file into memory
num_rows = os.path.getsize(data_dir) // 10 # 10 = 2 (size of int16) * 5 (entries per row)
data = np.memmap(data_dir, dtype=np.int16, mode='r', shape=(num_rows,5))

'''
Because the PETsys event processing only accepts per channel offset values, calibration will be done as follows:
1. For every channel, collect events for every LOR that includes said channel
2. Perform photopeak cuts per LOR and calculate mean offset (in time differences) per LOR (gaussian fit)
3. Fit offsets for each LOR per channel to calculate mean offset per CHANNEL, this will be the value written to tsv
'''

def gaussian(x, a, mu, c):
    return a * np.exp(-(x - mu)**2 / (2 * c**2))

# Function for finding photopeak, modified from ChannelPairBuilder.py
def findPhotopeak(lor_coincs):
    bins = np.linspace(0, 4500, 50) # Modification: int16s used from file_chunker.py require x100 on upper energy bound
    valuesL = np.histogram(lor_coincs[:,0], bins = bins)[0] # valuesL are to be fitted
    valuesR = np.histogram(lor_coincs[:,1], bins = bins)[0] # valuesR are to be fitted
    # Binning will be same for both so just define bin centers once w.r.t. binsR
    bin_centers = [(bins[q] + bins[q+1])/2 for q in range(0,len(bins)-1)]
                
    # Find the photopeak within 3 bins starting from the RHS of the graph
    for i in range(1, len(valuesR)-4):
        PhotopeakR = 0
        if valuesR[-1*i] > valuesR[-1*(i+1)] and valuesR[-1*i] > valuesR[-1*(i+2)] and valuesR[-1*i] > valuesR[-1*(i+3)] and valuesR[-1*(i+1)] != 0 and valuesR[-1*i] >= 4:
            PhotopeakR = len(valuesR)-1*i
            break
    for i in range(1, len(valuesL)-4):
        PhotopeakL = 0
        if valuesL[-1*i] > valuesL[-1*(i+1)] and valuesL[-1*i] > valuesL[-1*(i+2)] and valuesL[-1*i] > valuesL[-1*(i+3)] and valuesL[-1*(i+1)] != 0 and valuesL[-1*i] >= 4:
            PhotopeakL = len(valuesL)-1*i
            break
                
    fitR_p = [0,1,1] # Initialized to amplitude=0, mean=1, sigma=1
    fitL_p = [0,1,1]
    # Set the range of x and y values to fit: from 3 bins to the left of the photopeak, up to the end of the spectra
    fitR_x = list(bin_centers[PhotopeakR-3:])
    fitL_x = list(bin_centers[PhotopeakL-3:])
    fitR_y = list(valuesR[PhotopeakR-3:])
    fitL_y = list(valuesL[PhotopeakL-3:])

    # Make sure fit ranges are non-empty, perform the fits, keep parameters
    if len(fitR_y) != 0 and len(fitR_x) != 0:
        try:
            # curve_fit(function, x_data, y_data, p0=param_guesses, bounds=lower and upper param bounds)
            fitR_p, fitR_co = curve_fit(gaussian, fitR_x, fitR_y, p0=[max(fitR_y), fitR_x[fitR_y.index(max(fitR_y))], 0.5], bounds=[[max(fitR_y)-50, fitR_x[fitR_y.index(max(fitR_y))]-2, 0],[max(fitR_y)+50, fitR_x[fitR_y.index(max(fitR_y))]+2, 1.5*abs(stat.pstdev(fitR_x))]])
        except:
            #print('Left photpeak fit failed')
            pass
    # Repeat the fit for the left charge spectrum
    if len(fitL_x) != 0 and len(fitL_y) != 0:
        try:
            fitL_p, fitL_co = curve_fit(gaussian, fitL_x, fitL_y, p0=[max(fitL_y), fitL_x[fitL_y.index(max(fitL_y))], 0.5], bounds=[[max(fitL_y)-50, fitL_x[fitL_y.index(max(fitL_y))]-2, 0],[max(fitL_y)+50, fitL_x[fitL_y.index(max(fitL_y))]+2, 1.5*abs(stat.pstdev(fitL_x))]])
        except:
            #print('Right photopeak fit failed')
            pass
    return fitR_p, fitL_p

# Function for cutting on photopeak, modified from ChannelPairBuilder.py
# returns an array of only time differences (discard energy information after use)
def cutOnPhotopeak(lor_coincs, fitR_p, fitL_p):
    
    temp = np.asarray(lor_coincs)
    # Use the fitted mean and sigma to define an energy cut based on the photopeaks
    # Cut value is define to be 2.5 sigma below (to the left of) the mean
    R_cut = fitR_p[1] - 2.5*fitR_p[2]
    L_cut = fitL_p[1] - 2.5*fitL_p[2]
    # Only keep events with energies equal to or above the cut values
    temp = temp[temp[:,0] >= L_cut]
    temp = temp[temp[:,1] >= R_cut]
    
    return np.asarray(temp[:,2])

# Function for calculating mean offset for one LOR
# returns mean offset, number of passing events
def calculateMeanOffset(lor_time_diffs):
    # Create histogram
    counts, bins = np.histogram(lor_time_diffs, bins = 50)
    bin_centers = [(bins[q] + bins[q+1])/2 for q in range(0,len(bins)-1)]
    # Initial guess for amplitude
    inita = 0.8 * np.max(counts)
    # Initial guess for mean
    maxindex = np.where(counts == np.max(counts))[0][0]
    initmu = bin_centers[maxindex]
    # Fit range
    bl, br = maxindex - 10, maxindex + 10
    try:
        # try fit, with bad stats the fit sometimes doesn't work
        popt = curve_fit(gaussian, bin_centers[bl:br], counts[bl:br], p0 = [inita, initmu, 500])[0]
        # Check that mean is reasonable 
        if abs(popt[1]) < 10000: # 10000
            # Return mean offset, number of passing events
            return popt[1] #, sum(counts[bl:br])
    except:
        return None

# Function for calculating mean of mean offsets for one channel
def calculateMeanMeanOffset(mean_offsets):
    # Create histogram
    counts, bins = np.histogram(mean_offsets, bins = 50)
    bin_centers = [(bins[q] + bins[q+1])/2 for q in range(0,len(bins)-1)]
    # Initial guess for amplitude
    inita = 0.8 * np.max(counts)
    # Initial guess for mean
    maxindex = np.where(counts == np.max(counts))[0][0]
    initmu = bin_centers[maxindex]
    try:
        # try fit, with bad stats the fit sometimes doesn't work
        popt = curve_fit(gaussian, bin_centers, counts, p0 = [inita, initmu, 100])[0]
        # Check that mean is reasonable 
        if abs(popt[1]) < 10000: # 10000
            # Return mean offset, number of passing events
            return popt[1] #, sum(counts[bl:br])
    except:
        return None

# Read map
map = pd.read_csv(TPPT_geo_map_dir, usecols=[0,1,2], header=None)
map = map.to_numpy()
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

# Lookup arrays to store final values to be written to tsv
lookup_l = np.asarray([0.0 for idl in range(3072)])
lookup_r = np.asarray([0.0 for idr in range(3072)])

# Calibration function left
def calibrateleft():
    for idl in tqdm(range(3072)):

        # Create list to store coincidence events
        coinc_list = [[] for idr in range(3072)]

        # Loop through data to find coincidences for this channel
        for i in range(0, num_rows, 1000000):
            chunk = np.asarray(data[i:i+1000000]) # Read as chunks
            # Append events to coinc list
            for row in chunk:
                if row[1] == idl:
                    coinc_list[row[3]].append(np.asarray((row[0], row[2], row[4] + lookup_l[idl] + lookup_r[row[3]])))
    
        # Loop through collected data by LOR and perform photopeak cuts
        for idr in range(3072):
            # Only attempt if > 100 events
            if len(coinc_list[idr]) > 100:
                coinc_list[idr] = np.asarray(coinc_list[idr])
                #init_event_count = len(coinc_list[idr]) # Keep track of initial event count
                fitR_p, fitL_p = findPhotopeak(coinc_list[idr]) # Finds photopeak parameters
                coinc_list[idr] = cutOnPhotopeak(coinc_list[idr], fitR_p, fitL_p) # Cut on photopeak, returns array of only time differences
                # Only calculate mean if > 50 events passed
                if len(coinc_list[idr]) > 50:
                    coinc_list[idr] = calculateMeanOffset(coinc_list[idr]) # Calculate mean offset, returns only the mean
                    coinc_list[idr] = coinc_list[idr] + geometric_time_offset(idl, idr) # Apply geometric offset
                else:
                    coinc_list[idr] = None
            else: 
                coinc_list[idr] = None
    
        # Now, coinc_list should be a list of singular mean values or Nones, so we may calculate mean of mean offsets
        coinc_list = np.asarray(coinc_list)
        coinc_list = coinc_list[coinc_list != None]
        if len(coinc_list) > 100: # Only calculate if number of filled LORs > 100
            offset = calculateMeanMeanOffset(coinc_list)
            if offset == None:
                print(f'Offset for left channel {idl} could not be calculated')
            else:
                lookup_l[idl] += offset / 2
    
# Calibration function right
def calibrateright():
    for idr in tqdm(range(3072)):

        # Create list to store coincidence events
        coinc_list = [[] for idl in range(3072)]

        # Loop through data to find coincidences for this channel
        for i in range(0, num_rows, 1000000):
            chunk = np.asarray(data[i:i+1000000]) # Read as chunks
            # Append events to coinc list
            for row in chunk:
                if row[3] == idr:
                    coinc_list[row[1]].append(np.asarray((row[0], row[2], row[4] + lookup_l[row[1]] + lookup_r[idr])))
    
        # Loop through collected data by LOR and perform photopeak cuts
        for idl in range(3072):
            # Only attempt if > 100 events
            if len(coinc_list[idl]) > 100:
                coinc_list[idl] = np.asarray(coinc_list[idl])
                #init_event_count = len(coinc_list[idr]) # Keep track of initial event count
                fitR_p, fitL_p = findPhotopeak(coinc_list[idl]) # Finds photopeak parameters
                coinc_list[idl] = cutOnPhotopeak(coinc_list[idl], fitR_p, fitL_p) # Cut on photopeak, returns array of only time differences
                # Only calculate mean if > 50 events passed
                if len(coinc_list[idl]) > 50:
                    coinc_list[idl] = calculateMeanOffset(coinc_list[idl]) # Calculate mean offset, returns only the mean
                    coinc_list[idl] = coinc_list[idl] + geometric_time_offset(idl, idr) # Apply geometric offset
                else:
                    coinc_list[idl] = None
            else: 
                coinc_list[idl] = None
    
        # Now, coinc_list should be a list of singular mean values or Nones, so we may calculate mean of mean offsets
        coinc_list = np.asarray(coinc_list)
        coinc_list = coinc_list[coinc_list != None]
        if len(coinc_list) > 100: # Only calculate if number of filled LORs > 100
            offset = calculateMeanMeanOffset(coinc_list)
            if offset == None:
                print(f'Offset for left channel {idr} could not be calculated')
            else:
                lookup_r[idr] += -1 * offset / 2


# Prepare table to write to output tsv
output = pd.DataFrame()

for r in range(3072):
    output = pd.concat([output, pd.DataFrame([[0, r // 1024, (r % 1024) // 64, r % 64, lookup_r[r]]])], ignore_index=True)
for l in range(3072):
    output = pd.concat([output, pd.DataFrame([[1, l // 1024, (l % 1024) // 64, l % 64, lookup_l[l]]])], ignore_index=True)

# Write to tsv
output.to_csv(r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\tsvtest1.tsv', sep='\t', header=None, float_format = '%.7e', index=False)

del data
