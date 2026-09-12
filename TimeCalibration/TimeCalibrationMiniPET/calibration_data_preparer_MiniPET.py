# Calibration data preparer

'''Need to specify:

• calibration_data_file     - coinc.dat file from PETsys, no processing done outside of the DAQ software
• filled_lors_file          - file containing time differences for filled LORs, that is, LORs for which there are enough events and are ready for time calibration, these are no longer changed, only appended to
• incomplete_lors_file      - file containing time differences and energies for events in incomplete LORs, will be added to until sufficient events post cuts, then events saved to filled lors
• filled_indicators_file    - file containing filled indicators, a list of booleans indicating which lors are filled or incomplete, saves time by skipping already filled lors
• info_out_file             - file for output information to be written, including number of events in incomplete lors or 'FILLED' for filled lors

'''

calibration_data_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationMiniPET/VerticalMotion_HighIntSource_800sec_2_coinc.dat'
filled_lors_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationMiniPET/a/filled_lors.dat'
incomplete_lors_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationMiniPET/a/incomplete_lors.dat'
filled_indicators_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationMiniPET/a/filled_indicators.dat'
info_out_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationMiniPET/a/info_out_file.csv'


import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
import statistics as stat
import random
import time

start = time.time()
# Read which LORs are incomplete
if os.path.exists(filled_indicators_file):
    filled_indicators = np.fromfile(filled_indicators_file, dtype = np.bool_)
else:
    filled_indicators = np.asarray([False for lr in range(128 * 128)])
print(f'{time.time() - start} {len(filled_indicators)}')

# Group coincidences by LOR, to do this, read in incomplete LORs
coinc_list = [[] for lr in range(128 * 128)] # idL = // 128, idR = %128
if os.path.exists(incomplete_lors_file):
    num_rows = os.path.getsize(incomplete_lors_file) // 6
    incomplete_lors = np.memmap(incomplete_lors_file, dtype=np.float16, mode='r', shape=(num_rows,3))
    temp = []
    lr = 0
    for i in range(0, num_rows, 100000):
        chunk = np.asarray(incomplete_lors[i:i+100000])
        for row in chunk:
            if row[0] > 0 and row[0] < 1:
                if len(temp) != 0:
                    coinc_list[lr] = temp
                    temp = []
                lr = int(row[1] * 128 + row[2])
            else:
                temp.append(list(row))
    del incomplete_lors
print(f'{time.time() - start} read incomplete lors')

# Read calibration data
data = pd.DataFrame()
for chunk in pd.read_csv(calibration_data_file, sep = '\t', header = None, dtype = np.float64, usecols = [2,3,4,7,8,9], chunksize=1000000):
    data = pd.concat([data, chunk], ignore_index = True)
data = np.asarray(data)
print(f'{time.time() - start} calibration data read')

# Append new events to incomplete LORs 
for row in data:
    if not filled_indicators[int((row[2] % 128) * 128 + row[5])]:
        coinc_list[int((row[2] % 128) * 128 + row[5])].append([row[0] - row[3], row[1], row[4]])
print(f'{time.time() - start} coinc list filled')
del data

# Simple gaussian
def gauss(x,A,mu,sigma):
    # Simple Gaussian function to be used by the fitter
    y = A*np.exp(-(x - mu)**2 / (2 * sigma**2))
    return y


# Code to cut LORs with enough events and keep time differences of relevant (photopeak+timediff cut) events
# These cuts are taken from the ChannelPairBuilder.py and ChannelPairHeader.py scripts
filled_lors = []
count_loops = 0
count_filled = 0
for lr in range(128 * 128):
    if not filled_indicators[lr]:
        count_loops += 1
        # Perform cut only if number of events > 500
        if len(coinc_list[lr]) >= 500:
            temp = np.asarray(coinc_list[lr])

            '''FOR FINDING PHOTOPEAK'''
            
            valuesL, binsL = np.histogram(temp[:,1], bins = np.linspace(0,45,160)) # valuesL are to be fitted
            valuesR, binsR = np.histogram(temp[:,2], bins = np.linspace(0,45,160)) # valuesR are to be fitted
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
            
            '''FIT PHOTOPEAK'''

            fitR_p = [0,1,1] # Initialized to amplitude=0, mean=1, sigma=1
            fitL_p = [0,1,1]
            # Set the range of x and y values to fit: from 7 bins to the left of the photopeak, up to the end of the spectra
            fitR_x = list(bin_centers[PhotopeakR-7:])
            fitL_x = list(bin_centers[PhotopeakL-7:])
            fitR_y = list(valuesR[PhotopeakR-7:])
            fitL_y = list(valuesL[PhotopeakL-7:])

            # Make sure fit ranges are non-empty, perform the fits, keep parameters
            if len(fitR_y) != 0 and len(fitR_x) != 0:
                try:
                    # curve_fit(function, x_data, y_data, p0=param_guesses, bounds=lower and upper param bounds)
                    fitR_p, fitR_co = curve_fit(gauss, fitR_x, fitR_y, p0=[max(fitR_y), fitR_x[fitR_y.index(max(fitR_y))], 0.5], bounds=[[max(fitR_y)-50, fitR_x[fitR_y.index(max(fitR_y))]-2, 0],[max(fitR_y)+50, fitR_x[fitR_y.index(max(fitR_y))]+2, 1.5*abs(stat.pstdev(fitR_x))]])
                except:
                    #print('Left photpeak fit failed')
                    continue
            # Repeat the fit for the left charge spectrum
            if len(fitL_x) != 0 and len(fitL_y) != 0:
                try:
                    fitL_p, fitL_co = curve_fit(gauss, fitL_x, fitL_y, p0=[max(fitL_y), fitL_x[fitL_y.index(max(fitL_y))], 0.5], bounds=[[max(fitL_y)-50, fitL_x[fitL_y.index(max(fitL_y))]-2, 0],[max(fitL_y)+50, fitL_x[fitL_y.index(max(fitL_y))]+2, 1.5*abs(stat.pstdev(fitL_x))]])
                except:
                    #print('Right photopeak fit failed')
                    continue

            '''CUT ON PHOTOPEAK'''

            # Use the fitted mean and sigma to define an energy cut based on the photopeaks
            # Cut value is define to be 2.5 sigma below (to the left of) the mean
            R_cut = fitR_p[1] - 2.5*fitR_p[2]
            L_cut = fitL_p[1] - 2.5*fitL_p[2]
            # Only keep events with energies equal to or above the cut values
            temp = temp[temp[:,1] >= L_cut]
            temp = temp[temp[:,2] >= R_cut]

            '''FIT TIME DIFF'''

            # Histogram the time differences for photopeak events and fit to Gaussian
            fit_p = [0,1,1] # time difference fit parameters
            # Set ranges for the data and histograms based on the max and min values
            try:
                low = temp[:,0].min()
                high = temp[:,0].max()
            except:
                continue
            # Calculate the bin numbers and sizes so they are constant to help with comparisons
            num_bins = int((high-low)/50)
            bins = np.linspace(low,high,num_bins)
            #sometimes the time differences are real close together, ya see
            if len(bins) < 20:
                bins = np.linspace(low, high, 20)
            # Create the time difference histogram
            values, bins = np.histogram(temp[:,0], bins=bins)
            # Prepare arguments for the fitter
            peak = np.argmax(values) # Estimate the peak to be the bin with the most counts
            bin_centers = bins + ((bins.max() - bins.min()) / (2*num_bins))
            bin_centers.resize((len(bin_centers) -1))
            # Set the x and y range for the fitter to be 10 bins around the peak
            fit_x = bin_centers[max(0,peak-10):peak+11]
            fit_y = values[max(0,peak-10):peak+11]
            # Make sure fit ranges are non-empty, perform the fits, keep parameters
            if len(fit_x) != 0 and len(fit_y) != 0:
                try:
                    fit_p, fit_co = curve_fit(gauss, fit_x, fit_y, p0=[max(fit_y), stat.mean(fit_x), 0.5*stat.pstdev(fit_x)], bounds=[[max(fit_y)-10,stat.mean(fit_x)-100,0],[max(fit_y)+10,stat.mean(fit_x)+100,1.5*abs(stat.pstdev(fit_x))]])
                except:
                    #print("Error - CTR curve fit failed!!!")
                    continue
            
            '''CUT ON TIME DIFF'''

            time_cut_min = fit_p[1] - 5 * fit_p[2]
            time_cut_max = fit_p[1] + 5 * fit_p[2]
            temp = temp[temp[:,0] >= time_cut_min]
            temp = temp[temp[:,0] <= time_cut_max]

            '''CHECK ENOUGH EVENTS'''
            
            if len(temp) >= 100:
                if len(temp) > 200:
                    temp = np.asarray(random.sample(list(temp), 200))
                coinc_list[lr] = []
                filled_indicators[lr] = 1
                count_filled += 1
                # Add kept events to filled lors
                filled_lors.append(0.1)
                filled_lors.append(lr // 128)
                filled_lors.append(lr % 128)
                for timediff in temp[:,0]:
                    filled_lors.append(timediff)
                print(f'successfully filled L:{lr // 128} R:{lr % 128}')
                

print(f'{time.time() - start}, {count_loops} loops, {count_filled} filled')

# Datafile to keep filled LORs
filled_lors = np.asarray(filled_lors, dtype = np.float16)
if os.path.exists(filled_lors_file):
    with open(filled_lors_file, 'ab') as f:
        filled_lors.tofile(f)
else:
    filled_lors.tofile(filled_lors_file)

# Datafile to keep incomplete LORs
out_data = []
for lr in range(128 * 128):
    if not filled_indicators[lr]:
        out_data.append([0.1, lr // 128, lr % 128])
        for row in coinc_list[lr]:
            out_data.append(row)
out_data = np.asarray(out_data, dtype = np.float16)
out_data.tofile(incomplete_lors_file)

# Datafile to keep track of filled indicators
filled_indicators.tofile(filled_indicators_file)

# Export information 
out_info = [['IDL', 'IDR', 'Events']]
for lr in range(128 * 128):
    if filled_indicators[lr]:
        out_info.append([lr // 128, lr % 128 + 896, 'FILLED'])
    else:
        out_info.append([lr // 128, lr % 128 + 896, len(coinc_list[lr])])
out_df = pd.DataFrame(out_info)
out_df.to_csv(info_out_file, sep = ',', index = False, header = False)