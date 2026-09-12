# Calibration data preparer for TPPT
# Compresses .dat files to be used for time calibration, cuts on photopeak, keeps only necessary information

'''Need to specify:

• calibration_data_file     - coinc.dat file after chunking
• data_storage_dir          - location to store temporary files

'''

#calibration_data_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationTPPT/MovingSource_LowInt_60s_Bottom_HWTrigOff_coinc.dat'
#calibration_data_file = 'C:/Users/burri/Documents/PET/TimeCalibration/TimeCalibrationMiniPET/VerticalMotion_HighIntSource_800sec_coinc.dat'
calibration_data_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\run1'
data_storage_dir = 'C:/Users/burri/Documents/PET/TimeCalibration/CalibrationDataPreparer/storage/'

# Process data in chunks of LORs
num_lor_chunks = 8
lorchunksize = int(3072 * 3072 / num_lor_chunks)

import numpy as np
import os
from scipy.optimize import curve_fit
import statistics as stat
import random
import time
import matplotlib.pyplot as plt
from tqdm import tqdm

start = time.time() # For keeping track of time

# Simple gaussian
def gauss(x,A,mu,sigma):
    return A*np.exp(-(x - mu)**2 / (2 * sigma**2))

# Convert between chunk lr and true lr:
# When reading in chunks, use a chunk lr index for referencing in arrays, etc.
# but need to convert this back to true lr when writing to files
def tolrc(lrt): # To chunk lr
    return int(lrt % lorchunksize)
def tolrt(lrc, lorchunk): # To true lr
    return int(lrc + lorchunk * lorchunksize)

#out_info = [['IDL', 'IDR', 'Events']]

filled_lor_counts = []
for lorchunk in tqdm(range(num_lor_chunks)):
    
    # Read filled_indicators boolean file that keeps track of which LORs are filled
    filled_indicators_file = os.path.join(data_storage_dir, f'filled_indicators_file_{str(lorchunk).zfill(2)}.bin')
    if os.path.exists(filled_indicators_file): 
        filled_indicators = np.fromfile(filled_indicators_file, dtype = np.bool_) # Read from file if exists
    else:
        filled_indicators = np.asarray([False for lrc in range(lorchunksize)]) # Otherwise start anew
    print(f'{round(time.time() - start, 3)} read filled indicators {len(filled_indicators)}')

    # Read in incomplete LORs, skip already filled LORs to save computation and storage space
    incomplete_lors_file = os.path.join(data_storage_dir, f'incomplete_lors_file_{str(lorchunk).zfill(2)}.bin')
    coinc_list = [[] for lrc in range(lorchunksize)] # idL = lr // 3072, idR = lr % 3072
    if os.path.exists(incomplete_lors_file): 
        count = 0
        num_rows = os.path.getsize(incomplete_lors_file) // 12 # 12 = 4 (size of float32) * 3 (entries per row: timediff, energyL, energyR)
        incomplete_lors = np.memmap(incomplete_lors_file, dtype=np.float32, mode='r', shape=(num_rows,3))
        lrc = 0 # For matching file read to proper LOR
        for i in range(0, num_rows, 1000000):
            chunk = np.asarray(incomplete_lors[i:i+1000000]) # Read as chunks
            for row in chunk:
                if row[0] > 0 and row[0] < 1: # Indicator for new LOR is a non integer
                    lrc = tolrc(row[1] * 3072 + row[2]) # Match proper LOR
                else:
                    coinc_list[lrc].append(list(row))
                    count += 1
        del incomplete_lors, chunk, row # Delete memmap and dependent objects after use
        print(f'{round(time.time() - start, 3)} read incomplete lors {count}')
    else: 
        print(f'{round(time.time() - start, 3)} read incomplete lors 0')
    
    # Read chunked data, np.memmap allows flexibility in RAM usage
    idrange = range(lorchunk * lorchunksize, (lorchunk + 1) * lorchunksize) # Define range of compositeids to be handled on this iteration
    for chunknum in range(len(os.listdir(calibration_data_dir))):
        num_rows = os.path.getsize(os.path.join(calibration_data_dir, f'chunk{chunknum}.bin')) // 20 # num_rows needed for shape
        datammap = np.memmap(os.path.join(calibration_data_dir, f'chunk{chunknum}.bin'), dtype=np.float32, mode='r', shape=(num_rows,5))
        for i in range(0, num_rows, 1000000):
            chunk = np.asarray(datammap[i:i+1000000]) # Read as chunks
            # Append new events from calibration data to incomplete LORs
            for row in chunk:
                idcomp = int(row[1] * 3072 + row[3]) # Composite ID
                if idcomp in idrange:
                    if not filled_indicators[idcomp - lorchunk * lorchunksize]:
                        coinc_list[idcomp - lorchunk * lorchunksize].append([row[4], row[0], row[2]]) # time difference, energyL, energyR 
        del datammap, chunk, row
        print(f'{round(time.time() - start, 3)} read cal data chunk {chunknum}')   
    print(f'{round(time.time() - start, 3)} read calibration data')    

    
    # Code to cut LORs with enough events and keep time differences of relevant (photopeak cut) events
    # These cuts are taken from the ChannelPairBuilder.py and ChannelPairHeader.py scripts
    filled_lors = []
    count_loops = 0
    count_filled = 0
    for lrc in range(lorchunksize):
        if not filled_indicators[lrc]:
            count_loops += 1
            # Perform cut only if number of events > 100
            if len(coinc_list[lrc]) >= 100:
                temp = np.asarray(coinc_list[lrc])

                '''FOR FINDING PHOTOPEAK'''
                
                valuesL, binsL = np.histogram(temp[:,1], bins = np.linspace(0,4500,50)) # valuesL are to be fitted
                valuesR, binsR = np.histogram(temp[:,2], bins = np.linspace(0,4500,50)) # valuesR are to be fitted
                # Binning will be same for both so just define bin centers once w.r.t. binsR
                bin_centers = [(binsR[q] + binsR[q+1])/2 for q in range(0,len(binsR)-1)]
                
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
                
                '''FIT PHOTOPEAK'''

                fitR_p = [0,1,1] # Initialized to amplitude=0, mean=1, sigma=1
                fitL_p = [0,1,1]
                # Set the range of x and y values to fit: from 7 bins to the left of the photopeak, up to the end of the spectra
                fitR_x = list(bin_centers[PhotopeakR-3:])
                fitL_x = list(bin_centers[PhotopeakL-3:])
                fitR_y = list(valuesR[PhotopeakR-3:])
                fitL_y = list(valuesL[PhotopeakL-3:])

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

                '''CHECK ENOUGH EVENTS AFTER PHOTOPEAK CUT'''
                
                # Check at least 50
                if len(temp) >= 50:
                    filled_lor_counts.append(len(temp))
                    if len(temp) > 200:
                        temp = np.asarray(random.sample(list(temp), 200))
                    coinc_list[lrc] = []
                    filled_indicators[lrc] = 1
                    count_filled += 1
                    # Add kept events to filled lors
                    filled_lors.append(20000) # 20000 arbitrary number greater than 10000 used for indicator
                    filled_lors.append(int(tolrt(lrc, lorchunk) // 3072)) # IDL 
                    filled_lors.append(int(tolrt(lrc, lorchunk) % 3072)) # IDR
                    for timediff in temp[:,0]:
                        filled_lors.append(int(timediff))
                    print(f'successfully filled L:{int(tolrt(lrc, lorchunk) // 3072)} R:{int(tolrt(lrc, lorchunk) % 3072)}')
            
            # Do a cut past 1000 events
            elif len(coinc_list[lrc]) > 1000:
                filled_lor_counts.append(1000)
                temp = np.asarray(random.sample(list(np.asarray(coinc_list[lrc])), 200))
                coinc_list[lrc] = []
                filled_indicators[lrc] = 1
                count_filled += 1
                filled_lors.append(20000) # 20000 arbitrary number greater than 10000 used for indicator
                filled_lors.append(int(tolrt(lrc, lorchunk) // 3072)) # IDL 
                filled_lors.append(int(tolrt(lrc, lorchunk) % 3072)) # IDR
                for timediff in temp[:,0]:
                    filled_lors.append(int(timediff))
                print(f'successfully filled L:{int(tolrt(lrc, lorchunk) // 3072)} R:{int(tolrt(lrc, lorchunk) % 3072)} from overflow')

    print(f'{round(time.time() - start, 3)}, {count_loops} loops, {count_filled} filled')

    # Datafile to keep filled LORs
    filled_lors = np.asarray(filled_lors, dtype = np.int16)
    filled_lors_file = data_storage_dir + 'filled_lors_file_' + str(lorchunk).zfill(2) + '.bin'
    if os.path.exists(filled_lors_file):
        with open(filled_lors_file, 'ab') as f:
            filled_lors.tofile(f)
    else:
        filled_lors.tofile(filled_lors_file)
    print(f'{round(time.time() - start, 3)}, filled_lors written')

    # Datafile to keep incomplete LORs
    out_data = []
    for lrc in range(lorchunksize):
        if not filled_indicators[lrc]:
            out_data.append([20000 + len(coinc_list[lrc]), tolrt(lrc, lorchunk) // 3072, tolrt(lrc, lorchunk) % 3072])
            out_data.extend(coinc_list[lrc])
            #out_info.append([tolrt(lrc, lorchunk) // 3072, tolrt(lrc, lorchunk) % 3072, len(coinc_list[lrc])])
        else:
            #out_info.append([tolrt(lrc, lorchunk) // 3072, tolrt(lrc, lorchunk) % 3072, 'FILLED'])
            pass
    out_data = np.asarray(out_data, dtype = np.float32)
    if os.path.exists(incomplete_lors_file):
        os.remove(incomplete_lors_file)
        out_data.tofile(incomplete_lors_file)
    else:
        out_data.tofile(incomplete_lors_file)
    print(f'{round(time.time() - start, 3)}, incomplete_lors written')

    # Datafile to keep track of filled indicators
    if os.path.exists(filled_indicators_file):
        os.remove(filled_indicators_file)
        filled_indicators.tofile(filled_indicators_file)
    else:
        filled_indicators.tofile(filled_indicators_file)
    print(f'{round(time.time() - start, 3)}, filled_indicators written')
    
    break

# Export information 
#out_df = pd.DataFrame(out_info)
#out_df.to_csv(data_storage_dir + 'out_info.csv', sep = ',', index = False, header = False)
#print(f'{time.time() - start}, info_out written')

# LOR filling completion histogram
#events_col = list(zip(*out_info))[2][1:]
#incomplete_lor_counts = []
#for item in events_col:
#    if item != 'FILLED':
#        incomplete_lor_counts.append(item)
#print(f'{time.time() - start}', len(incomplete_lor_counts), len(filled_lor_counts), len(incomplete_lor_counts) + len(filled_lor_counts))

#plt.figure(1)
#plt.hist(incomplete_lor_counts, bins = np.linspace(0, 2000, 200), color = 'red', alpha = 0.4)
#plt.hist(filled_lor_counts, bins = np.linspace(0, 2000, 200), color = 'blue', alpha = 0.4)
#plt.grid()
#plt.title('Distribution of LOR populations')
#plt.xlabel('Events')
#plt.ylabel('Number of LORs')
#plt.show()

