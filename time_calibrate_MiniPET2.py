# Center source calibration data and miniPET geometry directories
center_source_data_dir = 'C:/Users/burri/Downloads/VerticalMotion_HighIntSource_800sec_coinc.dat'
miniPET_geo_map_dir = 'C:/Users/burri/Downloads/testMiniPETmap.csv' # Only used if generate_coincidence_heatmaps = True
output_dir = 'C:/Users/burri/Downloads/calibrationtest1'

# Generate time difference histograms
generate_time_difference_histograms = True
histogramlim = 10000
binwidth = 50
# Generate coincidence heatmaps for visual check
generate_coincidence_heatmaps = False
# Generate coincidence scatterplots for visual check
generate_coincidence_scatterplots = False


# Fitted bins for lookup values
fitted_bins = 200

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt



# Read calibration data into pd DataFrame
df = pd.DataFrame()
for chunk in pd.read_csv(center_source_data_dir, usecols=[2, 4, 7, 9], delimiter="\t", header=None, chunksize=10000):
    # Name columns
    chunk.columns = ['TimeL', 'IDL', 'TimeR', 'IDR']
    # Convert Left channel to mod 128
    chunk['IDL'] = chunk['IDL'] % 128
    # Time difference column, left minus right lets positive time differences be on the right side
    chunk['TimeDiff'] = chunk['TimeL'] - chunk['TimeR']
    
    df = pd.concat([df, chunk], ignore_index=True)

# Convert to numpy array, faster
arr = df.to_numpy()



# Geometric time offset matrix, stores time difference corresponding to the geometric center for each LOR
# MiniPET should be all 0
# each row corresponds to a left channel, each col to a right channel, such that each entry is a unique LOR
geometric_time_offsets = [[0 for i in range(128)] for j in range(128)]

# Read map
mapdf = pd.read_csv(miniPET_geo_map_dir, usecols=[1,2,3], header=None)
map = mapdf.to_numpy()

# For testing:
#for x in range(128):
#    map[x][0] = map[x][0] + 50

for l in range(128):
    for r in range(128):

        xL, yL, zL = map[l]
        xR, yR, zR = map[r+128]
        # Get midpoint between two channels
        mX, mY, mZ = (xL+xR)/2, (yL+yR)/2, (zL+zR)/2

        # If midpoint lies on x = 0, write to table
        if mX == 0:
            geometric_time_offsets[l][r] = 0
        # else need to calculate distance between x = 0 and midpoint along LOR
        else:
            # Parameterize, below are slopes with x being the parameter
            ay, az = (yR-yL)/(xR-xL), (zR-zL)/(xR-xL)
            # Find y and z for positive distance xL away from xL (x = 0)
            y0, z0 = yL + ay * -1 * xL, zL + az * -1 * xL
            # Distance, converted to light picoseconds
            d = ((mX)**2 + (y0-mY)**2 + (z0-mZ)**2)**0.5 / 0.299792458 #mm/ps
            # Find direction (midpoint left or right from x=0)
            if -1 * xL > xR:
                geometric_time_offsets[l][r] = -1 * d
            else:
                geometric_time_offsets[l][r] = d



# Initial setup

# Matrix for storing ARRAYS of coincidences, adjusted for geometric offset
# each row corresponds to a left channel, each col to a right channel, 
# so then each entry is a unique LOR 
coincidences_matrix = [[[] for i in range(128)] for j in range(128)]

# Creates coincidence matrix containing TimeDiffs (relative to x=0) for every LOR 
for row in arr:
    geometric_time_offset = geometric_time_offsets[int(row[1])][int(row[3])]
    coincidences_matrix[int(row[1])][int(row[3])].append(row[4] + geometric_time_offset)

# Array for storing left and right module channel time offsets, these will be used to calculate per channel offset values
coincidences_arr_l = [[] for i in range(128)]
coincidences_arr_r = [[] for i in range(128)]

# Fill arrays
for l in range(128):
    for r in range(128):
        coincidences_arr_l[l].extend(coincidences_matrix[l][r])
        coincidences_arr_r[r].extend(coincidences_matrix[l][r])

# Create lookup table for each channel in l and r
lookup_l = np.asarray([0.0 for i in range(128)])
lookup_r = np.asarray([0.0 for i in range(128)])

# Define gaussian fit
def gaussian(x, a, mu, c):
    return a * np.exp(-(x - mu)**2 / (2 * c**2))

def lookup(l, r):
    return -1 + lookup_l[l] + lookup_r[r]

# Updates coincidence matrix and left/right arrays
def update_coincidences():
    # Update coincidence matrix
    for l in range(128):
        for r in range(128):
            temp = []
            for t in coincidences_matrix[l][r]:
                temp.append(t + lookup(x[1], x[3]))
            coincidences_matrix[l][r] = temp

    # Update left/right arrays
    global coincidences_arr_l
    coincidences_arr_l = [[] for i in range(128)]
    global coincidences_arr_r
    coincidences_arr_r = [[] for i in range(128)]

    # Fill arrays
    for l in range(128):
        for r in range(128):
            coincidences_arr_l[l].extend(coincidences_matrix[l][r])
            coincidences_arr_r[r].extend(coincidences_matrix[l][r])

# Calibration function
def calibrate():
    # Left side
    for l in range(128):
        # Create histogram
        counts, bins = np.histogram(coincidences_arr_l[l], bins = fitted_bins)
        # Initial guess for amplitude
        inita = 0.8 * np.max(counts)
        # Initial guess for mean
        initmu = bins[:-1][np.where(counts == np.max(counts))[0][0]]
                    
        try:
            # try fit, with bad stats the fit sometimes doesn't work
            popt = curve_fit(gaussian, bins[:-1], counts, p0 = [inita, initmu, 500])[0]
            # Check that mean is reasonable 
            if abs(popt[1]) < 10000: # 10000
                # write time offset to table as negative so that it can be ADDED when used
                lookup_l[l] = -1 * popt[1] / 2
                #lookup_l[l] = -1 * popt[1] + lookup_l[l]
        except:
            pass

    # Right side
    for r in range(128):
        # Create histogram
        counts, bins = np.histogram(coincidences_arr_r[r], bins = fitted_bins)
        # Initial guess for amplitude
        inita = 0.8 * np.max(counts)
        # Initial guess for mean
        initmu = bins[:-1][np.where(counts == np.max(counts))[0][0]]
                    
        try:
            # try fit, with bad stats the fit sometimes doesn't work
            popt = curve_fit(gaussian, bins[:-1], counts, p0 = [inita, initmu, 500])[0]
            # Check that mean is reasonable 
            if abs(popt[1]) < 10000: # 10000
                # write time offset to table as negative so that it can be ADDED when used
                lookup_r[r] = popt[1] / 2
                #lookup_r[r] = -1 * popt[1] + lookup_r[r]
                #lookup_r[r] = popt[1] + lookup_r[r]
        except:
            pass

# We can try calibrating multiple times to tighten distribution, doens't seem to work
calibrate()
#update_coincidences()
#calibrate()


# Create master time difference array
master_time_diffs = []
for x in arr:
    master_time_diffs.append(x[4])

# Calibrated time difference array
master_time_diffs_calibrated = []
for x in arr:
    master_time_diffs_calibrated.append(x[0] + lookup_l[x[1]] - x[2] - lookup_r[x[3]])

# Prepare table to write to output file
output = pd.DataFrame()
for l in range(128):
    output = pd.concat([output, pd.DataFrame([[0, 0, l // 64, l % 64, lookup_l[l]]])], ignore_index=True)
for r in range(128):
    output = pd.concat([output, pd.DataFrame([[0, 0, r // 64 + 14, r % 64, lookup_r[r]]])], ignore_index=True)

# Write to csv
output.to_csv(output_dir + '/time_offset_calibration.tsv', sep='\t', header=None, float_format = '%.7e', index=False)


# Generate time difference histograms
if generate_time_difference_histograms:

    def formatlabel(a, mu, c):
        return f'a={round(a)}, mu={round(mu,1)}, sig={round(abs(c),1)}'

    # Histogram before time calibration
    hist_before = plt.hist(master_time_diffs, bins=range(int(min(master_time_diffs)), int(max(master_time_diffs)) + binwidth, binwidth), color = 'red', alpha = 0.4)
    bincenters0 = [(hist_before[1][x]+hist_before[1][x+1])/2 for x in range(len(hist_before[0]))]
    popt0 = curve_fit(gaussian, bincenters0, hist_before[0], p0 = [6000, 0, 1000])[0]
    #plt.plot(np.linspace(-1*histogramlim, histogramlim, 300), gaussian(np.linspace(-1*histogramlim, histogramlim, 300), *popt0), color = 'red', label = formatlabel(*popt0))
    plt.plot(np.linspace(-1*histogramlim, histogramlim, 300), gaussian(np.linspace(-1*histogramlim, histogramlim, 300), *popt0), color = 'red', label = 'before')
    print(popt0)

    # Histogram after time calibration
    hist_after = plt.hist(master_time_diffs_calibrated, bins=range(int(min(master_time_diffs_calibrated)), int(max(master_time_diffs_calibrated)) + binwidth, binwidth), color = 'blue', alpha = 0.4)
    bincenters1 = [(hist_after[1][x]+hist_after[1][x+1])/2 for x in range(len(hist_after[0]))]
    popt1 = curve_fit(gaussian, bincenters1, hist_after[0], p0 = [6000, 0, 1000])[0]
    #plt.plot(np.linspace(-1*histogramlim, histogramlim, 300), gaussian(np.linspace(-1*histogramlim, histogramlim, 300), *popt1), color = 'blue', label = formatlabel(*popt1))
    plt.plot(np.linspace(-1*histogramlim, histogramlim, 300), gaussian(np.linspace(-1*histogramlim, histogramlim, 300), *popt1), color = 'blue', label = 'after')
    print(popt1)

    plt.legend()
    plt.grid()
    plt.xlim(-1*histogramlim, histogramlim)
    plt.title('Time differences calibrated')
    plt.xlabel('Time difference (ps)')
    plt.ylabel('Counts')
    plt.savefig(output_dir + '/time_diffence_histograms.png')
    plt.close()

# Generate heatmaps before and after for visual check
if generate_coincidence_heatmaps or generate_coincidence_scatterplots:

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
                temp = getcoinclocation(l, r, timediff + lookup(l, r))
                coincxs1.append(temp[0])
                coincys1.append(temp[1])
                coinczs1.append(temp[2])

    if generate_coincidence_heatmaps:
        plt.hist2d(coincxs, coincys, range=[[-210, 210], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(output_dir + '/heatmap_before_xy.png')
        plt.close()
        plt.hist2d(coincxs, coinczs, range=[[-210, 210], [-15, 15]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.savefig(output_dir + '/heatmap_before_xz.png')
        plt.close()
        plt.hist2d(coincys, coinczs, range=[[-20, 20], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.savefig(output_dir + '/heatmap_before_yz.png')
        plt.hist2d(coincxs1, coincys1, range=[[-210, 210], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(output_dir + '/heatmap_after_xy.png')
        plt.close()
        plt.hist2d(coincxs1, coinczs1, range=[[-210, 210], [-15, 15]], bins=400, cmap='viridis')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.savefig(output_dir + '/heatmap_after_xz.png')
        plt.close()
        plt.hist2d(coincys1, coinczs1, range=[[-20, 20], [-20, 20]], bins=400, cmap='viridis')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.savefig(output_dir + '/heatmap_after_yz.png')
        plt.close()
    
    if generate_coincidence_scatterplots:
        plt.plot(coincxs, coincys, 'o', markersize = 0.5, color = 'black')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xlim(-20, 20)
        plt.ylim(-50, 50)
        plt.savefig(output_dir + '/scatterplot_before_xy.png')
        plt.close()
        plt.plot(coincxs, coinczs, 'o', markersize = 0.5, color = 'black')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.xlim(-20, 20)
        plt.ylim(-30, 30)
        plt.savefig(output_dir + '/scatterplot_before_xz.png')
        plt.close()
        plt.plot(coincxs, coincys, 'o', markersize = 0.5, color = 'black')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.xlim(-20, 20)
        plt.ylim(-20, 20)
        plt.savefig(output_dir + '/scatterplot_before_yz.png')
        plt.close()
        plt.plot(coincxs1, coincys1, 'o', markersize = 0.5, color = 'black')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xlim(-20, 20)
        plt.ylim(-50, 50)
        plt.savefig(output_dir + '/scatterplot_after_xy.png')
        plt.close()
        plt.plot(coincxs1, coincys1, 'o', markersize = 0.5, color = 'black')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.xlim(-20, 20)
        plt.ylim(-30, 30)
        plt.savefig(output_dir + '/scatterplot_after_xz.png')
        plt.close()
        plt.plot(coincxs, coincys, 'o', markersize = 0.5, color = 'black')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.xlim(-20, 20)
        plt.ylim(-20, 20)
        plt.savefig(output_dir + '/scatterplot_after_yz.png')
        plt.close()
