from ChannelPairHeader import *
# def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
#     new_cmap = colors.LinearSegmentedColormap.from_list(
#         'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
#         cmap(np.linspace(minval, maxval, n)))
#     return new_cmap

import sys

# sys.argv[1] is the first command-line argument
#data_file = sys.argv[1]

#data = genfromtxt(data_file, delimiter="\t", usecols=(2,3,4,7,8,9))
#set this variable to False if running on TACC, else this is a good thing to have on.
debug = False
#"HWTriggerTest9.2FullScanner10min44_Slave1Chip5.dat"
#file_name = input("Filename: ")

data_dir = 'C:/Users/burri/Downloads/'
file_name = 'comp'
data = genfromtxt(data_dir + file_name + '.dat', delimiter="\t", usecols=(2,3,4,7,8,9))

#slow but correct way to remove any consecutive coincidences within 10 ns of eachother
# newdata = []
# i = 1
# while i < (len(data) - 1):
#     if abs(data[i][0] - data[i+1][0]) < 5000:
#         i += 2
#     elif abs(data[i][0] - data[i-1][0]) < 5000:
#         i += 1
#     else:
#         np.append(newdata,data[i,:], axis = 0)
#         i+=1
# print(len(data))
# print(len(newdata))

# Input data column order is TimeL, ChargeL, ChannelIDL, TimeR, ChargeR, ChannelIDR

# Count the number of unique channel pairs and keep only those with more than 300  events
CPIDs = data[:,[2,5]] # Take a slice of "data" containing just channel pairs (ChannelIDL and ChannelIDR)
UniqueCPs = np.unique(CPIDs, axis=0, return_counts=True) # Creates new array of all the unique channel pairs
CP_num = len(UniqueCPs[1]) # Number of unique channel pair
UniqueCPs = np.hstack((UniqueCPs[0],UniqueCPs[1].reshape(CP_num, 1)))
# plt.hist(UniqueCPs[:,2], bins=np.linspace(-0.5, 500.5, 102))
# plt.show()

num_counts = 60 # Threshold for number of events we want per channel pair (LOR)
#print(UniqueCPs)
UniqueCPs = UniqueCPs[UniqueCPs[:,2] > num_counts] # Keep only channel pairs with occupancy above threshold
if debug:
    print(" *** Of the {} original channle pairs, {} survived the population cut of {} events ***".format(CP_num,len(UniqueCPs),num_counts))

# The almighty data frame to be filled in the main loop below
DATA = []

#counting the number of prompt gammas
prompt = 0

# Loop over the post-cut unique channel pairs to perform analysis on
for i in tqdm(range(len(UniqueCPs[:,1]))):
    try:
        # Define channel_pair array as subset of data matching the UniqueCPs array
        channel_pair = data[data[:,2] == UniqueCPs[i][0]]
        channel_pair = channel_pair[channel_pair[:,5] == UniqueCPs[i][1]]
        LChannel = UniqueCPs[i][0]
        RChannel = UniqueCPs[i][1]
        # Generate the charge spectrum histograms
        valuesL, valuesR, bin_centers = HistogramCharge(channel_pair)
        if debug:
            print(" Looking for photopeaks in channel pair {} {}".format(toGeoChannelID(LChannel),toGeoChannelID(RChannel)))
        # Estimate the photopeak locations in each spectrum
        PhotopeakL, PhotopeakR = findPhotopeak(valuesL, valuesR, debug)
        # Use the above estimate to fit the photopeaks and return the FWHM as EFitParam
        EFitParam = fitPhotopeak(valuesL, valuesR, bin_centers, PhotopeakL, PhotopeakR)
        if debug:
            # Plot and save the fitted energy spectra for the current channel pair
            PlotEnergySpectrum(channel_pair, EFitParam, LChannel, RChannel, display=True)

        # Check if the photopeaks were found and are both above the threshold value
        if PPThresholdCut(EFitParam):
            if debug:
                print("     Photopeak is above threshold, cutting on events within 2.5 sigma of mean")
            # If so, keep only events above 2.5 sigma below the photopeak mean
            channel_pair = CutOnEnergy(EFitParam,channel_pair)
            #print(channel_pair)
            # Now check that photopeaks contain enough events, defined by the occupancy threshold below
            PP_num_evts = num_counts * 0.1
            time_diff = np.subtract(channel_pair[:,0], channel_pair[:,3])
            if time_diff.size >= PP_num_evts:
    #            print(time_diff)
                time_diff.sort()
    #           print(time_diff)
                CTRFitParam = fitTimeDiff(time_diff)
    #            print(CTRFitParam[2])
                PP_Cut, time_diff_data = PPOccupancyCut(PP_num_evts, channel_pair, CTRFitParam)
                if PP_Cut and time_diff_data.size >= PP_num_evts:
                    if debug:
                        print("        Sufficient statistics (counts > {}) in photopeak".format(PP_num_evts))
                    # Now ready to perform the CTR fit on the remaining events and extract the CTR parameter
                    CTRFitParam = fitTimeDiff(time_diff_data)
                    if debug:
                        # Plot and save the fitted CTR spectra for the current channel pair
                        PlotCTR(time_diff_data, CTRFitParam, LChannel, RChannel, display=True)
                    # Output the relevant channel pair info to the output data frame
                    # Output lines vary depending on which cuts were passed
                    # case 1: both cuts passed - outputline = left channel ID, right channel ID, num_evts, Photpeak Left, Photopeak Right, ERes left, ERes right, CTR, true, true
                    line = [LChannel, RChannel, len(channel_pair), EFitParam[1][1],EFitParam[0][1], 100*(2.355*EFitParam[1][2] / EFitParam[1][1]), 100*(2.355*EFitParam[0][2] / EFitParam[0][1]), 2.355*CTRFitParam[2], 1, 1]
                    DATA.append(line)
                # If there weren't enough events in the photopeaks...
                else:
                    # case 2: first cut passed - outputline = left channel ID, right channel ID, num_evts, ERes left, ERes right, true, false
                    line = [LChannel, RChannel, len(channel_pair), EFitParam[1][1],EFitParam[0][1],100*(2.355*EFitParam[1][2] / EFitParam[1][1]), 100*(2.355*EFitParam[0][2] / EFitParam[0][1]),0, 1, 0]
                    DATA.append(line)
            else:
                # case 3: first cut failed - outputline = left channel ID, right channel ID, num_evts, ERes left, ERes right, true, false
                line = [LChannel, RChannel, len(channel_pair), EFitParam[1][1],EFitParam[0][1],100*(2.355*EFitParam[1][2] / EFitParam[1][1]), 100*(2.355*EFitParam[0][2] / EFitParam[0][1]),0, 1, 0]
                DATA.append(line)
            # If the photopeaks were not above threshold...
        elif not PPThresholdCut(EFitParam):
            # case 4: no cuts passed - outputline = left channel ID, right channel ID, num_evts, ERes left, ERes right, false, false
            line = [LChannel, RChannel, len(channel_pair), EFitParam[1][1],EFitParam[0][1],100*(2.355*EFitParam[1][2] / EFitParam[1][1]), 100*(2.355*EFitParam[0][2] / EFitParam[0][1]),0, 0, 0]
            DATA.append(line)
    except:
        print('skipped')

# Write the output to the desired .csv file
DATA = np.asarray(DATA)
np.savetxt(data_dir + file_name + '_test.csv', DATA, delimiter = "\t" )

#np.savetxt(data_file[:-4] + '_test.csv', DATA, delimiter = "\t" )
