from ChannelPairHeader import *
# def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
#     new_cmap = colors.LinearSegmentedColormap.from_list(
#         'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
#         cmap(np.linspace(minval, maxval, n)))
#     return new_cmap

#set this variable to False if running on TACC, else this is a good thing to have on.
debug = False
#file_name = "NormalizationRun14_210mmSep_45mmSourceHeight_coinc"
file_name = "VerticalMotion_HighIntSource_800sec_coinc.dat"
file_name_2 = "VerticalMotion_HighIntSource_800sec_2_coinc.dat"
dir = "NewRuns/MiniPET/" # Data file directory, if different from working directory
def read_large_file(file_path, usecols=(2,3,4,7,8,9), dtype=np.float64):
    return np.loadtxt(file_path, delimiter="\t", usecols=usecols, dtype=dtype)
data1 = read_large_file(dir + file_name)
data2 = read_large_file(dir + file_name_2)

data = np.vstack((data1, data2))
# Input data column order is TimeL, ChargeL, ChannelIDL, TimeR, ChargeR, ChannelIDR

# Count the number of unique channel pairs and keep only those with more than 300  events
CPIDs = data[:,[2,5]] # Take a slice of "data" containing just channel pairs (ChannelIDL and ChannelIDR)
UniqueCPs = np.unique(CPIDs, axis=0, return_counts=True) # Creates new array of all the unique channel pairs
CP_num = len(UniqueCPs[1]) # Number of unique channel pair
UniqueCPs = np.hstack((UniqueCPs[0],UniqueCPs[1].reshape(CP_num, 1)))
#num_counts = 50 # Threshold for number of events we want per channel pair (LOR)
#print(UniqueCPs)
#UniqueCPs = UniqueCPs[UniqueCPs[:,2] > num_counts] # Keep only channel pairs with occupancy above threshold
#if debug:
#    print(" *** Of the {} original channle pairs, {} survived the population cut of {} events ***".format(CP_num,len(UniqueCPs),num_counts))

# The almighty data frame to be filled in the main loop below
# DATA = []

#counting the number of prompt gammas
prompt = 0
def reduce_channel_pairs(UniqueCPs, debug=False):
    """
    Reduces the number of unique channel pairs by lowering the threshold 
    until at most 1000 channel pairs remain.
    """
    num_counts = 100 # Initial threshold
    temp = [1 for i in range(1)]
    while len(temp) < 1000 and num_counts > 30:
        temp = UniqueCPs[UniqueCPs[:, 2] > num_counts]
        num_counts -= 10 # Reduce threshold by 100
        print(len(temp))


    print(num_counts)
    return temp

UniqueCPs = reduce_channel_pairs(UniqueCPs)
DATA = []
# BAD_CHANNELS = []
f=0

# Loop over the post-cut unique channel pairs to perform analysis on
for i in tqdm(range(len(UniqueCPs[:,1]))):
    # Define channel_pair array as subset of data matching the UniqueCPs array
    channel_pair = data[data[:,2] == UniqueCPs[i][0]]
    channel_pair = channel_pair[channel_pair[:,5] == UniqueCPs[i][1]]
    LChannel = UniqueCPs[i][0]
    RChannel = UniqueCPs[i][1]
    # Generate the charge spectrum histograms
    valuesL, valuesR, bin_centers = HistogramCharge(channel_pair)
    if debug:
        print(" Looking for photopeaks in channel pair {} {}".format(LChannel,RChannel))
    # Estimate the photopeak locations in each spectrum
    PhotopeakL, PhotopeakR = findPhotopeak(valuesL, valuesR, debug)
    # print(f'PhotopeakL: {PhotopeakL}, \nPhotopeakR: {PhotopeakR}')
    # Use the above estimate to fit the photopeaks and return the FWHM as EFitParam
    EFitParam = fitPhotopeak(valuesL, valuesR, bin_centers, PhotopeakL, PhotopeakR)
    if debug:
        # Plot and save the fitted energy spectra for the current channel pair
        plot_energy_spectrum(channel_pair, EFitParam, LChannel, RChannel, display=False)


    # Check if the photopeaks were found and are both above the threshold value
    if PPThresholdCut(EFitParam):
        if debug:
            print("     Photopeak is above threshold, cutting on events within 2.5 sigma of mean")
        # If so, keep only events above 2.5 sigma below the photopeak mean
        channel_pair = CutOnEnergy(EFitParam,channel_pair)
#        print(channel_pair)

        # Now check that photopeaks contain enough events, defined by the occupancy threshold below
        PP_num_evts = 10
        # If the photopeaks are above threshold, and have enough stats (above PP_num_evts)...
        time_diff = np.subtract(channel_pair[:,0], channel_pair[:,3])
        if time_diff.size >= PP_num_evts:
            # print(time_diff)
            if debug:
                print("        Sufficient statistics (counts > {}) in photopeak".format(PP_num_evts))
            time_diff.sort()
        # If there weren't enough events in the photopeaks...
        else:
            if debug:
                print("        Insufficient statistics (counts < {}) in photopeak".format(PP_num_evts))
            # BAD_CHANNELS.append(channel_pair)
            continue
        CTRFitParam = fitTimeDiff(time_diff)
        PP_Cut, time_diff_data = PPOccupancyCut(PP_num_evts, channel_pair, CTRFitParam)
        if PP_Cut and time_diff_data.size >= PP_num_evts:
            filtered_arr = np.array([x for x in time_diff if x !=0])
        
            line = [LChannel, RChannel, len(channel_pair), EFitParam[1][1],EFitParam[0][1],EFitParam[1][2] , EFitParam[0][2],CTRFitParam[1], CTRFitParam[0]]
            if debug:
                plot_energy_spectrum(channel_pair,  EFitParam, LChannel, RChannel, display=True)
                plot_ctr(time_diff, CTRFitParam, LChannel, RChannel, display=True)
            DATA.append(line)
    # If the photopeaks were not above threshold...
    elif not PPThresholdCut(EFitParam):
        # case 3: no cuts passed - outputline = ChanID left, ChanID right, PP_num_evts, PP left, PP right, ERes left, ERes right, false, false
        line = [LChannel, RChannel, len(channel_pair), EFitParam[1][1],EFitParam[0][1],100*(2.355*EFitParam[1][2] / EFitParam[1][1]), 100*(2.355*EFitParam[0][2] / EFitParam[0][2]), 0, 0]
        DATA.append(line)
        # BAD_CHANNELS.append(channel_pair)

# Write the output to the desired .csv file
DATA = np.asarray(DATA)
# BAD_CHANNELS = np.asarray(BAD_CHANNELS)
np.savetxt(dir+"Output_Combined_{}.csv".format(file_name), DATA, delimiter = "\t" )
# np.savetxt(dir+"Output_BAD_{}.csv".format(file_name), BAD_CHANNELS, delimiter = "\t")
