#program meant to provide graphs for on site quick evaluation of FLASH Data, not a comprehensive graph profile
#Kyle Klein, John Cesar
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import scipy as sp
from iminuit import cost, Minuit
import iminuit as im
import os
from tqdm import tqdm
plt.style.use('ChannelPairStyles.mplstyle')

ln2 = 0.6391

def generateTextY(num,max,subplots=False):
    if max > 50 and subplots==False:
        y = [max - 0.04*(i+1)*max for i in range(0,num+1)]
    elif max <= 50 and subplots == False:
        y = [max - 0.06*(i+1)*max for i in range(0,num+1)]
    elif subplots == True:
        y = [0.96*max - 0.08*(i+1)*max for i in range(0,num+1)]
    return y

def toGeo(x):
    #converts PETSys ID to geometric ID
    y = 8*indices.get(x)[0] + indices.get(x)[1]
    return y

def PlotEnergySpectrum(ChannelPairData, LChannel, RChannel, display=False):
    # plotting energy spectra for every channel pair and saving them all to .png files
    # Plot the two histograms and their fitted functions
    minL = 0
    maxL = max(ChannelPairData[:,4])
    minR = 0
    maxR = max(ChannelPairData[:,1])
    if max(ChannelPairData[:,4]) < 0:
        minL = max(ChannelPairData[:,4])
        maxL = 0
    if max(ChannelPairData[:,1]) < 0:
        minR = max(ChannelPairData[:,1])
        maxR = 0
    # print(minL,maxL,minR,maxR)
    x_max = max(maxR, maxL)# Hard coded value for the maximum x-range
    x_f = np.linspace(0,x_max,300) # Sets the binning: 300 bins from 0 to x_max
    #making sure the giant ass text is in the figure when it's saved
    fig = plt.figure(figsize=(12,9))
    try:
        #Changed from Edge Color ec to color and made alpha 0.4
        valuesR,binsR,stuffR = plt.hist(ChannelPairData[:,4], bins = np.linspace(minL,maxL,abs(int((160*max(ChannelPairData[:,4]))/45))), color='C0',fill=True,alpha=0.5, label='Channel ID {}'.format(int(toGeoChannelID(LChannel))))
        valuesL,binsL,stuffL = plt.hist(ChannelPairData[:,1], bins = np.linspace(minR,maxR,abs(int((160*max(ChannelPairData[:,1]))/45))), color='salmon',fill=True,alpha=0.5, label='Channel ID {}'.format(int(toGeoChannelID(RChannel))))
    except Exception:
        print("Couldn't plot this hist for some reason")
    # Code to help generate text overlays for plot info
    # Set the title and axis labels, draw the legend, and then save the figure to a .png
    plt.title("Energy Spectrum - LOR {}-{}".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
    plt.xlabel("Energy in DAQ units")
    plt.ylabel("Counts per unit charge")
    plt.legend()
    plt.savefig(save_dir+"Energy Spectrum - LOR {}-{}".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
    if display:
        plt.show()
        plt.clf()
        plt.close()
    else:
        plt.clf()
        plt.close()

indices = {
      0 : (4,7-7),
      1 : (4,7-6),
      2 : (7,7-5),
      3 : (5,7-7),
      4 : (5,7-4),
      5 : (5,7-5),
      6 : (4,7-4),
      7 : (7,7-7),
      8 : (6,7-6),
      9 : (7,7-4),
      10 : (5,7-6),
      11 : (6,7-4),
      12 : (4,7-5),
      13 : (6,7-5),
      14 : (6,7-7),
      15 : (7,7-6),
      16 : (3,7-7),
      17 : (3,7-6),
      18 : (2,7-7),
      19 : (2,7-6),
      20 : (0,7-7),
      21 : (1,7-7),
      22 : (0,7-6),
      23 : (1,7-6),
      24 : (3,7-5),
      25 : (1,7-5),
      26 : (2,7-5),
      27 : (4,7-3),
      28 : (0,7-5),
      29 : (3,7-4),
      30 : (0,7-4),
      31 : (1,7-4),
      32 : (2,7-4),
      33 : (3,7-3),
      34 : (2,7-3),
      35 : (0,7-3),
      36 : (1,7-3),
      37 : (0,7-2),
      38 : (5,7-3),
      39 : (1,7-2),
      40 : (2,7-2),
      41 : (3,7-2),
      42 : (1,7-1),
      43 : (0,7-1),
      44 : (0,7-0),
      45 : (3,7-1),
      46 : (1,7-0),
      47 : (2,7-1),
      48 : (3,7-0),
      49 : (2,7-0),
      50 : (6,7-2),
      51 : (6,7-1),
      52 : (7,7-1),
      53 : (4,7-1),
      54 : (5,7-1),
      55 : (6,7-0),
      56 : (7,7-0),
      57 : (7,7-2),
      58 : (7,7-3),
      59 : (4,7-2),
      60 : (5,7-0),
      61 : (5,7-2),
      62 : (6,7-3),
      63 : (4,7-0),
      64:(3+8,7),
      65:(3+8,6),
      66:(2+8,4),
      67:(2+8,6),
      68:(3+8,4),
      69:(1+8,7),
      70:(1+8,5),
      71:(0+8,7),
      72:(1+8,6),
      73:(3+8,3),
      74:(2+8,7),
      75:(2+8,3),
      76:(3+8,5),
      77:(0+8,5),
      78:(2+8,5),
      79:(0+8,6),
      80:(4+8,7),
      81:(6+8,7),
      82:(5+8,7),
      83:(7+8,7),
      84:(5+8,6),
      85:(4+8,6),
      86:(6+8,6),
      87:(7+8,6),
      88:(4+8,5),
      89:(6+8,5),
      90:(5+8,5),
      91:(1+8,4),
      92:(7+8,5),
      93:(7+8,4),
      94:(6+8,4),
      95:(4+8,4),
      96:(5+8,4),
      97:(5+8,3),
      98:(6+8,3),
      99:(4+8,3),
      100:(7+8,3),
      101:(7+8,2),
      102:(0+8,4),
      103:(6+8,2),
      104:(7+8,1),
      105:(5+8,2),
      106:(6+8,1),
      107:(4+8,2),
      108:(7+8,0),
      109:(5+8,1),
      110:(6+8,0),
      111:(4+8,1),
      112:(5+8,0),
      113:(4+8,0),
      114:(0+8,2),
      115:(2+8,1),
      116:(0+8,1),
      117:(3+8,1),
      118:(1+8,1),
      119:(1+8,0),
      120:(0+8,0),
      121:(1+8,2),
      122:(1+8,3),
      123:(3+8,2),
      124:(2+8,0),
      125:(2+8,2),
      126:(0+8,3),
      127:(3+8,0)}

geo_channels = []
for i in range(128):
    geo_channels.append([i,toGeo(i)])
geo_channels = np.asarray(geo_channels)

def toGeoChannelID(AbsChannelID):
    # Convert PETSys absolute channel IDs to geomteric IDs
    slaveID = AbsChannelID // 4096
    chipID = (AbsChannelID - slaveID*4096) // 64
    channelID = AbsChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    AbsPCB_ChanID = geo_channels[geo_channels[:,0] == PCB_ChanID][0][1]

    #General formula can be found in above function "to AbsChannelID"
    GeoChannelID = 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64
    return GeoChannelID


file_name = "HorizontalSource_AlignedWithLowerArray_MidSeparationV2_10min_coinc"
dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/FLASH_Prep_Runs/"
save_dir = dir+file_name[:-6]+"/"
try:
    os.mkdir(save_dir)
except FileExistsError:
    print("The directory for storing figures for this run already exists")

# file_name = "Phantom6_coinc"
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/"

#For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
data = np.genfromtxt(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
# Input data column order is TimeL, ChargeL, ChannelIDL, TimeR, ChargeR, ChannelIDR

data[:,0] = data[:,0] / 1000000000000
data[:,3] = data[:,3] / 1000000000000

#First we must determine when the first spill ocurred
num_bins = 2*int(input("How long was the run in s? "))
fig = plt.figure(figsize=(16,9))
values_time,bins_time,params_time = plt.hist(data[:,0], bins=num_bins)
# plt.xlim(-50,175)
# plt.yscale("log")
plt.xlabel("Time in s")
plt.ylabel("No. Coincidences per 500 ms")
# plt.title("Number of Coincidences detected as a function of time \n Ge-68 Source Data Run 3")
plt.savefig(save_dir+"Coincidences-over-time")
plt.show()

# Histogram of coincidences per channel
# Left module
np.vectorize(toGeoChannelID)(data[:,2])
fig = plt.figure(figsize=(16,9))
values,bins,params = plt.hist(np.vectorize(toGeoChannelID)(data[:,2]), bins=np.linspace(799.5,964.5,166))
plt.ylim(1,1.2*max(values))
plt.xlabel("Geometric Channel ID Left")
plt.ylabel("No. of Coincidences")
plt.title("Number of Coincidences per Channel - Left Module")
plt.savefig(save_dir+"Coincidences-per-ChannelLeft")
plt.show()
# Right module
fig = plt.figure(figsize=(16,9))
values,bins,params = plt.hist(np.vectorize(toGeoChannelID)(data[:,5]), bins=np.linspace(599.5,764.5,166))
plt.ylim(1,1.2*max(values))
plt.xlabel("Geometric Channel ID Right")
plt.ylabel("No. of Coincidences")
plt.title("Number of Coincidences per Channel - Right Module")
plt.savefig(save_dir+"Coincidences-per-ChannelRight")
plt.show()

#We can zoom in on the max bin in the above coarse coincidences over time plot and look at the spill
bin_centers_time = 0.5*(bins_time[1:] + bins_time[:-1])
max_bin = np.argmax(values_time)
#getting the time of the max value (presumably the spill)
time_spill = bin_centers_time[max_bin]
#only looking plus or minus 250 ms around the spill
spill_data = data[(data[:,0] >= time_spill-0.25) & (data[:,0] <= time_spill+0.25)]
#500 ms / 5000 = 100 microseconds
fig = plt.figure(figsize=(16,9))
plt.hist(spill_data[:,0],bins=5000)
plt.xlabel("Time [s]")
plt.ylabel("No. Coincidences per 100 $\mu$s")
plt.title("Coincidences over time: Zoomed around the spill")
plt.savefig(save_dir+"Coincidences-over-time-spill")
plt.show()


# Count the number of unique channel pairs and keep only those with more than 300  events
CPIDs = data[:,[2,5]] # Take a slice of "data" containing just channel pairs (ChannelIDL and ChannelIDR)
UniqueCPs = np.unique(CPIDs, axis=0, return_counts=True) # Creates new array of all the unique channel pairs
CP_num = len(UniqueCPs[1]) # Number of unique channel pair
UniqueCPs = np.hstack((UniqueCPs[0],UniqueCPs[1].reshape(CP_num, 1)))
#Sort the unique CPs by the number of coincidences in them
UniqueCPs = UniqueCPs[UniqueCPs[:,2].argsort()]
#reverse the ordering so the most populated CPs are at the begining of the list. God I love numpy indexing
UniqueCPs = UniqueCPs[::-1]

#recording the time of earliest detection within a channel pair
CP_times = []
for i in tqdm(range(len(UniqueCPs[:,1]))):
    # Define channel_pair array as subset of data matching the UniqueCPs array
    channel_pair = data[data[:,2] == UniqueCPs[i][0]]
    channel_pair = channel_pair[channel_pair[:,5] == UniqueCPs[i][1]]
    LChannel = UniqueCPs[i][0]
    RChannel = UniqueCPs[i][1]
    #recording the earliest detection time wtihin a channel pair
    CP_times.append(channel_pair[channel_pair[:,0].argsort()][0][0])

    if i <10:
        PlotEnergySpectrum(channel_pair, LChannel, RChannel, display=False)
fig = plt.figure(figsize=(16,9))
plt.hist(CP_times,bins=num_bins//2)
plt.xlabel("Time [s]")
plt.ylabel("No. of Channel Pairs [$s^{-1}$]")
plt.title("Time of first detection in each channel pair")
fig.tight_layout()
plt.savefig(save_dir+"First-Detection-Time",dpi=2*96)
plt.show()
