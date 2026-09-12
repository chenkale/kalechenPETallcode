from ChannelPairHeader import *
import pandas as pd
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import axes3d, Axes3D
import math
import matplotlib
import statistics as stat
pd.set_option("display.max_rows",1000)

# For each CP, we need to find the total number of coincidences it saw and divide it by the average number of coincidence over all CPs
# This factor should then be written out as the normalization factor for that particular CP (LOR)

#set this variable to False if running on TACC, else this is a good thing to have on.
debug = True
toReco = False

file_name = "Output_Combined_VerticalMotion_HighIntSource_800sec_coinc"
#file_name = "Output_NormalizationRun8_210mmSep_33mmSourceHeight_coinc"
#file_name = "Test_Data"
dir = "NewOutputs/MiniPET/" # Data file directory, if different from working directory
data = genfromtxt(dir+"{}.csv".format(file_name), delimiter="\t")
# Input data column order is ChanID left, ChanID right, PP_num_evts, PP_pos left, PP_pos right, ERes left, ERes right, PP_thrush_cut, PP_occ_cut


# Get the average occupancy accross all CPs
AvgHits = np.mean(data[:,2]) # Take the average over the last column of data in UniqueCPs, the column of CP occupancy
print(" *** Average PP Occupancy = " + str(AvgHits))

# Get normalization constants by dividing CP occupancy by the average occupancy accross al CPs
NormVals = data[:,2]/AvgHits
#print(NormVals)

# Add normalization constants to array of unique CPs
data = np.hstack((data,NormVals.reshape(len(data[:,0]), 1)))
#print(data)

# Change from AbsChanID to NewGeoID (or to RecoID)
if (toReco):
    data[:,0] = data[:,0] 
    data[:,1] = data[:,1] 
    data[:,0] = np.vectorize(toRecoChannelID)(data[:,0])
    data[:,1] = np.vectorize(toRecoChannelID)(data[:,1])
else:
    data[:,0] = data[:,0] - 896 # Port 8 in FEB/D, "Left" module when facing in beam direction (downstream) -> PCB #43
    data[:,1] = data[:,1] # Port 1 in FEB/D, "Right" module when facing in beam direction (downstream) -> PCB #4
    data[:,0] = np.vectorize(toNewGeoChannelID)(data[:,0])
    data[:,1] = np.vectorize(toNewGeoChannelID)(data[:,1])
#print(data)

# Sort by left and then right channelIDs
LOR_Norm = data[:,[0,1,9]]
ind = np.lexsort((LOR_Norm[:,1],LOR_Norm[:,0]))
#print(LOR_Norm[ind])

LOR_Norm = LOR_Norm[ind]
LOR_Norm = np.asarray(LOR_Norm)

# Drop rows where norm came out to zero
df = pd.DataFrame(LOR_Norm, columns = ['ChanIDL','ChanIDR','Norm'])
df = df[df['Norm'] != 0]
print(df)

# Write the output to the desired .csv file
outputdir = "NewOutputs/MiniPET/"
if (toReco):
    np.savetxt(outputdir+"LOR_Reco_Norm_Consts.csv", df, delimiter = "\t" )
else:
    np.savetxt(outputdir+"LOR_Geo_Norm_Consts.csv", df, delimiter = "\t" )


# Output to sorted debug file for more easily checking results
LOR_Norm_Sorted = LOR_Norm[LOR_Norm[:,2].argsort()]
df2 = pd.DataFrame(LOR_Norm_Sorted, columns = ['ChanIDL','ChanIDR','Norm'])
df2 = df2.astype({'ChanIDL':'int32','ChanIDR':'int32'})
df2 = df2[df2['Norm'] != 0]
print(df2)
np.savetxt(dir+"Normalizer_Debug.txt", df2, fmt = '%03i %03i %1.6f', delimiter = "\t")
