#Will Matava, Kyle Klein, John Cesar, etc. 

import os
import glob
from ChannelPairIncludes import *
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import axes3d, Axes3D
import math
import matplotlib
import statistics as stat
#from tqdm import tqdm
import time

# UNCOMMENT FOR FIRST PART OF PROGRAM
# matplotlib.use('Agg')

plt.style.use('ChannelPairStyles.mplstyle')
pd.set_option("display.max_columns",100)

#converts PETSys ID to geometric ID
def toGeo(x):
    y = 8*indices.get(x)[0] + indices.get(x)[1]
    return y

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


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

geo_channels = pd.DataFrame(columns = ["AbsChannelID","GeoChannelID"])
for i in range(128):
    geo_channels = pd.concat([geo_channels, pd.DataFrame({"AbsChannelID": [i], "GeoChannelID": [toGeo(i)]})],ignore_index=True)

def toGeoChannelID(AbsChannelID):
    #old
    # AbsChannelID = AbsChannelID - 131072*2
    slaveID = AbsChannelID // 4096
    chipID = (AbsChannelID - slaveID*4096) // 64
    channelID = AbsChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    AbsPCB_ChanID = geo_channels[geo_channels[:, 0] == PCB_ChanID][0][1]
    return 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64

# Read in from one data file
dir = "NewOutputs/MiniPET/"
df = pd.read_csv(dir+"Output_VerticalMotion_HighIntSource_800sec_coinc.csv", sep="\t", low_memory=False)
#df.columns = ["mhn1", "mhn2", "ppcount1", "ERmean", "ELmean", "ERres", "ELres", "CTR", "true", "true"] 
df.columns = ["ChannelIDL", "ChannelIDR", "PPCount", "ERmean", "ELmean", "ERres", "ELres", "CTR", "E_Cut", "PP_Cut"] 
# # merging the files
# joined_files = os.path.join("./", "*.csv")  
# # A list of all joined files is returned
# joined_list = glob.glob(joined_files)
# #print(joined_list)
# # Finally, the files are joined
# df = pd.concat(map(pd.read_csv, joined_list), ignore_index=True)
# print(df)


#Try using this
df["GeoChannelIDL"] = df["ChannelIDL"].apply(lambda x: toGeoChannelID(x))
df["GeoChannelIDR"] = df["ChannelIDR"].apply(lambda x: toGeoChannelID(x))

df["CTR"] = df["CTR"].abs()

#only keeping the data that passed the cuts
data = df[df["E_Cut"] == 1].reset_index(drop=True)
data = data[data["PP_Cut"] == 1].reset_index(drop=True)
data = data[data["CTR"] >= 0.1].reset_index(drop=True)
ChipData = getChipData(data,1,10)
#data = df

cmap = plt.get_cmap('binary')
cmap = truncate_colormap(cmap, 0.2, 0.8)
Scatter(getChipData(data,1,10),"CTR","ERres")
plt.show()

#Total CTR of scanner
Hist1D(getChipData(data, 1, 10)[0].CTR, units = "ps", fit= gaussFit(data.CTR,num_bins=50))
plt.title("Full Scanner CTR",fontsize=15)
plt.show()

# Best CTR for the whole scanner
CTRL = data.groupby("GeoChannelIDL").CTR.min().reset_index()
CTRL.columns = ["ChannelIDL", "CTR"]
CTRR = data.groupby("GeoChannelIDR").CTR.min().reset_index()
CTRR.columns = ["ChannelIDL", "CTR"]
CTR_b = pd.concat([CTRL,CTRR])

Hist1D(CTR_b.CTR,units='ps', fit= gaussFit(CTR_b.CTR,num_bins=30))
plt.title("LOR Optimized \"Best\" CTR")
plt.show()


ChipScatter(ChipData, "ERres", OneToOne=True, weighted=False)
plt.show()
ChipScatter(ChipData, "ERres", OneToOne=False, weighted=True)
plt.show()
ChipScatter(ChipData, "CTR", OneToOne=True, weighted=False)
plt.show()
ChipScatter(ChipData, "CTR", OneToOne=False, weighted=True)
plt.show()

#Coincidence Channel Sums
PP_P = pd.DataFrame()
# THIS IS NOT A MISTAKE THE TEXT PLACEMENT IS LEGACY CODE THAT REQUIRES PETSYS CHANNEL ID
PP_P["ChannelIDL"] = data["ChannelIDL"]
PP_P["ChannelIDR"] = data["ChannelIDR"]
PP_P["Count"] = data["PPCount"]
PP_L = PP_P.groupby("ChannelIDL").Count.sum().reset_index()
PP_L.columns = ["ChannelIDL", "Counts"]
PP_R = PP_P.groupby("ChannelIDR").Count.sum().reset_index()
PP_R.columns = ["ChannelIDR", "Counts"]

# plotting the photopeak counts
T_y = [9000, 6000,2000]
plt.hist((getChipData(data, 1, 10)[0].PPCount), bins=np.linspace(73,140,23))
plt.title("Number of Coincident Photoelectric Events in a Channel Pair: \n")
plt.xlabel("Number of Coincident Events", fontsize=size)
plt.ylabel("Counts", fontsize=size)
plt.yscale("log")
# plt.xscale("log")
plt.xlim(0,200)
#first bin goes to 13500
plt.ylim(1,10000)
T_x = generateTextX(6,150)
T_y = generateTextY(13,2000)
plt.text(T_x[0],T_y[0], "Mean: "+format(df.PPCount.mean(), "6.2f"))
plt.text(T_x[0],T_y[7], "Std: "+format(df.PPCount.std(), "5.2f"))
plt.text(T_x[0],T_y[12], "Count: "+format(df.PPCount.count(), "6"))
plt.show()

#Charge Weighted CTR on a heatmap
wCTR = pd.DataFrame(columns=["ChannelID", "wCTR"])
for i in range(0,128):
    run = data[data["GeoChannelIDR"] == i]
    run["CTR"] = 1000*run["CTR"]
    if run.ChannelIDR.count() == 0:
        wCTR = pd.concat([wCTR, pd.DataFrame({"ChannelID": [i], "wCTR":[0]})], ignore_index=True)
        #wCTR = wCTR.append({"ChannelID":i, "wCTR":0},ignore_index=True)
    else:
        run["CTR_weight"] = run["CTR"]*run["PPCount"]
        wCTR_val = run.CTR_weight.sum()/run.PPCount.sum()
        wCTR = pd.concat([wCTR, pd.DataFrame({"ChannelID": [i], "wCTR":[wCTR_val]})], ignore_index=True)
