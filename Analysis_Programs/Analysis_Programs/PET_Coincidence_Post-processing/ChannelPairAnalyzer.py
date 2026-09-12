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

geo_channels = []
for i in range(128):
    geo_channels.append([i,toGeo(i)])
geo_channels = np.asarray(geo_channels)

def toGeoChannelID(AbsChannelID):
    # Convert PETSys absolute channel IDs to geomteric IDs
    portID = AbsChannelID // 131072
    slaveID = (AbsChannelID - 131072*portID) // 4096
    chipID = (AbsChannelID - slaveID*4096 - 131072*portID) // 64
    channelID = AbsChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    AbsPCB_ChanID = geo_channels[geo_channels[:,0] == PCB_ChanID][0][1]

    #General formula can be found in above function "to AbsChannelID"
    GeoChannelID = 10**6 * portID + 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64
    return GeoChannelID



# Read in from one data file
OV = 3.5
Thresh = 15
# dir = "/home/kilroy101/Experimental_Data/TPPT_DATA/Scanner_Commissioning/High_Intensity_Source/Threshold Optimization/OV = {}/".format(OV)
# df = pd.read_csv(dir+"Output_OV_{}_T1_{}.csv".format(OV, Thresh), sep="\t", low_memory=False)

dir = "/home/kilroy101/Experimental_Data/TPPT_DATA/Scanner_Commissioning/Threshold_Optimization/"
df = pd.read_csv(dir+"Output_LowIntensity2.5hour.csv", sep="\t", low_memory=False)
df.columns = ["ChannelIDL", "ChannelIDR", "PPCount", "ERmean", "ELmean", "ERres", "ELres", "CTR", "E_Cut", "PP_Cut"]


#Try using this
df["GeoChannelIDL"] = df["ChannelIDL"].apply(lambda x: toGeoChannelID(x))
df["GeoChannelIDR"] = df["ChannelIDR"].apply(lambda x: toGeoChannelID(x))

df["CTR"] = df["CTR"].abs()

# df = df.loc[(df["GeoChannelIDL"] >= 1011300) &  (df["GeoChannelIDL"] <= 1011400)]

# #cutting out all the geometrically localized outliers
# badindicesL = df.loc[(df["GeoChannelIDL"] >= 1011000) & (df["GeoChannelIDL"] <= 1011363)].index
# df.drop(badindicesL, inplace=True)
# df = df.reset_index(drop=True)
#
# badindicesR = df.loc[(df["GeoChannelIDR"] >= 11000) & (df["GeoChannelIDR"] <= 11363)].index
# df.drop(badindicesR, inplace=True)
# df = df.reset_index(drop=True)

#only keeping the data that passed the cuts
data = df[df["E_Cut"] == 1].reset_index(drop=True)
data = data[data["PP_Cut"] == 1].reset_index(drop=True)
data = data[data["CTR"] >= 100].reset_index(drop=True)
# data = data[data["CTR"] <= 400].reset_index(drop=True)
ChipData = getChipData(data,1,10)
#data = df
FullScannerHeatmap(data)

cmap = plt.get_cmap('binary')
cmap = truncate_colormap(cmap, 0.2, 0.8)
# Scatter(getChipData(data,1,10),"CTR","ERres")
# plt.show()

data["CWCTR"] = data["CTR"]*data["PPCount"]
print("Charge Weighted CTR: ", (data["CWCTR"].sum() / data["PPCount"].sum()))
#Total CTR of scanner
Hist1D(data.CTR, units = "ps",fit = gaussFit(data.CTR, num_bins=160))
plt.title("Full Scanner CTR") #OV = {} T1 = {} \n 35 mm shield High Intenstiy Source".format(OV, Thresh))
plt.xlim(150,400)
# T_x = generateTextX(6,150)
# T_y = generateTextY(14,400)
# plt.text(T_x[0],T_y[0], "Mean: "+format(df.PPCount.mean(), "6.2f"))
# plt.text(T_x[0],T_y[1], "Std: "+format(df.PPCount.std(), "5.2f"))
# plt.text(T_x[0],T_y[2], "Count: "+format(df.PPCount.count(), "6"))
plt.show()

# CTR for a single PCB
#ShowCTRArrayHist(data,0)

# # Best CTR by pixel for a PCB
# ShowBestCTRArrayHist(data,0, 0)

# # Best CTR for the whole scanner
# CTRL = data.groupby("GeoChannelIDL").CTR.min().reset_index()
# CTRL.columns = ["ChannelIDL", "CTR"]
# CTRR = data.groupby("GeoChannelIDR").CTR.min().reset_index()
# CTRR.columns = ["ChannelIDL", "CTR"]
# CTR_b = pd.concat([CTRL,CTRR])
#
# Hist1D(CTR_b.CTR,units='ps', fit= gaussFit(CTR_b.CTR,num_bins=30))
# plt.title("LOR Optimized \"Best\" CTR")
# plt.show()
#
#
# ChipScatter(ChipData, "ERres", OneToOne=True, weighted=False)
# plt.show()
# ChipScatter(ChipData, "ERres", OneToOne=False, weighted=True)
# plt.show()
# ChipScatter(ChipData, "CTR", OneToOne=True, weighted=False)
# plt.show()
# ChipScatter(ChipData, "CTR", OneToOne=False, weighted=True)
# plt.show()

# #ERes for the whole scanner
# ShowTotalERes(data)
# #ERes by array
# ShowEResArray(data,0)

# #Best ERes for a PCB
# EResB_data, EResBcolumn = getChipData(data,0,0)
# EResBT = EResB_data[EResB_data["Geo"+EResBcolumn] < 64].reset_index(drop=True)
# EResBT = EResBT.groupby("Geo"+EResBcolumn)["E_Res"+EResBcolumn[-1]].min().reset_index()
# EResBT.columns = [EResBcolumn, "ERes"]
# EResBB = EResB_data[EResB_data["Geo"+EResBcolumn] >= 64].reset_index(drop=True)
# EResBB = EResBB.groupby("Geo"+EResBcolumn)["E_Res"+EResBcolumn[-1]].min().reset_index()
# EResBB.columns = [EResBcolumn, "ERes"]
# ArrayHist(EResBT.ERes, EResBB.ERes, ChipID = 0, units = "%", fitT= gaussFit(EResBT.ERes, num_bins=15), fitB = gaussFit(EResBB.ERes,num_bins=15))
# plt.show()

# #Pixel map of coincident channel pairs
# for i in range(len(a)):
#     fig = plt.figure()
#     ax = fig.add_subplot(221)
#     ax1 = fig.add_subplot(222)
#     ax2 = fig.add_subplot(223)
#     ax3 = fig.add_subplot(224)
#     ax.set_title("PCB 6")
#     ax1.set_title("PCB 5")
#     for k in range(0,9):
#         ax.hlines(k,0,8)
#     for j in range(0,9):
#         ax.vlines(j,0,8)
#     if a[i] < 64:
#         ax.text(indices.get(a[i])[1]+0.05,8-indices.get(a[i])[0]-0.7, str(tot_phot[i]))

#     for k in range(0,9):
#         ax2.hlines(k,0,8)
#     for j in range(0,9):
#         ax2.vlines(j,0,8)
#     if a[i] >= 64:
#         ax2.text(indices.get(a[i])[1]+0.05,16-indices.get(a[i])[0]-0.7, str(tot_phot[i]))

#     for k in range(0,9):
#         ax1.hlines(k,0,8)
#     for j in range(0,9):
#         ax1.vlines(j,0,8)
#     if b[i] < 64:
#         ax1.text(indices.get(b[i])[1]+0.05,8-indices.get(b[i])[0]-0.7, str(tot_phot[i]))

#     for k in range(0,9):
#         ax3.hlines(k,0,8)
#     for j in range(0,9):
#         ax3.vlines(j,0,8)
#     if b[i] >= 64:
#         ax3.text(indices.get(b[i])[1]+0.05,16-indices.get(b[i])[0]-0.7, str(tot_phot[i]))
#     plt.savefig("Coincidence_Pair_{}.png".format(i))
#     plt.clf()
#     plt.close()



# Coincidence Channel Sums
# PP_P = pd.DataFrame()
# # THIS IS NOT A MISTAKE THE TEXT PLACEMENT IS LEGACY CODE THAT REQUIRES PETSYS CHANNEL ID
# PP_P["ChannelIDL"] = data["ChannelIDL"]
# PP_P["ChannelIDR"] = data["ChannelIDR"]
# PP_P["Count"] = data["PPCount"]
# PP_L = PP_P.groupby("ChannelIDL").Count.sum().reset_index()
# PP_L.columns = ["ChannelIDL", "Counts"]
# PP_R = PP_P.groupby("ChannelIDR").Count.sum().reset_index()
# PP_R.columns = ["ChannelIDR", "Counts"]
# #
#
# fig = plt.figure()
# ax = fig.add_subplot(221)
# ax1 = fig.add_subplot(222)
# ax2 = fig.add_subplot(223)
# ax3 = fig.add_subplot(224)
# ax.set_title("PCB 6")
# ax1.set_title("PCB 5")
#
#
# for k in range(0,9):
#     ax.hlines(k,0,8)
#     ax1.hlines(k,0,8)
#     ax2.hlines(k,0,8)
#     ax3.hlines(k,0,8)
# for j in range(0,9):
#     ax.vlines(j,0,8)
#     ax1.vlines(j,0,8)
#     ax2.vlines(j,0,8)
#     ax3.vlines(j,0,8)
#
# for i in range(PP_L.ChannelIDL.count()):
#     if PP_L.ChannelIDL[i] < 64:
#         ax.text(indices.get(PP_L.ChannelIDL[i])[1]+0.1,8-indices.get(PP_L.ChannelIDL[i])[0]-0.6, str(PP_L.Counts[i]))
#     elif PP_L.ChannelIDL[i] >= 64:
#         ax2.text(indices.get(PP_L.ChannelIDL[i])[1]+0.1,16-indices.get(PP_L.ChannelIDL[i])[0]-0.6, str(PP_L.Counts[i]))
#
# for i in range(PP_R.ChannelIDR.count()):
#     if PP_R.ChannelIDR[i] < 64:
#         ax1.text(indices.get(PP_R.ChannelIDR[i])[1]+0.1,8-indices.get(PP_R.ChannelIDR[i])[0]-0.6, str(PP_R.Counts[i]))
#     if PP_R.ChannelIDR[i] >= 64:
#         ax3.text(indices.get(PP_R.ChannelIDR[i])[1]+0.1,16-indices.get(PP_R.ChannelIDR[i])[0]-0.6, str(PP_R.Counts[i]))
#
# for i in range(df.ChannelIDL.count()):
#     if df.E_Cut[i] == 0:
#         if df.ChannelIDL[i] < 64:
#             ax.plot(indices.get(df.ChannelIDL[i])[1]+0.1,8-indices.get(df.ChannelIDL[i])[0]-0.6,'bs')
#         elif df.ChannelIDL[i] >= 64:
#             ax2.plot(indices.get(df.ChannelIDL[i])[1]+0.1,16-indices.get(df.ChannelIDL[i])[0]-0.6,'bs')
#     if df.PP_Cut[i] == 0:
#         if df.ChannelIDL[i] < 64:
#             ax.plot(indices.get(df.ChannelIDL[i])[1]+0.1,8-indices.get(df.ChannelIDL[i])[0]-0.6,'rs')
#         elif df.ChannelIDL[i] >= 64:
#             ax2.plot(indices.get(df.ChannelIDL[i])[1]+0.1,16-indices.get(df.ChannelIDL[i])[0]-0.6,'rs')
# #its not finished yet,don't yell at me
# for i in range(df.ChannelIDL.count()):
#     if df.E_Cut[i] == 0:
#         if df.ChannelIDR[i] < 64:
#             ax1.plot(indices.get(df.ChannelIDR[i])[1]+0.1,8-indices.get(df.ChannelIDR[i])[0]-0.6,'bs')
#         elif df.ChannelIDR[i] >= 64:
#             ax3.plot(indices.get(df.ChannelIDR[i])[1]+0.1,16-indices.get(df.ChannelIDR[i])[0]-0.6,'bs')
#     if df.PP_Cut[i] == 0:
#         if df.ChannelIDR[i] < 64:
#             ax1.plot(indices.get(df.ChannelIDR[i])[1]+0.1,8-indices.get(df.ChannelIDR[i])[0]-0.6,'rs')
#         elif df.ChannelIDR[i] >= 64:
#             ax3.plot(indices.get(df.ChannelIDR[i])[1]+0.1,16-indices.get(df.ChannelIDR[i])[0]-0.6,'rs')
#
# plt.show()
# # plt.savefig("Coincidence_Sums.png")
# plt.clf()
# plt.close()
#


# #Charge Integration
# fig = plt.figure(figsize=(12,8))
# ax = fig.add_subplot(211)
# ax1 = fig.add_subplot(212)
# ax.set_title("Reference PCB Charge Integral by Pixel")
# ax1.set_title("PCB 8 Charge Integral by Pixel")
# ax.set_ylim(0,6000)
# ax1.set_ylim(0,6000)
# #ax.scatter(E_Int.ChannelID, E_Int.Charge_INTR)
# #ax1.scatter(E_Int.ChannelID, E_Int.Charge_INTL)
# plt.subplots_adjust(hspace=0.5)
# plt.show()





# # plotting the photopeak counts
# T_y = [9000, 6000,2000]
# plt.hist((getChipData(data, 1, 10)[0].PPCount), bins=np.linspace(73,140,23))
# plt.title("Number of Coincident Photoelectric Events in a Channel Pair: \n")
# plt.xlabel("Number of Coincident Events", fontsize=size)
# plt.ylabel("Counts", fontsize=size)
# plt.yscale("log")
# # plt.xscale("log")
# plt.xlim(0,200)
# #first bin goes to 13500
# plt.ylim(1,10000)
# T_x = generateTextX(6,150)
# T_y = generateTextY(13,2000)
# plt.text(T_x[0],T_y[0], "Mean: "+format(df.PPCount.mean(), "6.2f"))
# plt.text(T_x[0],T_y[7], "Std: "+format(df.PPCount.std(), "5.2f"))
# plt.text(T_x[0],T_y[12], "Count: "+format(df.PPCount.count(), "6"))
# plt.show()

# ShowPPCount((getChipData(data, 1, 10)[0].PPCount), bins=np.linspace(73,140,23))
# plt.title("Number of Coincident Photoelectric Events in a Channel Pair: \n")
# plt.xlim(0,200)
# #first bin goes to 13500
# plt.ylim(1,10000)
# plt.show



# fig = plt.figure()
# ax = fig.add_subplot(221)
# ax1 = fig.add_subplot(222)
# ax2 = fig.add_subplot(223)
# ax3 = fig.add_subplot(224)
# ax.set_title("PCB 6")
# ax1.set_title("PCB 5")
#
# for k in range(0,9):
#     ax.hlines(k,0,8)
#     ax1.hlines(k,0,8)
#     ax2.hlines(k,0,8)
#     ax3.hlines(k,0,8)
# for j in range(0,9):
#     ax.vlines(j,0,8)
#     ax1.vlines(j,0,8)
#     ax2.vlines(j,0,8)
#     ax3.vlines(j,0,8)
#
# for i in range(LQ.Quality.count()):
#     if LQ.GeoChannelIDL[i] < 64:
#         # ax.text((LQ.GeoChannelIDL[i]%8)+0.2,8-(LQ.GeoChannelIDL[i]//8)-0.6, str(LQ.Quality[i]))
#         if LQ.Quality[i] == "A":
#             ax.plot((LQ.GeoChannelIDL[i]%8)+0.4,8-(LQ.GeoChannelIDL[i]//8)-0.6,'gs', ms=10)
#         elif LQ.Quality[i] == "B":
#             ax.plot((LQ.GeoChannelIDL[i]%8)+0.4,8-(LQ.GeoChannelIDL[i]//8)-0.6,'bs', ms=10)
#         elif LQ.Quality[i] == "C":
#             ax.plot((LQ.GeoChannelIDL[i]%8)+0.4,8-(LQ.GeoChannelIDL[i]//8)-0.6,'ys', ms=10)
#         elif LQ.Quality[i] == "D":
#             ax.plot((LQ.GeoChannelIDL[i]%8)+0.4,8-(LQ.GeoChannelIDL[i]//8)-0.6,'rs', ms=10)
#
#
#     elif LQ.GeoChannelIDL[i] >= 64:
#         # ax2.text((LQ.GeoChannelIDL[i]%8)+0.4,16-(LQ.GeoChannelIDL[i]//8)-0.6, str(LQ.Quality[i]))
#         if LQ.Quality[i] == "A":
#             ax2.plot((LQ.GeoChannelIDL[i]%8)+0.4,16-(LQ.GeoChannelIDL[i]//8)-0.6,'gs', ms=10)
#         elif LQ.Quality[i] == "B":
#             ax2.plot((LQ.GeoChannelIDL[i]%8)+0.4,16-(LQ.GeoChannelIDL[i]//8)-0.6,'bs', ms=10)
#         elif LQ.Quality[i] == "C":
#             ax2.plot((LQ.GeoChannelIDL[i]%8)+0.4,16-(LQ.GeoChannelIDL[i]//8)-0.6,'ys', ms=10)
#         elif LQ.Quality[i] == "D":
#             ax2.plot((LQ.GeoChannelIDL[i]%8)+0.4,16-(LQ.GeoChannelIDL[i]//8)-0.6,'rs', ms=10)
#
# for i in range(RQ.Quality.count()):
#     if RQ.GeoChannelIDR[i] < 64:
#         # ax1.text((RQ.GeoChannelIDR[i]%8)+0.4,8-(RQ.GeoChannelIDR[i]//8)-0.6, str(RQ.Quality[i]))
#         if RQ.Quality[i] == "A":
#             ax1.plot((RQ.GeoChannelIDR[i]%8)+0.4,8-(RQ.GeoChannelIDR[i]//8)-0.6,'gs', ms=10)
#         elif RQ.Quality[i] == "B":
#             ax1.plot((RQ.GeoChannelIDR[i]%8)+0.4,8-(RQ.GeoChannelIDR[i]//8)-0.6,'bs', ms=10)
#         elif RQ.Quality[i] == "C":
#             ax1.plot((RQ.GeoChannelIDR[i]%8)+0.4,8-(RQ.GeoChannelIDR[i]//8)-0.6,'ys', ms=10)
#         elif RQ.Quality[i] == "D":
#             ax1.plot((RQ.GeoChannelIDR[i]%8)+0.4,8-(RQ.GeoChannelIDR[i]//8)-0.6,'rs', ms=10)
#
#     elif RQ.GeoChannelIDR[i] >= 64:
#         # ax3.text((RQ.GeoChannelIDR[i]%8)+0.4,16-(RQ.GeoChannelIDR[i]//8)-0.6, str(RQ.Quality[i]))
#         if RQ.Quality[i] == "A":
#             ax3.plot((RQ.GeoChannelIDR[i]%8)+0.4,16-(RQ.GeoChannelIDR[i]//8)-0.6,'gs', ms=10)
#         elif RQ.Quality[i] == "B":
#             ax3.plot((RQ.GeoChannelIDR[i]%8)+0.4,16-(RQ.GeoChannelIDR[i]//8)-0.6,'bs', ms=10)
#         elif RQ.Quality[i] == "C":
#             ax3.plot((RQ.GeoChannelIDR[i]%8)+0.4,16-(RQ.GeoChannelIDR[i]//8)-0.6,'ys', ms=10)
#         elif RQ.Quality[i] == "D":
#             ax3.plot((RQ.GeoChannelIDR[i]%8)+0.4,16-(RQ.GeoChannelIDR[i]//8)-0.6,'rs', ms=10)
#
# plt.show()
# # plt.savefig("Coincidence_Sums.png")
# plt.clf()
# plt.close()

# plt.scatter(data.PPCount,data.CTR, s=5)
# plt.title("CTR vs Number of Photopeak events in a Channel Pair")
# plt.ylabel("CTR in ps")
# plt.xlabel("Photopeak Counts")
# plt.show()




# ###Single PCB evaluation graphs###
# plt.rcParams.update({'axes.labelsize': 'x-large'})
# plt.rcParams.update({'axes.titlesize': 'x-large'})

# #Scatter Plot of all CTRs
# test = data.copy()
# test["CTR"] = 1000*test["CTR"]
# fig = plt.figure()
# ax1 = fig.add_subplot(111, ylim=(175,375),xlim=(0,128), xlabel="ChannelID", ylabel="CTR in ps", title="CTRs PCB {}".format(PCB_num))#ylim = (100,300)
# ax1.scatter(test.GeoChannelIDR, test.CTR, marker="_")
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.show()

# #Charge Weighted CTR on a heatmap
# wCTR = pd.DataFrame(columns=["ChannelID", "wCTR"])
# for i in range(0,128):
#     run = data[data["GeoChannelIDR"] == i]
#     run["CTR"] = 1000*run["CTR"]
#     if run.ChannelIDR.count() == 0:
#         wCTR = pd.concat([wCTR, pd.DataFrame({"ChannelID": [i], "wCTR":[0]})], ignore_index=True)
#         #wCTR = wCTR.append({"ChannelID":i, "wCTR":0},ignore_index=True)
#     else:
#         run["CTR_weight"] = run["CTR"]*run["PPCount"]
#         wCTR_val = run.CTR_weight.sum()/run.PPCount.sum()
#         wCTR = pd.concat([wCTR, pd.DataFrame({"ChannelID": [i], "wCTR":[wCTR_val]})], ignore_index=True)
        #wCTR = wCTR.append({"ChannelID":i, "wCTR":wCTR_val}, ignore_index=True)
#print(wCTR)

# CTRRT = wCTR[wCTR["ChannelID"] < 64]
# CTRRB = wCTR[wCTR["ChannelID"] >= 64]

# CTRRT = np.reshape(np.array(CTRRT.wCTR), (8,8))
# CTRRB = np.reshape(np.array(CTRRB.wCTR), (8,8))

# plt.figure()
# plt.subplot(121)
# plt.xticks(np.arange(8))
# plt.yticks(np.arange(8))
# for i in range(8):
#     for j in range(8):
#         plt.text(j, 7-i, format(int(CTRRT[i][j]),"3"), ha="center", va="center", color="b", fontsize='large')
# plt.title("PCB {} Top".format(PCB_num))
# plt.imshow(CTRRT,cmap=cmap,vmin=150, vmax=320)

# plt.subplot(122)
# plt.xticks(np.arange(8))
# plt.yticks(np.arange(8))
# for i in range(8):
#     for j in range(8):
#         plt.text(j, 7-i, format(int(CTRRB[i][j]),"3"), ha="center", va="center", color="b", fontsize='large')
# plt.title("PCB {} Bottom".format(PCB_num))
# plt.imshow(CTRRB,cmap=cmap,vmin=150, vmax=320)

# # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
# cax = plt.axes([0.95, 0.1, 0.025, 0.8])
# plt.colorbar(cax=cax)
# plt.show()


# #Scatter Plot of all ERes
# test = data.copy()
# fig = plt.figure()
# ax1 = fig.add_subplot(111, ylim=(4,15),xlim=(0,128), xlabel="ChannelID", ylabel="ERes in %", title="Energy Resolutions PCB {}".format(PCB_num))#ylim=(4,10)
# ax1.scatter(test.GeoChannelIDR, test.ERres, marker="_")
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.show()
# plt.close()

# # Charge Weighted ERes on a heatmap
# wERes = pd.DataFrame(columns=["ChannelID", "wERes"])
# for i in range(0,128):
#     run = data[data["GeoChannelIDR"] == i]
#     if run.ChannelIDR.count() == 0:
#         wERes = wERes.append({"ChannelID":i, "wERes":0},ignore_index=True)
#     else:
#         run["ERres_weight"] = run["ERres"]*run["PPCount"]
#         wERes_val = run.ERres_weight.sum()/run.PPCount.sum()
#         wERes = wERes.append({"ChannelID":i, "wERes":wERes_val}, ignore_index=True)

# CTRRT = wERes[wERes["ChannelID"] < 64]
# CTRRB = wERes[wERes["ChannelID"] >= 64]

# CTRRT = np.reshape(np.array(CTRRT.wERes), (8,8))
# CTRRB = np.reshape(np.array(CTRRB.wERes), (8,8))

# plt.figure()
# plt.subplot(121)
# plt.xticks(np.arange(8))
# plt.yticks(np.arange(8))
# for i in range(8):
#     for j in range(8):
#         plt.text(j, 7-i, format(CTRRT[i][j],"4.1f"), ha="center", va="center", color="b", fontsize='x-large')
# plt.title("PCB {} Top".format(PCB_num))
# plt.imshow(CTRRT,cmap=cmap,vmin=4, vmax=15)

# plt.subplot(122)
# plt.xticks(np.arange(8))
# plt.yticks(np.arange(8))
# for i in range(8):
#     for j in range(8):
#         plt.text(j, 7-i, format(CTRRB[i][j],"4.1f"), ha="center", va="center", color="b", fontsize='x-large')
# plt.title("PCB {} Bottom".format(PCB_num))
# plt.imshow(CTRRB,cmap=cmap,vmin=4, vmax=15)

# # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
# cax = plt.axes([0.95, 0.1, 0.025, 0.8])
# plt.colorbar(cax=cax)
# plt.show()

# # Charge Weighted Photopeak Location on a heatmap
# wPP = pd.DataFrame(columns=["ChannelID", "wPP"])
# for i in range(0,128):
#     run = data[data["GeoChannelIDR"] == i]
#     if run.ChannelIDR.count() == 0:
#         wPP = wPP.append({"ChannelID":i, "wPP":0},ignore_index=True)
#     else:
#         run["PP_weight"] = run["PeakR"]*run["PPCount"]
#         wPP_val = run.PP_weight.sum()/run.PPCount.sum()
#         wPP = wPP.append({"ChannelID":i, "wPP":wPP_val}, ignore_index=True)

# CTRRT = wPP[wPP["ChannelID"] < 64]
# CTRRB = wPP[wPP["ChannelID"] >= 64]

# CTRRT = np.reshape(np.array(CTRRT.wPP), (8,8))
# CTRRB = np.reshape(np.array(CTRRB.wPP), (8,8))

# plt.figure()
# plt.subplot(121)
# plt.xticks(np.arange(8))
# plt.yticks(np.arange(8))
# for i in range(8):
#     for j in range(8):
#         plt.text(j, 7-i, format(CTRRT[i][j],"4.1f"), ha="center", va="center", color="b", fontsize='x-large')
# plt.title("PCB {} Top".format(PCB_num))
# plt.imshow(CTRRT,cmap=cmap,vmin=10, vmax=40)

# plt.subplot(122)
# plt.xticks(np.arange(8))
# plt.yticks(np.arange(8))
# for i in range(8):
#     for j in range(8):
#         plt.text(j, 7-i, format(CTRRB[i][j],"4.1f"), ha="center", va="center", color="b", fontsize='x-large')
# plt.title("PCB {} Bottom".format(PCB_num))
# plt.imshow(CTRRB,cmap=cmap,vmin=10, vmax=40)

# # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
# cax = plt.axes([0.95, 0.1, 0.025, 0.8])
# plt.colorbar(cax=cax)
# plt.show()
