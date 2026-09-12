import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm




dir = "/home/kilroy101/path/Geant4/Expirimental_Data/Scanner_Commissioning/"
#convert the sub array channel ID to the PETsys convention
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

#converts PETSys ID to geometric ID for a single PCB of 128
def toGeo(x):
    x_channel = x % 128
    y = 8*indices.get(x_channel)[0] + indices.get(x_channel)[1]
    return y


geo_channels = pd.DataFrame(columns=["AbsChannelID", "GeoChannelID"])
for i in range(128):
    geo_channels = geo_channels.append({"AbsChannelID": i, "GeoChannelID": toGeo(i)},ignore_index=True)

#general form to convert Reco Channel ID to Geo Channel ID
def toGeoChannelID(RecoChannelID):
    global geo_channels
    slaveID = RecoChannelID // 1024
    chipID = (RecoChannelID - slaveID*1024) // 64
    #in this case, it is equivlent to take mod 64 (more traditional) and mod 100 because no numbers above 63 and below 100 will ever appear
    channelID = RecoChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    conv = geo_channels.loc[geo_channels["AbsChannelID"] == PCB_ChanID].reset_index()
    GeoPCB_ChanID = conv.GeoChannelID[0]

    GeoChannelID = 10000*slaveID + 100*chipID + GeoPCB_ChanID % 64

    return GeoChannelID

#getting the row of the channel ID to only plot one row at a time on the histograms
def ZRow(ChipID):
    if ChipID < 8:
        Col = ChipID % 2
    elif ChipID > 7:
        Col = 2 + ChipID % 2
    return Col

#properly formatting the channelID values so they increment by one every time for the purposes of graphing
def XCol(ChannelID):
    Slave = ChannelID // 10000
    Chip = (ChannelID - Slave*10000) // 100
    Channel = ChannelID % 100
    if Slave >= 3 and Chip < 8:
        XCol = 32*(5 - Slave) + 8*(Chip // 2 ) + (Channel%8)
        return XCol
    elif Slave < 3 and Chip < 8:
        XCol = 32*(Slave) + 8*((Chip // 2) ) + (7 - (Channel%8))
        return XCol
    elif Slave >= 3 and Chip >= 8:
        XCol = 32*(5 - Slave) + 8*((15-Chip) // 2 ) + (Channel%8)
        return XCol
    else:
        XCol = 32*(Slave) + 8*((15-Chip) // 2 ) + (7 - (Channel%8))
        return XCol


#Check Channel ID Range
chunks = pd.read_csv(dir+"Center_Source_8-8_reco.lm", sep=",", header=None, names= ["ChannelIDL", "ChannelIDR"],usecols=[0, 1],low_memory=False,chunksize=100000)

#takes a little less than thirty seconds
concat_df = pd.concat(chunks,ignore_index=True)

fig = plt.figure()
fig1 = plt.figure()
ax = fig.add_subplot(111)
ax1 = fig1.add_subplot(111)

ax.hist(concat_df.ChannelIDL, bins = np.linspace(3072,6145,300), color="r")
ax1.hist(concat_df.ChannelIDR, bins = np.linspace(0,3073,300), color="b")
ax.set_title("Channel IDs present in the Left Crescent")
ax1.set_title("Channel IDs present in the Right Crescent")
ax.set_xlim(3072,6145)
ax1.set_xlim(0,3073)
plt.show()

# #Check geometry file
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/Scanner_Commissioning/Geometry/"
# df = pd.read_csv(dir+"TPPT_Scanner_map.csv", sep=",", header = None)
# df.columns = ["ChannelID", "x", "y", "z"]
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter3D(df.x, df.y, df.z)
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# ax.set_title("Position of Channels in Geometry lookup file")
# plt.show()


# #plots x y projections of lines for debugging
# df = pd.read_csv(dir+"Center_Source_8-8_reco.lm", sep=",", header=None, names= ["ChannelIDL", "ChannelIDR"],usecols=[0, 1],low_memory=False,nrows=20000)
# df_geo = pd.read_csv(dir+"Geometry/TPPT_Scanner_map.csv",sep=",", header=None, names=["ChannelID", "x", "y","z"], low_memory=False)
#
# #projections onto my XY plane
# def XYLine(x,ChannelIDL, ChannelIDR, geo_file):
#     x1 = geo_file.loc[geo_file["ChannelID"] == ChannelIDL].reset_index().x[0]
#     y1 = geo_file.loc[geo_file["ChannelID"] == ChannelIDL].reset_index().y[0]
#     x2 = geo_file.loc[geo_file["ChannelID"] == ChannelIDR].reset_index().x[0]
#     y2 = geo_file.loc[geo_file["ChannelID"] == ChannelIDR].reset_index().y[0]
#
#     m = (y2 - y1) / (x2 - x1)
#     b = y1 - m*x1
#     y = m*x + b
#     return [y,m]
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title("First 20,000 lines with conversion to Reco ID")
# for i in tqdm(range(df.ChannelIDL.count())):
#     x1 = df_geo.loc[df_geo["ChannelID"] == df.iloc[i]["ChannelIDL"]].reset_index().x[0]
#     x2 = df_geo.loc[df_geo["ChannelID"] == df.iloc[i]["ChannelIDR"]].reset_index().x[0]
#     x = np.linspace(x1, x2, 3)
#     ax.plot(x,XYLine(x,df.iloc[i]["ChannelIDL"],df.iloc[i]["ChannelIDR"],df_geo)[0], color = 'b', alpha=0.01)
# plt.show()



#Going column by column in the scanner to look at activity per channel ID
df = pd.read_csv(dir+"Center_Source_8-8_reco.lm", sep=",", header=None, names= ["ChannelIDL", "ChannelIDR"],usecols=[0, 1],low_memory=False, nrows=200000)

dfL = df.ChannelIDL.value_counts().reset_index()
dfL.columns = ["ChannelIDL", "Counts"]

dfR = df.ChannelIDR.value_counts().reset_index()
dfR.columns = ["ChannelIDR", "Counts"]


dfL["ChannelIDL"] = dfL["ChannelIDL"].apply(toGeoChannelID)
dfR["ChannelIDR"] = dfR["ChannelIDR"].apply(toGeoChannelID)

dfL["ChipIDL"] = (dfL["ChannelIDL"] - 10000*(dfL["ChannelIDL"]//10000)) // 100
dfR["ChipIDR"] = (dfR["ChannelIDR"] - 10000*(dfR["ChannelIDR"]//10000)) // 100

dfL["Col"] = dfL["ChipIDL"].apply(ZRow)
dfR["Col"] = dfR["ChipIDR"].apply(ZRow)

dfL["XCol"] = dfL["ChannelIDL"].apply(XCol)
dfR["XCol"] = dfR["ChannelIDR"].apply(XCol)


# #row by row across the crescent
# for i in range(32):
#     fig = plt.figure()
#     fig1 = plt.figure()
#     ax = fig.add_subplot(111)
#     ax1 = fig1.add_subplot(111)
#     sliceL = dfL.loc[(8*dfL["Col"] + ((dfL["ChannelIDL"] % 100) // 8)) == i]
#     sliceL = sliceL.sort_values(by="ChannelIDL").reset_index(drop=True)
#     ax.bar(sliceL.XCol, sliceL.Counts, align='center')
#     ax.set_xlabel("Crystal Column")
#     ax.set_ylabel("Counts")
#     ax.set_xlim(-1,96)
#     ax.set_title("Left Side Row: "+str(i))
#
#     sliceR = dfR.loc[(8*dfR["Col"] + ((dfR["ChannelIDR"] % 100) // 8)) == i]
#     sliceR = sliceR.sort_values(by="ChannelIDR").reset_index(drop=True)
#     ax1.bar(sliceR.XCol, sliceR.Counts, align='center')
#     ax1.set_xlabel("Crystal Column")
#     ax1.set_ylabel("Counts")
#     ax1.set_xlim(-1,96)
#     ax1.set_title("Right Side Row: "+str(i))
#     plt.show()


# #heat map
# dfL["HMapIndex"] = 96*(8*dfL["Col"] + ((dfL["ChannelIDL"] % 100) // 8)) + dfL["XCol"]
# dfR["HMapIndex"] = 96*(8*dfR["Col"] + ((dfR["ChannelIDR"] % 100) // 8)) + dfR["XCol"]
#
# #filling in dead pixels, 3072 = 32*96
# for i in range(3072):
#     if i not in dfL.HMapIndex.unique():
#         fillerL = pd.DataFrame({"ChannelIDL":[0], "Counts":[0], "ChipIDL":[0], "Col":[0], "XCol":[0], "HMapIndex":[i]})
#         dfL = pd.concat([dfL, fillerL], ignore_index=True)
#     if i not in dfR.HMapIndex.unique():
#         fillerR = pd.DataFrame({"ChannelIDR":[0], "Counts":[0], "ChipIDR":[0], "Col":[0], "XCol":[0], "HMapIndex":[i]})
#         dfR = pd.concat([dfR, fillerR], ignore_index=True)
#
# dfL = dfL.sort_values(by="HMapIndex").reset_index(drop=True)
# dfR = dfR.sort_values(by="HMapIndex").reset_index(drop=True)
#
# crescent_mapL = np.reshape(np.array(dfL.Counts), (32,96))
# crescent_mapR = np.reshape(np.array(dfR.Counts), (32,96))
# index_test = np.reshape(np.array(dfL.index), (32,96))
#
# # for i in range(32):
# #     for j in range(96):
# #         plt.text(j, i, crescent_mapL[i][j], ha="center", va="center", color="k", fontsize='small')
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# im = ax.imshow(crescent_mapL,vmin=0, vmax=dfL.Counts.max())
# cax = plt.axes([0.95, 0.1, 0.025, 0.8])
# fig.colorbar(im, cax=cax)
# ax.set_title("Left Crescent")
# plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
# plt.show()
#
#
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# im1 = ax1.imshow(crescent_mapR,vmin=0, vmax=dfR.Counts.max())
# cax = plt.axes([0.95, 0.1, 0.025, 0.8])
# fig1.colorbar(im1, cax=cax)
# ax1.set_title("Right Crescent")
# plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
# plt.show()
