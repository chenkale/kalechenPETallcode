#Kyle Klein
#5-2-21
#Converts Position in TPPT scanner to both Geometric and PETsys ID
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from tqdm import tqdm

#for debugging, if you want to see all rows
pd.set_option('display.max_rows', None)

#definition of constants, distances in cm, angles in radians
pi = 3.14159265
radius = 167.540
angular_offset = 40.5 * pi / 180
starting_angle = 3*pi/2 - angular_offset

#for geometric layout, there will be 4 indices (3 in reality.
#The first index is a boolean index denoting the left, 0, or right, 1, that will determine the sign of the x coordinate. This is not a real index and is purely to make computation easier.
#The next one denotes the column of the module so 4 modules will have the same column, called ModCol
#The next index is the column within a module where column 0 will always have the smallest y coordinate, called ArrCol
#The next index is the row within a module, which directly correlates to the z coordinate, called ArrRow


def middle_column_x(ModCol):
    theta = starting_angle - (9*ModCol * pi / 180)
    x = radius*np.cos(theta)
    return x

def middle_column_y(ModCol):
    theta = starting_angle - (9*ModCol * pi / 180)
    y = radius*np.sin(theta)
    return y

#ModCol_theta denotes the angle from the y axis that the normal vector of the array is pointing
def ModCol_theta(ModCol):
    angle = abs((49.5 * pi / 180) - (9*ModCol * pi / 180))
    return angle



def nth_column_x(ModCol, ArrCol):
    if ModCol < 6:
        fourth_column_x = middle_column_x(ModCol) - 1.6*np.sin(ModCol_theta(ModCol))
        if ArrCol <= 4:
            x = fourth_column_x + 3.2*(4 - ArrCol)*np.sin(ModCol_theta(ModCol))
            return x
        elif ArrCol > 4:
            x = fourth_column_x - 3.2*(ArrCol - 4)*np.sin(ModCol_theta(ModCol))
            return x
    else:
        fourth_column_x = middle_column_x(ModCol) + 1.6*np.sin(ModCol_theta(ModCol))
        if ArrCol <= 4:
            x = fourth_column_x - 3.2*(4 - ArrCol)*np.sin(ModCol_theta(ModCol))
            return x
        elif ArrCol > 4:
            x = fourth_column_x + 3.2*(ArrCol - 4)*np.sin(ModCol_theta(ModCol))
            return x

def nth_column_y(ModCol, ArrCol):
    fourth_column_y = middle_column_y(ModCol) + 1.6*np.cos(ModCol_theta(ModCol))
    if ArrCol <= 4:
        y = fourth_column_y - 3.2*(4 - ArrCol)*np.cos(ModCol_theta(ModCol))
        return y
    elif ArrCol > 4:
        y = fourth_column_y + 3.2*(ArrCol - 4)*np.cos(ModCol_theta(ModCol))
        return y


# #testing if the above functions work
# cryst_x = []
# cryst_y = []
# #12 module columns
# for i in range(0,6):
#     #7 crystals within the array
#     for j in range(0,8):
#         cryst_x.append(nth_column_x(i,j))
#         cryst_y.append(nth_column_y(i,j))
#
# #displays a single row of the left crescent
# plt.scatter(cryst_x,cryst_y,s=10)
# #must keep square window or image will be distorted
# plt.xlim(-180,180)
# plt.ylim(-180,180)
# plt.show()

def nth_row_z(ArrRow):
    if ArrRow < 8:
        z = 1.6 + 3.2*ArrRow
        return z
    elif ArrRow >= 8 and ArrRow < 16:
        z = 1.6 + 3.2*7 + 4 + 3.2*(ArrRow - 8)
        return z
    elif ArrRow >= 16 and ArrRow < 24:
        z = 1.6 + 2*(3.2*7 + 4) + 3.2*(ArrRow - 16)
        return z
    elif ArrRow >= 24:
        z = 1.6 + 3*(3.2*7 + 4) + 3.2*(ArrRow - 24)
        return z




#Now that the geometric coordinates are defined for every crystal, we must begin matching them to their geometric channel IDs
#Geometric Channel ID is defined as slaveID*10^4 + chipID * 10^2 + channelID


module_dictL = {"0,0": 500, "0,1": 501, "0,2": 514, "0,3": 515,
"1,0": 502, "1,1": 503, "1,2": 512, "1,3": 513,
"2,0": 504, "2,1": 505, "2,2": 510, "2,3": 511,
"3,0": 506, "3,1": 507, "3,2": 508, "3,3": 509,

"4,0": 400, "4,1": 401, "4,2": 414, "4,3": 415,
"5,0": 402, "5,1": 403, "5,2": 412, "5,3": 413,
"6,0": 404, "6,1": 405, "6,2": 410, "6,3": 411,
"7,0": 406, "7,1": 407, "7,2": 408, "7,3": 409,

"8,0": 300, "8,1": 301, "8,2": 314, "8,3": 315,
"9,0": 302, "9,1": 303, "9,2": 312, "9,3": 313,
"10,0": 304, "10,1": 305, "10,2": 310, "10,3": 311,
"11,0": 306, "11,1": 307, "11,2": 308, "11,3": 309}

module_dictR = {"0,0": 0, "0,1": 1, "0,2": 14, "0,3": 15,
"1,0": 2, "1,1": 3, "1,2": 12, "1,3": 13,
"2,0": 4, "2,1": 5, "2,2": 10, "2,3": 11,
"3,0": 6, "3,1": 7, "3,2": 8, "3,3": 9,

"4,0": 100, "4,1": 101, "4,2": 114, "4,3": 115,
"5,0": 102, "5,1": 103, "5,2": 112, "5,3": 113,
"6,0": 104, "6,1": 105, "6,2": 110, "6,3": 111,
"7,0": 106, "7,1": 107, "7,2": 108, "7,3": 109,

"8,0": 200, "8,1":201, "8,2": 214, "8,3": 215,
"9,0": 202, "9,1": 203, "9,2": 212, "9,3": 213,
"10,0": 204, "10,1": 205, "10,2":210, "10,3": 211,
"11,0": 206, "11,1": 207, "11,2":208, "11,3": 209}



def GeoChannelID(LR,ModCol, ArrCol, ArrRow):
    global module_dictL
    global module_dictR
    if LR == 0:
        new_ArrCol = ArrCol
        new_ArrRow = 31 - ArrRow
        int_new_ArrRow = new_ArrRow // 8
        mod_new_ArrRow = new_ArrRow % 8
        # #old FEM layout
        # module_dict = {"0,0": 502, "0,1": 503, "0,2": 514, "0,3": 515,
        #             "1,0": 500, "1,1": 501, "1,2": 512, "1,3": 513,
        #             "2,0": 414, "2,1": 415, "2,2": 510, "2,3": 511,
        #             "3,0": 412, "3,1": 413, "3,2": 508, "3,3": 509,
        #             "4,0": 410, "4,1": 411, "4,2": 506, "4,3": 507,
        #             "5,0": 408, "5,1": 409, "5,2": 504, "5,3": 505,
        #             "6,0": 310, "6,1": 311, "6,2": 406, "6,3": 407,
        #             "7,0": 308, "7,1": 309, "7,2": 404, "7,3": 405,
        #             "8,0": 306, "8,1": 307, "8,2": 402, "8,3": 403,
        #             "9,0": 304, "9,1": 305, "9,2": 400, "9,3": 401,
        #             "10,0": 302, "10,1": 303, "10,2": 314, "10,3": 315,
        #             "11,0": 300, "11,1": 301, "11,2": 312, "11,3": 313}
        GeoChannelID = 100 *module_dictL.get(str(ModCol)+","+str(int_new_ArrRow)) + 8*mod_new_ArrRow + new_ArrCol
        return GeoChannelID
    else:
        new_ArrCol = 7 - ArrCol
        new_ArrRow = 31 - ArrRow
        int_new_ArrRow = new_ArrRow // 8
        mod_new_ArrRow = new_ArrRow % 8

        # #old FEM layout
        # module_dict = {"0,0": 0, "0,1": 1, "0,2": 12, "0,3": 13,
        #             "1,0": 2, "1,1": 3, "1,2": 14, "1,3": 15,
        #             "2,0": 4, "2,1": 5, "2,2": 100, "2,3": 101,
        #             "3,0": 6, "3,1": 7, "3,2": 102, "3,3": 103,
        #             "4,0": 8, "4,1": 9, "4,2": 104, "4,3": 105,
        #             "5,0": 10, "5,1": 11, "5,2": 106, "5,3": 107,
        #             "6,0": 108, "6,1": 109, "6,2": 204, "6,3": 205,
        #             "7,0": 110, "7,1": 111, "7,2": 206, "7,3": 207,
        #             "8,0": 112, "8,1":113, "8,2": 208, "8,3": 209,
        #             "9,0": 114, "9,1": 115, "9,2": 210, "9,3": 211,
        #             "10,0": 200, "10,1": 201, "10,2":212, "10,3": 213,
        #             "11,0": 202, "11,1": 203, "11,2":214, "11,3": 215}
        GeoChannelID = 100 *module_dictR.get(str(ModCol)+","+str(int_new_ArrRow)) + 8*mod_new_ArrRow + new_ArrCol
        return GeoChannelID

#The covnvesion between GeoChannelID and AbsChannelID

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

    #converts PETSys ID to geometric ID
def toGeo(x):
    x_channel = x % 128
    y = 8*indices.get(x_channel)[0] + indices.get(x_channel)[1]
    return y


geo_channels = pd.DataFrame(columns=["AbsChannelID", "GeoChannelID"])
for i in range(128):
    geo_channels = geo_channels.append({"AbsChannelID": i, "GeoChannelID": toGeo(i)},ignore_index=True)

def toAbsChannelID(GeoChannelID):
    slaveID = GeoChannelID // 10**4
    chipID = (GeoChannelID - slaveID*10000) // 10**2
    #in this case, it is equivlent to take mod 64 (more traditional) and mod 100 because no numbers above 63 and below 100 will ever appear
    channelID = GeoChannelID % 100

    PCB_ChanID = 64*(chipID % 2) + channelID
    conv = geo_channels.loc[geo_channels["GeoChannelID"] == PCB_ChanID].reset_index()
    AbsPCB_ChanID = conv.AbsChannelID[0]

    #General formula can be found in the software user's guide. All FEB/Ds are plugged into Port 2.
    AbsChannelID = 4096*slaveID + 64*chipID + AbsPCB_ChanID % 64

    return AbsChannelID

#The covnvesion between GeoChannelID and AbsChannelID
def toGeoChannelID(AbsChannelID):
    global geo_channels
    # old
    # AbsChannelID = AbsChannelID - 131072*2
    slaveID = AbsChannelID // 4096
    chipID = (AbsChannelID - slaveID*4096) // 64
    channelID = AbsChannelID % 64


    PCB_ChanID = 64*(chipID % 2) + channelID
    conv = geo_channels.loc[geo_channels["AbsChannelID"] == PCB_ChanID].reset_index()
    AbsPCB_ChanID = conv.GeoChannelID[0]

    #General formula can be found in above function "to AbsChannelID"
    GeoChannelID = 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64
    return GeoChannelID


def toRecoChannelID(AbsChannelID):
    slaveID = AbsChannelID // 4096
    chipID = (AbsChannelID - slaveID*4096) // 64
    channelID = AbsChannelID % 64
    RecoChannelID = 1024*slaveID + 64*chipID + channelID
    return RecoChannelID

def toGeoChannelIDfromRecoChannelID(RecoChannelID):
    global geo_channels
    slaveID = RecoChannelID // 1024
    chipID = (RecoChannelID - 1024*slaveID) // 64
    channelID = RecoChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    conv = geo_channels.loc[geo_channels["AbsChannelID"] == PCB_ChanID].reset_index()
    AbsPCB_ChanID = conv.GeoChannelID[0]

    #General formula can be found in above function "to AbsChannelID"
    GeoChannelID = 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64
    return GeoChannelID

#need to do GeoChannelID as well if you want to have that intuition
df = pd.DataFrame(columns=["RecoChannelID","x", "y", "z"])
for i in range(2):#2
    for j in tqdm(range(12)):#12
        for k in range(8):#8
            for m in range(32):#32
                if i == 0:
                    df = df.append({"RecoChannelID":toRecoChannelID(toAbsChannelID(GeoChannelID(i,j,k,m))),"x":nth_column_x(j,k), "y":nth_column_y(j,k), "z":nth_row_z(m)},ignore_index=True)
                else:
                    df = df.append({"RecoChannelID":toRecoChannelID(toAbsChannelID(GeoChannelID(i,j,k,m))),"x":-1*nth_column_x(j,k), "y":nth_column_y(j,k), "z":nth_row_z(m)},ignore_index=True)

df = df.sort_values(by="RecoChannelID").reset_index(drop=True)
df.to_csv("~/path/Geant4/Expirimental_Data/Scanner_Commissioning/Geometry/TPPT_Scanner_map.csv", header=False, index=False)

#Debug plots

# # full 3D display of coordinates as points
# cryst_x = []
# cryst_y = []
# cryst_z = []
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# for i in range(2):#2
#     for j in range(12):#12
#         for k in range(8):#8
#             for m in range(32):#32
#                 if i == 0:
#                     cryst_x.append(nth_column_x(j,k))
#                     cryst_y.append(nth_column_y(j,k))
#                     cryst_z.append(nth_row_z(m))
#                 else:
#                     cryst_x.append(-1*nth_column_x(j,k))
#                     cryst_y.append(nth_column_y(j,k))
#                     cryst_z.append(nth_row_z(m))
#
# ax.scatter(cryst_x, cryst_y, cryst_z, s=10)
# plt.show()

# #Testing if GeoChannelID works visually, also too cumbersome
# cryst_x = []
# cryst_y = []
# cryst_z = []
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# for i in range(1,2):#2
#     for j in range(1):#12
#         for k in range(8):#8
#             for m in range(8):#32
#                 if i == 0:
#                     cryst_x.append(nth_column_x(j,k))
#                     cryst_y.append(nth_column_y(j,k))
#                     cryst_z.append(nth_row_z(m))
#                     ax.text(nth_column_x(j,k), nth_column_y(j,k), nth_row_z(m), str(GeoChannelID(i,j,k,m)), fontsize="small")
#                 else:
#                     cryst_x.append(-1*nth_column_x(j,k))
#                     cryst_y.append(nth_column_y(j,k))
#                     cryst_z.append(nth_row_z(m))
#                     ax.text(-1*nth_column_x(j,k), nth_column_y(j,k), nth_row_z(m), str(GeoChannelID(i,j,k,m)), fontsize="small")
#
# ax.set_xlim(-170,170)
# ax.set_ylim(-170,170)
# ax.set_zlim(0,100)
# plt.show()

# # full 3D display of coordinates marked with indicies. TOO CUMBERSOME DO NOT RUN
# cryst_x = []
# cryst_y = []
# cryst_z = []
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# for i in range(1,2):#2
#     for j in range(6):#12
#         for k in range(8):#8
#             for m in range(32):#32
#                 if i == 0:
#                     cryst_x.append(nth_column_x(j,k))
#                     cryst_y.append(nth_column_y(j,k))
#                     cryst_z.append(nth_row_z(m))
#                     ax.text(nth_column_x(j,k), nth_column_y(j,k), nth_row_z(m), "({}, {}, {}, {})".format(i,j,k,m), fontsize="small")
#                 else:
#                     cryst_x.append(-1*nth_column_x(j,k))
#                     cryst_y.append(nth_column_y(j,k))
#                     cryst_z.append(nth_row_z(m))
#                     ax.text(-1*nth_column_x(j,k), nth_column_y(j,k), nth_row_z(m), "({}, {}, {}, {})".format(i,j,k,m), fontsize="small")
#
# ax.set_xlim(-170,170)
# ax.set_ylim(-170,170)
# ax.set_zlim(0,100)
# plt.show()

#Plotting the Geometric Slave/ChipIDs, to match up with back of scanner
fig = plt.figure()
fig1 = plt.figure()
ax = fig.add_subplot(111)
ax1 = fig1.add_subplot(111)
for i in range(12):
    for j in range(4):
        slave_chipL = module_dictL.get(str(i)+","+str(j))
        slave_chipR = module_dictR.get(str(i)+","+str(j))
        slaveL = slave_chipL // 100
        chipL = slave_chipL % 100
        slaveR = slave_chipR // 100
        chipR = slave_chipR % 100
        ax.text(i+0.05, 3-j + 0.5, "Slave: "+str(slaveL)+" Chip: "+str(chipL))
        ax1.text(i+0.05, 3-j + 0.5, "Slave: "+str(slaveR)+" Chip: "+str(chipR))
ax.grid()
ax1.grid()
ax.set_xlabel("x", fontsize=20)
ax.set_ylabel("z", fontsize =20)
ax1.set_xlabel("x", fontsize=20)
ax1.set_ylabel("z", fontsize =20)
ax.set_xticks(np.linspace(0,12,13))
ax.set_yticks(np.linspace(0,4,5))
ax1.set_xticks(np.linspace(0,12,13))
ax1.set_yticks(np.linspace(0,4,5))
ax.set_title("Left Crescent Layout")
ax1.set_title("Right Crescent Layout")
# ax.set_xlim(0,12)
# ax.set_ylim(0,4)
# ax1.set_xlim(0,12)
# ax1.set_ylim(0,4)
plt.show()

#Reading in the geometry file and showing the XZ projection of a particular array
def PlotChipID(slaveID, chipID):
    fig = plt.figure()
    RecoChannelID = 1024*slaveID + 64*chipID
    GeoChannelID = 10000*slaveID + 100*chipID
    geo_df = pd.read_csv("~/path/Geant4/Expirimental_Data/Scanner_Commissioning/Geometry/TPPT_Scanner_map.csv", sep=",", names=["ChannelID", "x", "y", "z"])
    for i in range(RecoChannelID, RecoChannelID+64):
        slice = geo_df.loc[geo_df["ChannelID"] == i].reset_index()
        plt.text(slice.y[0], slice.z[0], str(int(toGeoChannelIDfromRecoChannelID(slice.ChannelID[0]))), size="small")
    geo_df["GeoChannelID"] = geo_df["ChannelID"].apply(toGeoChannelIDfromRecoChannelID)
    ymin = geo_df.loc[geo_df["GeoChannelID"] == GeoChannelID].reset_index().y[0]
    ymax = geo_df.loc[geo_df["GeoChannelID"] == GeoChannelID+63].reset_index().y[0]
    zmax = geo_df.loc[geo_df["GeoChannelID"] == GeoChannelID].reset_index().z[0]
    zmin = geo_df.loc[geo_df["GeoChannelID"] == GeoChannelID+63].reset_index().z[0]

    #depending on what side of the scanner the chip is on, x can either increase or decrease with channel ID, so the maxes and mins are needed
    plt.xlim(min(ymin, ymax)-5, max(ymin, ymax)+5)
    plt.ylim(min(zmin, zmax) -5, max(zmin, zmax)+5)
    plt.title("Geometric Channel IDs \n"+"Slave: "+str(slaveID)+" ChipID: "+str(chipID))
    plt.xlabel("y", fontsize=20)
    plt.ylabel("z", fontsize =20)
    plt.show()


PlotChipID(1,0)
