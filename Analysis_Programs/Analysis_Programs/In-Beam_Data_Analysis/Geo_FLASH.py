#Kyle Klein
#1-17-23
#Lays out the geometry for two horizontal array
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from tqdm import tqdm

#Everything is in mm
x_distance = int(input("What is the distance from the center of the source to the crystal faces in mm? "))

#The coordinate system is set up such that positive z is the direction the beam is traveling. FAcing this direction, positive y is up and positive x is to the left
#as a clarification, because the row is horizontal, what you would refer to as a row in the crystal array is positioned like a column and vice versa.

#arrays seperated by 3.2 mm per column (except between arrays of 64) and z=0 is set to be the middle of the arrays
def nth_column_z(ArrCol):
    #the -0.4 is to center the array around z=0
    if ArrCol < 8:
        z = - 1.6 + 3.2*(ArrCol-7) -0.4
        return z
    elif ArrCol >= 8 and ArrCol < 16:
        z = -1.6+ 4 + 3.2*(ArrCol-8) -0.4
        return z

#middle of the arrays is 1.6+3.2+3.2+3.2+1.6
def nth_row_y(ArrRow):
    return 3.2*ArrRow - 11.2

def toGeoChannelID(LR, ArrCol, ArrRow):
    if LR == 0:
        if ArrCol < 8:
            ChipID = 8
            new_ArrCol = ArrCol
        else:
            ChipID = 9
            new_ArrCol = ArrCol-8

        return 100*ChipID + 8*new_ArrCol + ArrRow
    else:
        if ArrCol < 8:
            ChipID = 6
            new_ArrCol = ArrCol
        else:
            ChipID = 7
            new_ArrCol = ArrCol-8
        return 100*ChipID + 8*new_ArrCol + (7-ArrRow)

def detectorSide(LR):
    global x_distance
    return -2*x_distance*LR + x_distance

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

def toRecoChannelID(AbsChannelID):
    slaveID = AbsChannelID // 4096
    chipID = (AbsChannelID - slaveID*4096) // 64
    channelID = AbsChannelID % 64
    RecoChannelID = 1024*slaveID + 64*chipID + channelID
    return RecoChannelID

df = pd.DataFrame(columns=["GeoChannelID","x", "y", "z", "nx", "ny", "nz"])
#both sides of the detector
for i in range(2):
    #16 rows
    for j in range(16):
        #8 columns
        for k in range(8):
            #[toRecoChannelID(toAbsChannelID(toGeoChannelID(i,j,k)))]
            slice = pd.DataFrame({"GeoChannelID": [toRecoChannelID(toAbsChannelID(toGeoChannelID(i,j,k)))], "x":[detectorSide(i)], "y":[nth_row_y(k)], "z":[nth_column_z(j)], "nx":[-1*np.sign(detectorSide(i))], "ny":[0], "nz":[0]})
            df = pd.concat([df,slice], ignore_index=True)

df = df.sort_values(by="GeoChannelID").reset_index(drop=True)
df.to_csv("~/path/Geant4/Expirimental_Data/Scanner_Commissioning/Geometry/Alex_FLASH_Scanner_map_3-5-23_{}mm.csv".format(x_distance), header=False, index=False)
