import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import scipy as sp
from iminuit import cost, Minuit
import iminuit as im
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

def Exp2(x,A,tau, B, tau2):
    e = 2.718281828
    y = A*e**(-x/tau) + B*e**(-x/tau2)
    return y

def Exp(x,A,tau):
    e = 2.718281828
    y = A*e**(-x/tau)
    return y

def Exp3(x,A,tau1, B, tau2, C, tau3):
    e = 2.718281828
    y = A*e**(-x/tau1) + B*e**(-x/tau2) + C*e**(-x/tau3)
    return y

def ExpC(x,A, tau1, B, tau2, C):
    e = 2.718281828
    y = A*e**(-x/tau1) + B*e**(-x/tau2) + C
    return y

def ExpT(x,A,B,C):
    e = 2.718281828
    y = A*e**(-x/27.411) + B*e**(-x/173.123) + C
    return y

def ExpT2L(x,A,B,C):
    e = 2.718281828
    y = A*e**(-(x-6.95)/1765.123) + B*e**(-(x-6.95)/865.617) + C
    return y

def ExpT4(x,A,B,C,D,E):
    e = 2.718281828
    y = A*e**(-x/27.411) + B*e**(-x/173.123) + C*e**(-x/1765.859) + D*e**(-x/865.617) + E
    return y

def ExpOHGODWHY(x,C10,C11,N13,O15):
    e = 2.718281828
    y = C10*(e**(-(x)/27.411)+e**(-(x-4.1)/27.411) + e**(-(x-6.4)/27.411))+ C11*(e**(-(x)/1764.13)+e**(-(x-4.1)/1764.13) + e**(-(x-6.4)/1764.13)) +N13*(e**(-(x)/862.15)+e**(-(x-4.1)/862.15) + e**(-(x-6.4)/862.15))+O15*(e**(-(x-0.55)/176)+e**(-(x-4.65)/176) + e**(-(x-6.95)/176)) + 143.2452
    return y

def ExpBack(x,C10,C11,N13,O15):
    e = 2.718281828
    y = C10*e**(-(x-0.55)/27.411)+ C11*e**(-(x-0.55)/1764.13) +N13*e**(-(x-0.55)/862.15)+O15*e**(-(x-0.55)/176) + Y
    return y

def Constant(x,A):
    return A+(10**-7)*x

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


# file_name = "HorizontalSource_AlignedWithLowerArray_MaxSeparation_10min_coinc"
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/FLASH_Prep_Runs/"
file_name = "Phantom6_coinc"
dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/"

#For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
df = pd.read_csv(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
df.columns = ["TimeL", "ChargeL", "ChannelIDL", "TimeR", "ChargeR", "ChannelIDR"]

df["GeoChannelIDL"] = df["ChannelIDL"].apply(toGeoChannelID)
df["GeoChannelIDR"] = df["ChannelIDR"].apply(toGeoChannelID)
df["TimeL"] = df["TimeL"] / 1000000000000
df["TimeR"] = df["TimeR"] / 1000000000000
#djusting the time so the start time of the first spill is at 0:
df["TimeL"] = df["TimeL"] - 44.66

#Comparing Hollistic Data before 2 minutes and after

# #First we must determine when the first spill ocurred
# values,bins,params = plt.hist(df["TimeL"], bins=400)
# plt.xlim(-50,175)
# plt.yscale("log")
# plt.ylim(1,20000)
# plt.xlabel("Time [s]")
# plt.ylabel("No. Coincidences per 500 ms")
# plt.subplots_adjust(left=0.11, right=0.89, top=0.98,bottom=0.13)
# # plt.title("Number of Coincidences detected as a function of time \n Ge-68 Source Data Run 3")
# plt.show()

#Then we cut to right before the first spill occurred and we see that 156 seconds of data are recorded from when the first spill starts
# df = df[df["TimeL"] >= 44]
# #after the third spill
# #51.6 seconds after the start of the run
# decay = df[df["TimeL"] >= 6.95]


# #Now we can look at the left and right charge spectra
# plt.hist(df["ChargeR"], bins=2000)
# plt.xlabel("Charge Units")
# plt.ylabel("Counts")
# plt.xlim(-2,55)
# plt.title("Charge Values of Right Array, All Detectors")
# plt.show()



#Total Process for fitting the data to 4 exponentials
num_bins = 200
#fitting the background to a constant and then subtracting it from the data for scaling purposes
values,bins,params = plt.hist(df["TimeL"], bins=num_bins, fill=False, ec="C0")
plt.clf()
values = np.array(values)
bins = np.array(bins)
bin_centers = 0.5*(bins[1:] + bins[:-1])

values_const = values[(bin_centers <= -1)]
bin_centers_const = bin_centers[(bin_centers <= -1)]
c = cost.LeastSquares(bin_centers_const, values_const,np.sqrt(values_const), Constant)
fitter = Minuit(c,A=100)
fitter.limits = [(0,None)]
fitter.migrad()
fitter.hesse()
y_shift = fitter.values["A"]
print("Constant Background Fit")
print(fitter.values)
print(fitter.errors)
print(fitter.fval /(len(bins) - 2) )
print("Minuit isn't robust enough for you to set a value, put it in the big fitting function")

# #Now we can look at the counts right around the beam spills. The last spill ends around 52 seconds
values,bins,parameters = plt.hist(df["TimeL"], bins=num_bins)
values = np.array(values)
bins= np.array(bins)
bin_centers = 0.5*(bins[1:] + bins[:-1])
values = values[bin_centers >= 8]
bin_centers = bin_centers[bin_centers >= 8]

c = cost.LeastSquares(bin_centers, values,np.sqrt(values), ExpOHGODWHY)
fitter = Minuit(c,C10=1000,C11=1000,N13=1000, O15=1000)
fitter.limits = [(0,None),(0,None),(0,None), (0,None)]
fitter.migrad()
fitter.hesse()
print("4 exponential 3 spill Fit")
print(fitter.values)
print(fitter.errors)
print(fitter.fval / (len(values) - 5))

x_f = np.linspace(8,156,400)
plt.plot(x_f, ExpOHGODWHY(x_f, *fitter.values), color='k')
plt.xlabel("Time [s]")
plt.ylabel("No. Coincidences per s")
plt.ylim(0,2000)
t_y = generateTextY(5,max(values))
t_x = 60
#plt.title("Number of Coincidences during the spills: \n fixed decay constant 4 exponential 3 spill plus constant fit")
# plt.text(0.7*t_x,t_y[0], "Fitting Function: f(t) = $A*e^{-t/27.4} + B*e^{-t/173.1} + C$")
# plt.text(t_x,t_y[1], "C10 Decay Amplitude: "+format(fitter.values["A"], "6.2f"))
# plt.text(t_x,t_y[2], "O15 Decay Amplitude: "+format(fitter.values["B"], "6.2f"))
# plt.text(t_x,t_y[3], "Y-Constant: "+format(fitter.values["C"], "6.2f"))
# plt.text(t_x,t_y[4], "Least Squares Residue: "+format(fitter.fval, "6.2f"))
plt.show()

# # #Now we can look at the counts right around the beam spills. The last spill ends around 52 seconds
# num_bins = 44
# values,bins,parameters = plt.hist(decay["TimeL"], bins=num_bins)
# bin_centers = 0.5*(bins[1:] + bins[:-1])
# c = cost.LeastSquares(bin_centers, values,np.sqrt(values), ExpOHGODWHY)
# fitter = Minuit(c,C10=400,C11=195,N13=4, O15=850,Y=100)
# fitter.limits = [(100,510),(100,500),(3.85,6.75), (775,None), (50,None)]
# fitter.migrad()
# x_f = np.linspace(6.95,156,400)
#
# plt.plot(x_f, ExpOHGODWHY(x_f, *fitter.values), color='k')
# plt.xlabel("Time in s")
# # plt.xlim(6.95,156)
# plt.ylabel("Number of Coincidences")
# t_y = generateTextY(8,max(values))
# t_x = 60
# # plt.ylim(0,1.15*max(values))
# plt.title("Number of Coincidences after the spills")#: \n fixed decay constant 2 exponential plus constant fit")
# # plt.text(0.3*t_x,t_y[0], "Fitting Function: f(t) = $A*e^{-t/27.4} + B*e^{-t/173.1}$ \n" +"$+ C*e**(-x/1765.859) + D*e**(-x/865.617) + E$")
# # plt.text(t_x,t_y[1], "C11 Decay Amplitude: "+format(fitter.values["A"], "8.4f"))
# # plt.text(t_x,t_y[2], "O15 Decay Amplitude: "+format(fitter.values["B"], "8.4f"))
# # plt.text(t_x,t_y[3], "Amplitude Ratio: "+format((fitter.values["B"] / fitter.values["A"]), "8.4f"))
# # plt.text(t_x,t_y[4], "Y Constant: "+format(fitter.values["C"], "8.4f"))
# plt.text(t_x,t_y[1], "Chi Square: "+format((fitter.fval / (num_bins-5)), "6.2f"))
# # plt.text(t_x,t_y[4], "C11 Decay Amplitude: "+format(fitter.values["A"], "6.2f"))
# # # plt.text(t_x,t_y[5], "N10 Decay Amplitude: "+format(fitter.values["D"], "6.2f"))
# print(fitter.values)
# print("O15", fitter.values["O15"]/fitter.values["N13"])
# print("C10", fitter.values["C10"]/fitter.values["N13"])
# print("C11", fitter.values["C11"]/fitter.values["N13"])
# plt.show()

# # #Now we can look at the counts right around the beam spills. The last spill ends around 52 seconds
# num_bins = 40
# # fit_data = decay.loc[decay["TimeL"] >= 106.95]["TimeL"]
# values,bins,parameters = plt.hist(decay["TimeL"], bins=num_bins)
# # plt.show()
# bin_centers = 0.5*(bins[1:] + bins[:-1])
# #x1_p, x1_co = curve_fit(ExpC, bin_centers, values, p0=[300,120,300,20,100], bounds = [[0,0,0,0,0], [1000000,1000000,1000000,1000000,100000]])
# c = cost.LeastSquares(bin_centers, values,np.sqrt(values), ExpOHGODWHY)
# fitter = Minuit(c,A=1000,B=1000,C=1000, D=1000, E=1000,F=1000,G=1000,H=1000,I=1000,J=1000,K=1000,L=1000,M=100)
# fitter.limits = [(0,200),(0,200),(0,200),(0,200),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None)]
# fitter.migrad()
# x_f = np.linspace(6.95,156,400)
# # plt.clf()
# # values,bins,parameters = plt.hist(decay["TimeL"], bins=num_bins)
# plt.plot(x_f, ExpOHGODWHY(x_f, *fitter.values), color='k')
# plt.xlabel("Time in s")
# # plt.xlim(6.95,156)
# plt.ylabel("Number of Coincidences")
# t_y = generateTextY(8,max(values))
# t_x = 60
# # plt.ylim(0,1.15*max(values))
# plt.title("Number of Coincidences after the spills")#: \n fixed decay constant 2 exponential plus constant fit")
# # plt.text(0.3*t_x,t_y[0], "Fitting Function: f(t) = $A*e^{-t/27.4} + B*e^{-t/173.1}$ \n" +"$+ C*e**(-x/1765.859) + D*e**(-x/865.617) + E$")
# # plt.text(t_x,t_y[2], "C10 Decay Amplitude: "+format(fitter.values["A"], "6.2f"))
# # plt.text(t_x,t_y[3], "O15 Decay Amplitude: "+format(fitter.values["B"], "6.2f"))
# # plt.text(t_x,t_y[4], "C11 Decay Amplitude: "+format(fitter.values["A"], "6.2f"))
# # # plt.text(t_x,t_y[5], "N10 Decay Amplitude: "+format(fitter.values["D"], "6.2f"))
# # plt.text(t_x,t_y[6], "Y Constant: "+format(fitter.values["C"], "6.2f"))
# # plt.text(t_x,t_y[7], "Least Squares Residue: "+format(fitter.fval, "6.2f"))
# print(fitter.values)
# plt.show()

# # First spill 44.55 to 45.2 - 44.66 = -0.075 to 0.15, 2250 bins
# # Second spill 48.7 to 49.3 - 44.66 = 4.04 to 4.25, 2100 bins
# # Thrid spill 51.0 to 51.6 - 44.66 = 6.35 to 6.55, 2000 bins
# #All 3
# spill1 = df[df["TimeL"] >= -0.075]
# spill1 = spill1[spill1["TimeL"] <= 7]
# plt.hist(spill1["TimeL"], bins=7074)
# # plt.title("All three spills: 1 ms per bin")
# plt.xlabel("Time [s]")
# plt.ylabel("No. Coinidences per ms")
# plt.ylim(0,250)
# plt.show()

# Before = df[df["TimeL"] < 171.6]
# After = df[df["TimeL"] >= 171.6]
#
# fig = plt.figure()
# ax = fig.add_subplot(121)
# ax1 = fig.add_subplot(122)
# ax.hist(Before["ChargeR"], bins=2000)
# ax1.hist(After["ChargeR"], bins=1100)
# ax.set_title("Charge Spectra during \n the first two minutes after beam spills")
# ax1.set_title("Charge Spectra after \n the first two minutes after beam spills")
# ax1.set_xlabel("Charge Units")
# ax.set_xlabel("Charge Units")
# ax.set_ylabel("Number of Coincidences")
# ax.set_xlim(-2,55)
# ax1.set_xlim(-2,55)
# plt.show()


# #Analyzing the group data
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/"
# file_name = "Phantom6_group"
# df = pd.read_csv(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4))
# df.columns = ["Time", "Charge", "ChannelID"]
# df["GeoChannelID"] = df["ChannelID"].apply(toGeoChannelID)
# df["Time"] = df["Time"] / 1000000000000
#
# #looking at any depositions during the time of no Coincidences
# spill1 = df[df["Time"] >= 44.5]
# spill1 = spill1[spill1["Time"] <= 45.2]
# plt.hist(spill1["Time"],bins=8000)
# plt.title("First Spill: \n 100 $\mu$s per bin")
# plt.xlabel("Time in s")
# plt.ylabel("Number of Coinidences")
# plt.show()

# #reading in Andrey's graphs
# C10 = pd.read_csv(dir+"Andrey_Plots/C10Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# # C10["Values"] = C10["Values"] * 10**6
#
# C11 = pd.read_csv(dir+"Andrey_Plots/C11Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# # C11["Values"] = C11["Values"] * 10**6
#
# N13 = pd.read_csv(dir+"Andrey_Plots/N13Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# # N13["Values"] = N13["Values"] * 10**6
#
# O15 = pd.read_csv(dir+"Andrey_Plots/O15Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# # O15["Values"] = O15["Values"] * 10**6
#
# Tot = pd.read_csv(dir+"Andrey_Plots/TotalValues.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# # Tot["Values"] = Tot["Values"] * 10**6
#
# bins = np.linspace(0,249.9,2500)
# plt.plot(bins,C10.Values, linestyle="--", color="b", label="C10")
# plt.plot(bins,C11.Values, linestyle="dotted", color="mediumblue", label="C11")
# plt.plot(bins,N13.Values, linestyle="dotted", color="cyan", label="N13")
# plt.plot(bins,O15.Values, linestyle="dotted", color="m", label="O15")
# plt.plot(bins,Tot.Values, linestyle="-", color="r", label="Sum")
# plt.ylabel("Activity per primary proton, $s^{-1}$")
# plt.xlabel("Time from first spill, s")
# plt.legend()
# plt.xlim(0,250)
# plt.yticks(ticks=[5*i*10**-6 for i in range(11)],labels=[str(5*i) for i in range(11)])
# plt.text(0,50.5*(10**-6), "x $10^{-6}$", size='large')
# plt.ylim(0,50*(10**-6))
# plt.show()
