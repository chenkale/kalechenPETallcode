#!/usr/bin/env python
# coding: utf-8

# Imports, Indices, helper functions, constants

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import scipy as sp
from iminuit import cost, Minuit
import iminuit as im
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
# plt.style.use('ChannelPairStyles.mplstyle')

ln2 = 0.6391
constant = 42.447360677913636
e = 2.718281828

def generateTextY(num,max,subplots=False):
    if max > 50 and subplots==False:
        y = [max - 0.04*(i+1)*max for i in range(0,num+1)]
    elif max <= 50 and subplots == False:
        y = [max - 0.06*(i+1)*max for i in range(0,num+1)]
    elif subplots == True:
        y = [0.96*max - 0.08*(i+1)*max for i in range(0,num+1)]
    return y

def halferror(x,xerr):
    return ln2*x * xerr/x

def toGeo(x):
    #converts PETSys ID to geometric ID
    y = 8*indices.get(x)[0] + indices.get(x)[1]
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


# Spill times/run lengths need to be adjusted

# In[2]:


spillTimes = [[52.60406546037813, 52.70576546037813],
 [7.658639091560348, 7.760139091560347],
 [21.685239651808505, 21.786839651808506],
 [62.611831181721065, 62.71383118172106],
 [10.560941430694248, 10.662941430694248],
 [68.50573713055009, 68.60743713055008],
 [52.370972189816804, 52.4725721898168],
 [13.310021040987447, 13.411821040987448],
 [69.17295778706057, 69.27475778706057],
 [23.0682897979999, 23.1698897979999],
 [39.9248759449842, 40.0264759449842], #nickel
 [76.49762062100974, 76.59942062100974], #copper
 [100.57416844934949, 100.67546844934948]]

runLengths = [180,180,180,1200,1200,900,900,900,300,1200,600,600,900] #nickel index 10, copper index 11

def chiSq(o,e):
    return (o-e)**2/e

def halferror(x,xerr):
    return ln2*x * xerr/x


# ## fit to background
# 
# Pass in file name/dir

# In[4]:


file_name = 'MDA-Ni-10min-210mm-NoCu_coinc'
dir = "Flash_Therapy/PET_3-5-23/"

#For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
df = pd.read_csv(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
df.columns = ["TimeL", "ChargeL", "ChannelIDL", "TimeR", "ChargeR", "ChannelIDR"]

df["GeoChannelIDL"] = df["ChannelIDL"].apply(toGeoChannelID)
df["GeoChannelIDR"] = df["ChannelIDR"].apply(toGeoChannelID)
df["TimeL"] = df["TimeL"] / 1000000000000
df["TimeR"] = df["TimeR"] / 1000000000000
#djusting the time so the start time of the first spill is at 0:
df["TimeL"] = df["TimeL"] - spillTimes[10][0] #44.66

#Total Process for fitting the data to 4 exponentials
num_bins = 100
binwidth = runLengths[10]/num_bins
#fitting the background to a constant and then subtracting it from the data for scaling purposes
values,bins,params = plt.hist(df["TimeL"], bins=num_bins, fill=False, ec="C0")
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
# bounds = list(ax.get_xlim())
# leftbound = bounds[0]
# rightbound = bounds[0]
plt.xlim(-30,200)
plt.ylim(0,4000)
x_f = np.linspace(-5,0,len(values_const))
plt.plot(x_f, Constant(x_f, *fitter.values))
print("Minuit isn't robust enough for you to set a value, put it in the big fitting function")


# Creation of the following box should be programmed in a way to take the next box with the specific isotopes. Likely create a class that can return a single in the form of
# def Cu58(t,A):
#     tau_Cu58 = 3.204/ln2
#     return A*e**(-t/tau_Cu58)
# 
# Probably have the tau or half life be an initialization parameter and the actual function above be a subroutine.

# In[6]:


def single(t,A,tau1):
    return A*e**(-t/tau1) + constant

def Cu58(t,A):
    tau_Cu58 = 3.204/ln2
    return A*e**(-t/tau_Cu58)

def double(t,A,tau1,B,tau2):
    return A*e**(-t/tau1) + B*e**(-t/tau2)  + constant #C*e**(-t/tau3)

def double_fit2constant(t,A,tau1,B,tau2,C):
    return A*e**(-t/tau1) + B*e**(-t/tau2) + C

def triple(t,A,tau1,B,tau2,C,tau3):
    return A*e**(-t/tau1) + B*e**(-t/tau2)  + C*e**(-t/tau3) + constant


# ## Helper fitting functions
# ## https://journals.aps.org/pr/abstract/10.1103/PhysRev.54.1021

# In[88]:


tau_Cu58 = 3.204/ln2
tau_Cu59 = 89.5/ln2
tau_Cu60 = 1422/ln2
tau_Cu62 = 9.672*60/ln2
tau_Zn60 = 2.38*60/ln2
tau_Zn61 = 89.1/ln2
tau_Zn63 = 38.47*60/ln2
tau_Co54m = 88.8/ln2

def double_Cu58_Cu62(t,A,B,C):
    return A*e**(-t/tau_Cu58) + B*e**(-t/tau_Cu62) + C

def triple_Cu58_Cu62_Cu59(t,A,B,C):
    return A*e**(-t/tau_Cu58) + B*e**(-t/tau_Cu62) + C*e**(-t/tau_Cu59)

def triple_Cu58_Cu62_Zn60(t,A,B,C):
    return A*e**(-t/tau_Cu58) + B*e**(-t/tau_Cu62) + C*e**(-t/tau_Zn60)

def quad_Cu58_Cu62_Cu59_Zn60(t,A,B,C,D):
    return A*e**(-t/tau_Cu58) + B*e**(-t/tau_Cu62) + C*e**(-t/tau_Cu59) + D*e**(-t/tau_Zn60)

def quad_Cu58_Cu60_Cu59_Cu62(t,A,B,C,D):
    return A*e**(-t/tau_Cu58) + B*e**(-t/tau_Cu60) + C*e**(-t/tau_Cu59) + D*e**(-t/tau_Cu62)

def all_Cu58_Cu62_Cu59_Cu60_Zn60_Zn61_Zn63_Co54m(t,A,B,C,D,E,F,G,H):
    return A*e**(-t/tau_Cu58) + B*e**(-t/tau_Cu62) + C*e**(-t/tau_Cu59) + D*e**(-t/tau_Cu60) + E*e**(-t/tau_Zn60) + F*e**(-t/tau_Zn61) + G*e**(-t/tau_Zn63) + H*e**(-t/tau_Co54m)

def Cu58(t,A):
    return A*e**(-t/tau_Cu58)

def Cu59(t,A):
    return A*e**(-t/tau_Cu59)

def Cu60(t,A):
    return A*e**(-t/tau_Cu60)

def Cu62(t,A):
    return A*e**(-t/tau_Cu62)

def Zn60(t,A):
    return A*e**(-t/tau_Zn60)

def Zn61(t,A):
    return A*e**(-t/tau_Zn61)

def Zn63(t,A):
    return A*e**(-t/tau_Zn63)

def Co54m(t,A):
    return A*e**(-t/tau_Co54m)


# In[98]:


###### fig,ax = plt.subplots(figsize = (16,9))

t = df["TimeL"][df["TimeL"] >= 1] - min(df["TimeL"][df["TimeL"] >= 1] )
new_num_bins = int((max(t) - min(t))/binwidth)
values,bins,parameters = plt.hist(t, bins=new_num_bins,label = 'PET Data',alpha = 0.8)
values = np.array(values)
bins= np.array(bins)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# values = values[bin_centers >= 8]
# bin_centers = bin_centers[bin_centers >= 8]

c = cost.LeastSquares(bin_centers, values,np.sqrt(values), triple_Cu58_Cu62_Cu59)
fitter = Minuit(c,A = 500,B = 250,C=68)#,C = 75, tau3 = 0.4)
fitter.limits = [(0,None),(0,None),(0,None)]#,(0,None),(0,None)]
fitter.migrad()
fitter.hesse()
print("Single Exponential Fit")
print(fitter.values)
print(fitter.errors)
redChiSq = fitter.fval / (len(values) - 3)
print(redChiSq)

x_f = np.linspace(2,max(t),len(values))
plt.plot(x_f,Cu58(x_f,fitter.values[0]),linewidth = 4,color = 'navy',linestyle = 'dotted',label = 'Cu58')
plt.plot(x_f,Cu62(x_f,fitter.values[1]),linewidth = 4,color = 'orange',linestyle = 'dashdot',label = 'Cu62')
plt.plot(x_f,Cu59(x_f,fitter.values[2]),linewidth = 4,color = 'green',linestyle = 'dashdot',label = 'Cu59')
#plt.plot(x_f,Cu62(x_f,fitter.values[3]),linewidth = 4,color = 'purple',linestyle = 'dashdot',label = 'Cu62')
plt.plot(x_f, triple_Cu58_Cu62_Cu59(x_f, *fitter.values),linewidth = 4,color="r",alpha=0.7,label = 'Sum')

plt.legend(ncol = 2,fontsize = 25)

plt.ylabel(r'PET Events [s$^{-1}$]')
plt.xlabel('Time [s]')
plt.rcParams['figure.figsize'] = [20, 9]

#plt.yscale('log')
plt.text(35,4.8*10**2,r'A = ' + str(np.round(fitter.values[0],1)) + ' ± ' + str(np.round(fitter.errors[0],1)),fontsize = 20)
plt.text(35,4.5*10**2,r'T$_{1/2,A}$ (assumed) = 3.204 s (T$_{1/2}$ of $^{58}$Cu)',fontsize = 20)
plt.text(35,4.2*10**2,r'B = ' + str(np.round(fitter.values[1],1)) + ' ± ' + str(np.round(fitter.errors[1],1)),fontsize = 20)
plt.text(35,3.9*10**2,r'T$_{1/2,B}$ (assumed) = 580.3 s (T$_{1/2}$ of $^{62}$Cu)',fontsize = 20)
plt.text(35,3.6*10**2,r'C = ' + str(np.round(fitter.values[2],1)) + ' ± ' + str(np.round(fitter.errors[2],1)),fontsize = 20)
plt.text(35,3.3*10**2,r'T$_{1/2,C}$ (assumed) = 89.5 s (T$_{1/2}$ of $^{59}$Cu)',fontsize = 20)
#plt.text(35,3*10**2,r'D = ' + str(np.round(fitter.values[3],1)) + ' ± ' + str(np.round(fitter.errors[3],1)),fontsize = 20)
#plt.text(35,2.7*10**2,r'T$_{1/2,D}$ (assumed) = 142.8 s (T$_{1/2}$ of $^{60}$Zn)',fontsize = 20)
plt.text(35,2.4*10**2,r'reduced $\chi^{2}$ = ' + str(np.round(redChiSq,3)),fontsize = 20)
plt.text(35,2*10**2,'Ni (21 cm)',fontsize = 35)
# plt.ylim(10,10**3)
plt.xlim(0,630)
plt.savefig('asdf.png')


# In[ ]:




