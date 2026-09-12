import matplotlib.pyplot as plt
import numpy as np
from iminuit import cost, Minuit
import iminuit as im
import pandas as pd

def Expo(x,A,lamda, C):
    y = A*(2.718281828)**(-1*lamda*x) + C
    return y

file_name = "MDA-H2O1-15min-210mm-NoCu_coinc"
dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/FLASH_Prep_Runs/"

#For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
df = pd.read_csv(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
df.columns = ["TimeL", "ChargeL", "ChannelIDL", "TimeR", "ChargeR", "ChannelIDR"]

df["TimeL"] = df["TimeL"] / 1000000000000
df["TimeR"] = df["TimeR"] / 1000000000000

# #determining cutoff point
# plt.hist(df["TimeL"], bins=np.linspace(10,15,500))
# plt.ylim(1,10**5)
# plt.yscale("log")
# plt.show()

cutoff = 13.5
fit_data = df.loc[df["TimeL"] >= cutoff]
fit_data["TimeL"] = fit_data["TimeL"] - 13.5

values,bins,params = plt.hist(fit_data["TimeL"], bins=np.linspace(0,886.5,887))
plt.show()
