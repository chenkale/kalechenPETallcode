#converting Andrey plots to python for the FLASH paper
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from iminuit import cost, Minuit
import scipy.stats as sp
plt.style.use('ChannelPairStyles.mplstyle')

def ReducedChiSquare(data,model):
    chi = 0
    for i in range(len(data)):
        chi += (data[i] - model[i])**2 / data[i]
    return chi/len(data)

def ExpandedChiSquare(data,model):
    differences = []
    chis = []
    for i in range(len(data)):
        differences.append(data[i] - model[i])
        chis.append(((data[i] - model[i])**2) / data[i])

    plt.hist(differences, bins=50)
    plt.xlabel("Data - Model")
    plt.ylabel("Differences per bin")
    plt.title("Differences between data and Model")
    plt.show()
    plt.hist(chis, bins=50,density=True)
    x_f = np.linspace(0,50,400)
    plt.plot(x_f, sp.chi2.pdf(x_f,1), color='k')
    plt.title("$\chi^2$ Values compared to $\chi^2$ PDF for {} dof".format(len(data)-1))
    plt.xlabel("Normalized $\chi^2$")
    plt.ylabel("$\chi^2$ per bin")
    plt.show()
    plt.plot()

    return 0

def Constant(x,A):
    return A+(10**-7)*x

PES_dict = {"C11":["navy", "dotted"], "C10":["crimson", "dashed"], "N13":["g", (5,(10,3))], "O15":["m", "dashdot"], "Sum":["r", "solid"]}

#adding expirimenal data to histogram to compare
file_name = "Phantom6_coinc"
dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/"

# #For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
# df = pd.read_csv(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
# df.columns = ["TimeL", "ChargeL", "ChannelIDL", "TimeR", "ChargeR", "ChannelIDR"]
#
# # df["GeoChannelIDL"] = df["ChannelIDL"].apply(toGeoChannelID)
# # df["GeoChannelIDR"] = df["ChannelIDR"].apply(toGeoChannelID)
# df["TimeL"] = df["TimeL"] / 1000000000000
# df["TimeR"] = df["TimeR"] / 1000000000000
# #djusting the time so the start time of the first spill is at 0:
# df["TimeL"] = df["TimeL"] - 44.66
#
# num_bins = 200
# fig,ax = plt.subplots(figsize=(16,9))
# values,bins,params = ax.hist(df["TimeL"], bins=num_bins, color="C0",alpha=0.7, label="PET Data")
# # ax_twin.set_xlim(-50,175)
# ax.set_ylim(0,1500)
# ax.set_ylabel("PET Events [$s^{-1}$]")
#
# #fitting the background to a constant and then subtracting it from the data for scaling purposes
# values = np.array(values)
# bins = np.array(bins)
# bin_centers = 0.5*(bins[1:] + bins[:-1])
#
# values_const = values[(bin_centers <= -1)]
# bin_centers_const = bin_centers[(bin_centers <= -1)]
# c = cost.LeastSquares(bin_centers_const, values_const,np.sqrt(values_const), Constant)
# fitter = Minuit(c,A=100)
# fitter.limits = [(0,None)]
# fitter.migrad()
# fitter.hesse()
# print("Constant Background Fit")
# print(fitter.values)
# print(fitter.errors)
# print(fitter.covariance)
#
#
# #PES population vs time with expirimental data added
# bins = np.linspace(-0.1,249.9,2501)
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/Andrey_Plots/75.8MeVGraphs/ActivityVsTime/"
# C10 = pd.read_csv(dir+"C10Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# C10["Values"] = C10["Values"] * 10**6
# C10["Bins"] = bins
# C10 = C10.loc[C10["Bins"] <= 156]
#
# C11 = pd.read_csv(dir+"C11Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# C11["Values"] = C11["Values"] * 10**6
# C11["Bins"] = bins
# C11 = C11.loc[C11["Bins"] <= 156]
#
# N13 = pd.read_csv(dir+"N13Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# N13["Values"] = N13["Values"] * 10**6
# N13["Bins"] = bins
# N13 = N13.loc[N13["Bins"] <= 156]
#
# O15 = pd.read_csv(dir+"O15Values.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# O15["Values"] = O15["Values"] * 10**6
# O15["Bins"] = bins
# O15 = O15.loc[O15["Bins"] <= 156]
#
# Tot = pd.read_csv(dir+"SumValues.csv", sep=',', header=None, low_memory=False,usecols=[0], names=["Values"])
# Tot["Values"] = Tot["Values"] * 10**6
# Tot["Bins"] = bins
# Tot = Tot.loc[Tot["Bins"] <= 156]
#
# values_norm = values - fitter.values["A"]
# sim_values_norm = Tot.loc[(Tot["Bins"] >= 135) & (Tot["Bins"] <= 155)].reset_index(drop=True)["Values"]
# values_norm = values_norm[(bin_centers >= 135) & (bin_centers <= 155)]
# norm_factor = values_norm.mean() / sim_values_norm.mean()
# Tot["Values"] = Tot["Values"]*norm_factor + fitter.values["A"]
# C10["Values"] = C10["Values"]*norm_factor + fitter.values["A"]
# C11["Values"] = C11["Values"]*norm_factor + fitter.values["A"]
# N13["Values"] = N13["Values"]*norm_factor + fitter.values["A"]
# O15["Values"] = O15["Values"]*norm_factor + fitter.values["A"]
# const = pd.DataFrame({"Bins": np.linspace(-44.6,-0.2,445), "Values":445*[fitter.values["A"]]})
# C10 = pd.concat([C10, const], ignore_index=True).sort_values(by="Bins")
# C11 = pd.concat([C11, const], ignore_index=True).sort_values(by="Bins")
# N13 = pd.concat([N13, const], ignore_index=True).sort_values(by="Bins")
# O15 = pd.concat([O15, const], ignore_index=True).sort_values(by="Bins")
# Tot = pd.concat([Tot, const], ignore_index=True).sort_values(by="Bins")
#
# ax_twin = ax.twinx()
# ax_twin.step(C10.Bins,C10.Values, linestyle=PES_dict.get("C10")[1], color=PES_dict.get("C10")[0], label="C10", alpha=0.7)
# ax_twin.step(C10.Bins,C11.Values, linestyle=PES_dict.get("C11")[1], color=PES_dict.get("C11")[0], label="C11", alpha=0.7)
# ax_twin.step(C10.Bins,N13.Values, linestyle=PES_dict.get("N13")[1], color=PES_dict.get("N13")[0], label="N13", alpha=0.7)
# ax_twin.step(C10.Bins,O15.Values, linestyle=PES_dict.get("O15")[1], color=PES_dict.get("O15")[0], label="O15", alpha=0.7)
# ax_twin.step(C10.Bins,Tot.Values, linestyle=PES_dict.get("Sum")[1], color=PES_dict.get("Sum")[0], label="Sum", alpha=0.7)
#
# # ax_twin.fill_between(C10.Bins,C10.Values, step="pre", alpha=0.3, color='b')
# # ax_twin.fill_between(C10.Bins,C11.Values, step="pre", alpha=0.3, color='mediumblue')
# # ax_twin.fill_between(C10.Bins,N13.Values, step="pre", alpha=0.3, color='g')
# # ax_twin.fill_between(C10.Bins,O15.Values, step="pre", alpha=0.3, color='m')
# # ax_twin.fill_between(C10.Bins,Tot.Values, step="pre", alpha=0.3, color='r')
# # ax_twin.fill_between(np.linspace(-44.6,-0.1,443),443*[fitter.values["A"]], step="pre", alpha=0.3,color='b')
# # ax_twin.fill_between(np.linspace(-44.6,-0.1,443),443*[fitter.values["A"]], step="pre", alpha=0.3,color='mediumblue')
# # ax_twin.fill_between(np.linspace(-44.6,-0.1,443),443*[fitter.values["A"]], step="pre", alpha=0.3,color='g')
# # ax_twin.fill_between(np.linspace(-44.6,-0.1,443),443*[fitter.values["A"]], step="pre", alpha=0.3,color='m')
# # ax_twin.fill_between(np.linspace(-44.6,-0.1,443),443*[fitter.values["A"]], step="pre", alpha=0.3,color='r')
#
# #ax_twin.hlines(y = fitter.values["A"], xmin = -44.66, xmax = 0, color='r', alpha=0.5)
# ax_twin.set_ylabel("Normalized Activity [Bq]")
# ax_twin.set_xlabel("Time [s]")
# lines, labels = ax.get_legend_handles_labels()
# lines2, labels2 = ax_twin.get_legend_handles_labels()
# ax_twin.legend(lines + lines2, labels + labels2, loc=0, ncol=2)
# ax_twin.set_xlim(-50,175)
# ax_twin.set_ylim(0,1500)
# ax.set_xlabel("Time [s]")
# plt.subplots_adjust(left=0.13, right=0.87, top=0.98,bottom=0.13)
# #ax_twin.set_title("Using the ratios of averages in 135-155 as norm factor")
# # plt.yticks([5*i for i in range(9)])
# # plt.tick_params(right = False, labelright = False)
#
# plt.show()

# #Dose and PES
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/Andrey_Plots/75.8MeVGraphs/DosePES/"
# Dose = pd.read_csv(dir+"Dose.csv", sep=',', header=None, low_memory=False, names=["Bin Centers","Values"])
# # C10["Values"] = C10["Values"] * 10**6
#
# PES = pd.read_csv(dir+"PES.csv", sep=',', header=None, low_memory=False,names=["Bin Centers","Values"])
# # C11["Values"] = C11["Values"] * 10**6
# plt.step(PES["Bin Centers"],PES.Values, color="r", label="PES")
# plt.step(Dose["Bin Centers"],Dose.Values, color="b", label="Dose")
# plt.ylabel("Linear density per primary proton [$mm^{-1}$]")
# plt.xlabel("Depth [mm]")
# plt.legend()
# plt.xlim(0,60)
# #plt.yticks(ticks=[5*i*10**-6 for i in range(11)],labels=[str(5*i) for i in range(11)])
# # plt.text(-50,40.5*(10**-6), "x $10^{-6}$")
# plt.ylim(0,1.2)
# plt.show()

# #PES for Reco
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/Andrey_Plots/75.8MeVGraphs/ActivityForPETrec/"
# PES = pd.read_csv(dir+"Activity-155s-1D.csv", sep=',', header=None, low_memory=False, names=["Bin Centers","Values"])
# PES["Values"] = PES["Values"]*10**6
# plt.step(PES["Bin Centers"],PES.Values, color="r")
# plt.ylabel("Total activity [$s^{-1}$ per $10^6$ protons]")
# plt.xlabel("Depth [mm]")
# plt.xlim(0,60)
# # plt.yticks(ticks=[2*i*10**-5 for i in range(6)],labels=[str(2*i*10**-2) for i in range(6)])
# # plt.text(0,0.11*10**(-3), "x $10^{-3}$")
# plt.ylim(0,120)
# plt.show()

# # #Geant4 and Expiriment verification
# #adding expirimenal data to histogram to compare
# file_name = "Phantom6_coinc"
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/"
#
# #For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
# df = pd.read_csv(dir+"{}.dat".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
# df.columns = ["TimeL", "ChargeL", "ChannelIDL", "TimeR", "ChargeR", "ChannelIDR"]
#
# # df["GeoChannelIDL"] = df["ChannelIDL"].apply(toGeoChannelID)
# # df["GeoChannelIDR"] = df["ChannelIDR"].apply(toGeoChannelID)
# df["TimeL"] = df["TimeL"] / 1000000000000
# df["TimeR"] = df["TimeR"] / 1000000000000
# #djusting the time so the start time of the first spill is at 0:
# df["TimeL"] = df["TimeL"] - 44.66
# values,bins,params = plt.hist(df["TimeL"], bins=200, alpha=0.5, color="C0", label="PET Data")
# bin_centers = 0.5*(bins[1:] + bins[:-1])
# Exp = pd.DataFrame({"Bin Centers":bin_centers, "Values":values})
# plt.close()
#
#
#
# dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/Andrey_Plots/75.8MeVGraphs/ExpFramewGeant/"
# Framework = pd.read_csv(dir+"Framework.csv", sep=',', header=None, low_memory=False,names=["Bin Centers","Values"])
# Framework["Values"] = Framework["Values"] * 10**6
# Framework["Bin Centers"] = Framework["Bin Centers"] - 44.66
#
# G4 = pd.read_csv(dir+"Geant4.csv", sep=',', header=None, low_memory=False,names=["Bin Centers","Values"])
# G4["Values"] = G4["Values"] * 10**6
# G4["Bin Centers"] = G4["Bin Centers"] - 44.66
#
# t_min = [135,120,120]
# t_max = [155,155,140]
#
# for i in range(len(t_min)):
#     normal_exp = Exp.loc[Exp["Bin Centers"] >= t_min[i]]
#     normal_exp = normal_exp.loc[normal_exp["Bin Centers"] <= t_max[i]]
#
#     normal_g4 = G4.loc[G4["Bin Centers"] >= t_min[i]]
#     normal_g4 = normal_g4.loc[normal_g4["Bin Centers"] <= t_max[i]]
#
#     normal_frame = Framework.loc[Framework["Bin Centers"] >= t_min[i]]
#     normal_frame = normal_frame.loc[normal_frame["Bin Centers"] <= t_max[i]]
#
#     factor = normal_exp["Values"] / normal_g4["Values"]
#     norm_factor = factor.mean()
#     G4["Values"] = G4["Values"]*norm_factor
#     factor = normal_exp["Values"] / normal_frame["Values"]
#     norm_factor = factor.mean()
#     Framework["Values"] = Framework["Values"]*norm_factor
#
#     chi_g4 = G4.loc[G4["Bin Centers"] >= 25]
#     chi_g475 = chi_g4.loc[chi_g4["Bin Centers"] <= 75]
#     chi_g4end = G4.loc[G4["Bin Centers"] >= 75]
#
#     chi_frame = Framework.loc[Framework["Bin Centers"] >= 25]
#     chi_frame75 = chi_frame.loc[chi_frame["Bin Centers"] <= 75]
#     chi_frameend = Framework.loc[Framework["Bin Centers"] >= 75]
#
#     chi_exp = Exp.loc[Exp["Bin Centers"] >= 25]
#     chi_exp75 = chi_exp.loc[chi_exp["Bin Centers"] <= 75]
#     chi_expend = Exp.loc[Exp["Bin Centers"] >= 75]
#
#     print("Normalization from {} to {}".format(t_min[i],t_max[i]))
#     print("Modified Geant as model")
#     print("Redcued Chi Square from 25 to end: ",ReducedChiSquare(chi_exp["Values"].reset_index(drop=True), chi_frame["Values"].reset_index(drop=True)))
#     print("Redcued Chi Square from 25 to 75: ",ReducedChiSquare(chi_exp75["Values"].reset_index(drop=True), chi_frame75["Values"].reset_index(drop=True)))
#     print("Redcued Chi Square from 75 to end: ",ReducedChiSquare(chi_expend["Values"].reset_index(drop=True), chi_frameend["Values"].reset_index(drop=True)))
#     fig,ax = plt.subplots()
#     ax.step(Exp["Bin Centers"],Exp.Values, color="C0", label="PET Data")
#     ax.step(G4["Bin Centers"],G4.Values, color="seagreen", label="Geant4", linestyle="dashed")
#     ax.step(Framework["Bin Centers"],Framework.Values, color="k", label="Modified Geant4", linestyle="dotted")
#     ax.fill_between(Exp["Bin Centers"],Exp.Values, step="pre", alpha=0.3,color='C0')
#     ax.fill_between(G4["Bin Centers"],G4.Values, step="pre", alpha=0.3,color='seagreen')
#     ax.fill_between(Framework["Bin Centers"],Framework.Values, step="pre", alpha=0.3,color='lightcoral')
#     ax.set_ylabel("PET Events [$s^{-1}$]")
#     ax.set_xlabel("Time [s]")
#     ax.legend()
#     ax.set_xlim(-50,175)
#     ax.set_xticks(ticks=[25*i - 50 for i in range(10)])
#     # plt.yticks(ticks=[5*i for i in range(11)])
#     ax.set_ylim(0,1500)
#     ax_twin = ax.twinx()
#     ax_twin.set_ylim(0,1500)
#     ax_twin.set_ylabel("Normalized Activity [Bq]")
#     plt.subplots_adjust(left=0.11, right=0.89, top=0.98,bottom=0.13)
#     plt.show()
#
#     # ExpandedChiSquare(chi_exp["Values"].reset_index(drop=True), chi_g4["Values"].reset_index(drop=True))
#     # ExpandedChiSquare(chi_exp75["Values"].reset_index(drop=True), chi_g475["Values"].reset_index(drop=True))

#PES population vs depth
dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/Andrey_Plots/75.8MeVGraphs/DosePES/"
Dose = pd.read_csv(dir+"Dose.csv", sep=',', header=None, low_memory=False, names=["Bin Centers","Values"])
Dose["Values"] = 1000*Dose["Values"] / (Dose["Values"].sum() / 0.6)

#due to Normalization introduced in previous data set (40th bin of dose and 26th bin of PES set to 1), dose needs to be readjusted to fit in this dataset
dir = "/home/kilroy101/path/Geant4/Expirimental_Data/FLASH/Andrey_Plots/75.8MeVGraphs/PES-spatial/"
C10 = pd.read_csv(dir+"C10Values.csv", sep=',', header=None, low_memory=False,names=["Bin Centers" ,"Values"])
C10["Values"] = C10["Values"] * 10**6

C11 = pd.read_csv(dir+"C11Values.csv", sep=',', header=None, low_memory=False,names=["Bin Centers" ,"Values"])
C11["Values"] = C11["Values"] * 10**6

N13 = pd.read_csv(dir+"N13Values.csv", sep=',', header=None, low_memory=False,names=["Bin Centers" ,"Values"])
N13["Values"] = N13["Values"] * 10**6

O15 = pd.read_csv(dir+"O15Values.csv", sep=',', header=None, low_memory=False,names=["Bin Centers" ,"Values"])
O15["Values"] = O15["Values"] * 10**6

Tot = pd.read_csv(dir+"SumValues.csv", sep=',', header=None, low_memory=False,names=["Bin Centers" ,"Values"])
Tot["Values"] = Tot["Values"] * 10**6

fig, ax = plt.subplots()

ax_twin = ax.twinx()
ax.step(Tot["Bin Centers"],Tot.Values, linestyle=PES_dict.get("Sum")[1], color=PES_dict.get("Sum")[0], label="Sum", alpha=0.3)
ax.step(C11["Bin Centers"],C11.Values, linestyle=PES_dict.get("C11")[1], color=PES_dict.get("C11")[0], label="C11", alpha=0.3)
ax.step(O15["Bin Centers"],O15.Values, linestyle=PES_dict.get("O15")[1], color=PES_dict.get("O15")[0], label="O15", alpha=0.3)
ax.step(C10["Bin Centers"],C10.Values, linestyle=PES_dict.get("C10")[1], color=PES_dict.get("C10")[0], label="C10", alpha=0.3)
ax.step(N13["Bin Centers"],N13.Values, linestyle=PES_dict.get("N13")[1], color=PES_dict.get("N13")[0], label="N13", alpha=0.3)

ax.fill_between(Tot["Bin Centers"],Tot.Values, step="pre",color="r",alpha=0.3)
ax.fill_between(C11["Bin Centers"],C11.Values, step="pre",linestyle="dotted", color="navy", alpha=0.3)
ax.fill_between(O15["Bin Centers"],O15.Values, step="pre",linestyle="dotted", color="m", alpha=0.3)
ax.fill_between(C10["Bin Centers"],C10.Values, step="pre",linestyle="--", color="b", alpha=0.4)
ax.fill_between(N13["Bin Centers"],N13.Values, step="pre",linestyle="dotted", color="g",  alpha=0.3)

ax_twin.step(Dose["Bin Centers"], Dose.Values, color='dimgray', label="Dose",alpha=0.5)
ax_twin.fill_between(Dose["Bin Centers"], Dose.Values, step="pre",color='dimgray', alpha=0.5)



ax.set_ylabel("PES Yield/mm [per $10^6$ protons]")
ax.set_xlabel("Depth [mm]")
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
ax_twin.legend(lines + lines2, labels + labels2, loc=0)
ax.set_xlim(0,60)
ax.set_ylim(0,500)
ax_twin.set_ylabel("Dose [$\mu$Gy per $10^6$ protons]")
ax_twin.set_ylim(0,50)
plt.subplots_adjust(left=0.15, right=0.87, top=0.98,bottom=0.13)
plt.show()
