import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from scipy import stats

def generateTextY(num,max,subplots=False):
    if max > 60 and subplots==False:
        y = [max - 0.04*(i+1)*max for i in range(0,num+1)]
    elif max <= 60 and subplots == False:
        y = [max - 0.06*(i+1)*max for i in range(0,num+1)]
    elif subplots == True:
        y = [0.96*max - 0.08*(i+1)*max for i in range(0,num+1)]
    return y

def gauss(x,A,mu,sigma):
    y = A*np.exp(-(x - mu)**2 / (2 * sigma**2))
    return y

def getChipData(df,ChipID):
    tempL  = pd.DataFrame(columns=["ChipID"])
    tempR  = pd.DataFrame(columns=["ChipID"])
    tempL["ChipID"] = df["ChannelIDL"] // 64
    tempR["ChipID"] = df["ChannelIDR"] // 64
    ChipL = tempL.ChipID.unique()
    ChipR = tempR.ChipID.unique()
    if ChipID in ChipL:
        column = "ChannelIDL"
    else:
        column = "ChannelIDR"

    PCB_data = df[df[column] >= 64*ChipID]
    PCB_data = PCB_data[PCB_data[column] < 64*(ChipID + 2)]
    PCB_data["Geo"+column] = PCB_data["Geo"+column] % 1000
    return PCB_data, column

def gaussFit(data,num_bins = 50):
    low = data.min()
    high = data.max() #high = low * -1
    bins = np.linspace(low,high,num_bins)
    x = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),num_bins-1)
    y1 = data.value_counts(bins=bins,sort=False).tolist()
    x1_p, x1_co = curve_fit(gauss,x,y1, p0=[75,data.mean(),data.std()])
    return [bins,x1_p, max(y1)]



def ArrayHist(Top, Bottom,ChipID,units="",fitT = [], fitB =[]):
    if fitT == [] and fitB == []:
        fitT = [np.linspace(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max()),100),5*[0],100]
        fitB = [np.linspace(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max()),100),5*[0],100]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(211,xlim=(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max())), ylim=(0,1.2*max(fitT[2], fitB[2])),ylabel="Counts",title=Top.name+": Chip ID {}".format(ChipID))
    ax2 = fig1.add_subplot(212,xlim=(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max())), ylim=(0,1.2*max(fitT[2], fitB[2])),xlabel=Top.name+" in "+units,ylabel="Counts")

    x_p = np.linspace(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max()),400)
    T_x = 0.95*min(Top.min(), Bottom.min()) + 0.05*max(Top.max(),Bottom.max())
    T_y = generateTextY(9,1.2*max(fitT[2], fitB[2]))

    ax1.hist(Top, bins=fitT[0],color='b')
    ax1.text(T_x,T_y[0],"Data Stats:")
    ax1.text(T_x,T_y[1],"mean: "+format(Top.mean(),"7.4f"))
    ax1.text(T_x,T_y[2],"std: "+format(Top.std(),"7.4f"))
    ax1.text(T_x,T_y[3],"count: "+str(Top.count()))
    ax1.plot(x_p,gauss(x_p,*fitT[1]), color='k')
    ax1.text(T_x,T_y[4],"Fit Paramters:")
    ax1.text(T_x,T_y[5],"mu: "+format(fitT[1][1],"7.4f"))
    ax1.text(T_x,T_y[6],"sigma: "+format(fitT[1][2],"7.4f"))
    ax1.text(T_x,T_y[7],"A: "+format(fitT[1][0],"5.2f"))
    ax1.text(T_x,T_y[8],"FWHM: "+format(2.3548* fitT[1][2],"7.4f"))
    ax1.text(T_x-0.002,T_y[8]-2*(T_y[0] - T_y[1]),Top.name+" Mean: "+format(fitT[1][1],"7.3f"), fontsize='x-large')

    ax2.hist(Bottom, bins=fitB[0],color='b')
    ax2.text(T_x,T_y[0],"Data Stats:")
    ax2.text(T_x,T_y[1],"mean: "+format(Bottom.mean(),"7.4f"))
    ax2.text(T_x,T_y[2],"std: "+format(Bottom.std(),"7.4f"))
    ax2.text(T_x,T_y[3],"count: "+str(Bottom.count()))
    ax2.plot(x_p,gauss(x_p,*fitB[1]), color='k')
    ax2.text(T_x,T_y[4],"Fit Paramters:")
    ax2.text(T_x,T_y[5],"mu: "+format(fitB[1][1],"7.4f"))
    ax2.text(T_x,T_y[6],"sigma: "+format(fitB[1][2],"7.4f"))
    ax2.text(T_x,T_y[7],"A: "+format(fitB[1][0],"5.2f"))
    ax2.text(T_x,T_y[8],"FWHM: "+format(2.3548* fitB[1][2],"7.4f"))
    ax2.text(T_x-0.002,T_y[8]-2*(T_y[0] - T_y[1]),Bottom.name+" Mean: "+format(fitB[1][1],"7.3f"), fontsize='x-large')

def Hist1D(data, units="", fit =[]):
    if fit == []:
        fit = [np.linspace(data.min(),data.max(),100), 5*[0], 100]

    #T_x = 0.97*data.min() + 0.03*data.max()
    T_x = 0.142
    T_y = generateTextY(9,1.2*fit[2])
    x_p = np.linspace(data.min(),data.max(),400)

    fig = plt.figure()
    ax1 = fig.add_subplot(111,xlim=(data.min(),data.max()), ylim=(0,1.2*fit[2]))
    ax1.set_xlabel(data.name+" in "+units,fontsize=15)
    ax1.set_ylabel("Counts", fontsize = 15)
    ax1.set_title(data.name,fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.hist(data, bins=fit[0],color='b')
    ax1.text(T_x,T_y[0],"Data Stats:")
    ax1.text(T_x,T_y[1],"mean: "+format(data.mean(),"7.4f"))
    ax1.text(T_x,T_y[2],"std: "+format(data.std(),"7.4f"))
    ax1.text(T_x,T_y[3],"count: "+str(data.count()))
    ax1.plot(x_p,gauss(x_p,*fit[1]), color='k')
    # ax1.text(T_x,T_y[4],"Fit Paramters:")
    # ax1.text(T_x,T_y[5],"mu: "+format(fit[1][1],"7.4f"))
    # ax1.text(T_x,T_y[6],"sigma: "+format(fit[1][2],"7.4f"))
    # ax1.text(T_x,T_y[7],"A: "+format(fit[1][0],"5.2f"))
    # ax1.text(T_x,T_y[8],"FWHM: "+format(2.3548* fit[1][2],"7.4f"))
    ax1.text(T_x-0.002,T_y[3]-1.75*(T_y[0] - T_y[1]),data.name+" Mean: "+format(fit[1][1],"6.3f")+" "+units, fontsize='x-large')




def ShowBestCTRArrayHist(df, ChipID):
    PCB_data, column = getChipData(df,ChipID)
    # Best CTR by pixel
    CTRL = PCB_data.groupby("Geo"+column).CTR.min().reset_index()
    CTRL.columns = [column, "CTR"]
    CTRL = CTRL.sort_values(by=column)

    CTRLT = CTRL[CTRL[column] < 64].reset_index(drop=True)
    CTRLB = CTRL[CTRL[column] >= 64].reset_index(drop=True)
    ArrayHist(CTRLT.CTR, CTRLB.CTR, ChipID = ChipID,units="ns",fitT = gaussFit(CTRLT.CTR,num_bins= 15), fitB = gaussFit(CTRLB.CTR,num_bins= 15))

    plt.show()
    return 0

def ShowCTRArrayHist(df, ChipID):
    PCB_data, column = getChipData(df,ChipID)
    # Best CTR by pixel

    CTRLT = PCB_data[PCB_data[column] < 64].reset_index(drop=True)
    CTRLB = PCB_data[PCB_data[column] >= 64].reset_index(drop=True)
    ArrayHist(CTRLT.CTR, CTRLB.CTR, ChipID = ChipID,units="ns",fitT = gaussFit(CTRLT.CTR,num_bins= 60), fitB = gaussFit(CTRLB.CTR,num_bins= 60))
    plt.show()
    return 0

def ShowEResArray(df, ChipID):
    PCB_data, column = getChipData(df,ChipID)
    # ERes by array
    EResT = PCB_data[PCB_data["Geo"+column] < 64].reset_index(drop=True)
    EResT = EResT["E_Res"+column[-1]]
    EResB = PCB_data[PCB_data["Geo"+column] >= 64].reset_index(drop=True)
    EResB = EResB["E_Res"+column[-1]]
    ArrayHist(EResT,EResB,ChipID = ChipID, units="%", fitT = gaussFit(EResT), fitB=gaussFit(EResB))


    # #adjusting the lower subplots so they do not get in the way of the upper plots
    # l, b, w, h = ax2.get_position().bounds
    # ax2.set_position([l,b-0.03,w,h-0.03])
    plt.show()

    return 0

def ShowTotalERes(data):
    # ERes Graph
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111,xlim=(2,10), ylim=(0,3400))
    ax1.set_xlabel("Energy Resolution",fontsize=15)
    ax1.set_ylabel("Counts", fontsize = 15)
    ax1.set_title("Energy Resolution of Active Channel Pairs",fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=13)

    jazz = data.E_ResR
    cabbage = data.E_ResL
    jazz = pd.concat([jazz,cabbage])
    low = 0
    high = low +15 #high = low * -1
    num_bins = 150
    bins = np.linspace(low,high,num_bins)
    x = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),num_bins-1)
    y1 = jazz.value_counts(bins=bins,sort=False).tolist()
    a, loc, scale = stats.skewnorm.fit(jazz)
    x_p = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),400)
    y_fit = 5900*stats.skewnorm.pdf(x_p,a,loc,scale)
    mean, var, skew, kurt = stats.skewnorm.stats(a,loc,scale, moments='mvsk')


    T_x = 7.5 #0.05
    T_y = generateTextY(9,3400)

    ax1.hist(jazz, bins=bins,color='b')
    ax1.plot(x_p,y_fit,color='k')
    ax1.text(T_x,T_y[0],"Data Stats:")
    ax1.text(T_x,T_y[1],"mean: "+format(jazz.mean(),"7.4f"))
    ax1.text(T_x,T_y[2],"std: "+format(jazz.std(),"7.4f"))
    ax1.text(T_x,T_y[3],"count: "+str(len(jazz)))
    # ax1.text(T_x,T_y[4],"Fit Paramters:")
    # ax1.text(T_x,T_y[5],"mean: "+format(stats.skewnorm.mean(a,loc,scale),"5.2f"))
    # ax1.text(T_x,T_y[6],"std: "+format(stats.skewnorm.std(a,loc,scale),"5.2f"))
    # ax1.text(T_x,T_y[7],"skewness: "+format(skew,"5.2f"))
    ax1.text(T_x-0.002,T_y[3]-1.75*(T_y[0] - T_y[1]),"ERes Mean: "+format(stats.skewnorm.mean(a,loc,scale),"7.3f"), fontsize='x-large')
    ax1.set_xlim(4.5,10)
    plt.show()
