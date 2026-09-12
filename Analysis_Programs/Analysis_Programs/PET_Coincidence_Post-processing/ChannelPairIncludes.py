import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from scipy import stats
from iminuit import cost, Minuit
import iminuit as im

# Helps old people see the code
size = "x-large"
plt.style.use('ChannelPairStyles.mplstyle')

def generateTextY(num,max, subplots=False):
    if 100 > max > 40 and subplots==False:
        y = [max - 0.045*(i+1)*max for i in range(0,num+1)]
    elif 1000 > max > 100 and subplots==False:
        y = [max - 0.03*(i+1)*max for i in range(0,num+1)]
    elif max > 1000 and subplots==False:
         y = [max - 0.045*(i+1)*max for i in range(0,num+1)]
    elif 20 <= max <= 40 and subplots == False:
        y = [max - 0.08*(i+1)*max for i in range(0,num+1)]
    elif max < 20 and subplots ==False:
        y = [max -0.03*(i+1)*max for i in range(0,num+1)]
    elif subplots == True:
        y = [0.96*max - 0.15*(i+1)*max for i in range(0,num+1)]
    return y

def generateTextX(num,min, subplots=False):
    if 1000 > min > 80 and subplots==False:
        x = [min + 0.08*(i+1)*min for i in range(0,num+1)]
    elif min <= 80 and subplots == False:
        x = [min + 0.15*(i+1)*min for i in range(0,num+1)]
    elif min > 1000 and subplots==False:
        x = [min + 0.0007*(i+1)*min for i in range(0, num+1)]
    elif subplots == True:
        x = [0.75*min + 0.08*(i+1)*min for i in range(0,num+1)]
    return x

def gauss(x,A,mu,sigma):
    y = A*np.exp(-(x - mu)**2 / (2 * sigma**2))
    return y

# def getChipData(df,ChipID):
#     tempL  = pd.DataFrame(columns=["ChipID"])
#     tempR  = pd.DataFrame(columns=["ChipID"])
#     tempL["ChipID"] = df["ChannelIDL"] // 64
#     tempR["ChipID"] = df["ChannelIDR"] // 64
#     ChipL = tempL.ChipID.unique()
#     ChipR = tempR.ChipID.unique()
#     if ChipID in ChipL:
#         column = "ChannelIDL"
#     else:
#         column = "ChannelIDR"

#     PCB_data = df[df[column] >= 64*ChipID]
#     PCB_data = PCB_data[PCB_data[column] < 64*(ChipID + 2)]
#     PCB_data["Geo"+column] = PCB_data["Geo"+column] % 1000
#     return PCB_data, column

def getChipData(df,SlaveID,ChipID):
    GeoChannelIDmin = 10**4 * SlaveID + 100*ChipID
    GeoChannelIDmax = GeoChannelIDmin + 63
    if SlaveID > 2:
        column = "ChannelIDL"
    else:
        column = "ChannelIDR"
    PCB_data = df[df["Geo"+column] >=GeoChannelIDmin]
    PCB_data = PCB_data[PCB_data["Geo"+column] <=GeoChannelIDmax]
    return [PCB_data, column]

def gaussFit(data,num_bins = 50):
    fit_data = data.copy()

    #This will force a second fit on the data cutting data that is not near the initial fit
    low = fit_data.min()
    high = fit_data.max() #high = low * -1
    #hard cut for the sake of getting the tail from the data removed
    fit_data = fit_data.loc[fit_data <=275]
    for i in range(4):
        bins = np.linspace(low,high,num_bins)
        x = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),num_bins-1)
        y1 = fit_data.value_counts(bins=bins,sort=False).tolist()
        c = cost.LeastSquares(x, y1,np.sqrt(y1), gauss)
        fitter = Minuit(c,A=max(y1), mu=fit_data.mean(), sigma=fit_data.std())
        fitter.limits = [(0,None), (fit_data.mean()-2*fit_data.std(), fit_data.mean()+2*fit_data.std()), (0, 3*fit_data.std())]
        fitter.migrad()
        fitter.hesse()
        params = np.array(fitter.values)
        #making a 3 sigma cut on the initial fit
        fit_data = fit_data.loc[fit_data >= (params[1] - 2.5*params[2])]
        fit_data = fit_data.loc[fit_data <= (params[1] + 2.5*params[2])]
        if i == 0:
            inital_A = fitter.values["A"]

    return [bins,[inital_A, params[1], params[2]], max(y1)]

def SkewedGaussFit(data,num_bins=50):
    low = data.min()
    high = data.max()
    bins = np.linspace(low,high,num_bins)
    x = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),num_bins-1)
    y1 = data.value_counts(bins=bins,sort=False).tolist()
    x_f = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),300)
    a, loc, scale = stats.skewnorm.fit(data)
    print("peak location: ", x[np.argmax(y1)])
    y_fit = np.trapz(y1, x)*stats.skewnorm.pdf(x_f, a, loc, scale)
    return [x_f, y_fit]



def ArrayHist(Top, Bottom,ChipID,units="",fitT = [], fitB =[]):
    if fitT == [] and fitB == []:
        fitT = [np.linspace(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max()),100),5*[0],100]
        fitB = [np.linspace(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max()),100),5*[0],100]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(211,xlim=(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max())), ylim=(0,1.2*max(fitT[2], fitB[2])),ylabel="Counts",title=Top.name+": Chip ID {}".format(ChipID))
    ax2 = fig1.add_subplot(212,xlim=(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max())), ylim=(0,1.2*max(fitT[2], fitB[2])),xlabel=Top.name+" in "+units,ylabel="Counts")

    x_p = np.linspace(min(Top.min(), Bottom.min()),max(Top.max(),Bottom.max()),400)
    T_x = 3*min(Top.min(), Bottom.min()) + 0.05*max(Top.max(),Bottom.max())
    #T_x = generateTextX(9,1.2*max(fitT[2], fitB[2]))
    T_y = generateTextY(9,1.2*max(fitT[2], fitB[2]))

    ax1.hist(Top, bins=fitT[0],color='b')
    ax1.text(T_x,T_y[0],"Data Stats:")
    ax1.text(T_x,T_y[2],"mean: "+format(Top.mean(),"7.4f"))
    ax1.text(T_x,T_y[4],"std: "+format(Top.std(),"7.4f"))
    ax1.text(T_x,T_y[6],"count: "+str(Top.count()))
    ax1.plot(x_p,gauss(x_p,*fitT[1]), color='k')
    ax1.text(T_x,T_y[8],"Fit Paramters:")
    ax1.text(T_x,T_y[10],"mu: "+format(fitT[1][1],"7.4f"))
    ax1.text(T_x,T_y[12],"sigma: "+format(fitT[1][2],"7.4f"))
    ax1.text(T_x,T_y[14],"A: "+format(fitT[1][0],"5.2f"))
    ax1.text(T_x,T_y[16],"FWHM: "+format(2.3548* fitT[1][2],"7.4f"))
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

def ShowPPCount(data, bins=()):
    fig3=plt.figure()
    ax3 = fig3.add_subplot(111,xlim=(data.min(),data.max()), ylim=(0, 100000))
    ax3.set_xlabel("Number of Coincident Events", fontsize=size)
    ax3.set_ylabel("Counts", fontsize=size)
    #bins=np.linespace(73, 140, 23)
    T_x = generateTextX(6,150)
    T_y = generateTextY(14,400)
    ax3.text(T_x[1],T_y[2],"mean: "+format(data.mean(),"7.4f"))
    ax3.text(T_x[1],T_y[3],"std: "+format(data.std(),"7.4f"))
    ax3.text(T_x[1],T_y[4],"count: "+str(data.count()))

def Hist1D(data, units="", fit = [], num_bins = 0):
    fig = plt.figure()

    #ax1 = fig.add_subplot(111,xlim=(data.min(),data.max()), ylim=(0,1.2*fit[2]))
    ax1 = fig.add_subplot(111,xlim=(data.min(),data.max()), ylim=(0,25))
    ax1.set_xlabel(data.name+" in "+units,fontsize=size)
    ax1.set_ylabel("Counts", fontsize = size)
    ax1.set_title(data.name,fontsize=size)

    # ax1.tick_params(axis='both', which='major', labelsize=13)

    if fit == []:
        numbins = num_bins
        fit = [np.linspace(data.min(),data.max(),100), 3*[0], 100]
        y,x,_ = ax1.hist(data, bins=numbins,color='b')
    else:
        x_p = np.linspace(data.min(),data.max(),400)
        y,x,_ = ax1.hist(data, bins=fit[0],color='b')
        xmax = max(x)
        ymax = max(y)
        #T_x = 0.2*data.min() + 0.75*data.max()
        T_x = [0,175, 300]
        #T_y = generateTextY(9,1.2*fit[2])
        T_y = generateTextY(10,max(y))
        ax1.plot(x_p,gauss(x_p,*fit[1]), color='k')
        # ax1.text(T_x[1],T_y[5],"Fit Paramters:")
        # ax1.text(T_x[1],T_y[6],"mu: "+format(fit[1][1],"7.4f"))
        # ax1.text(T_x[1],T_y[7],"sigma: "+format(fit[1][2],"7.4f"))
        # ax1.text(T_x[1],T_y[8],"A: "+format(fit[1][0],"5.2f"))
        # ax1.text(T_x[1],T_y[9],"FWHM: "+format(2.3548* fit[1][2],"7.4f"))
        ax1.text(T_x[2]-0.002,T_y[1],data.name+" Mean: "+format(fit[1][1],"5.1f")+" "+units, fontsize='x-large')


    xmax = max(x)
    ymax = max(y)
    #T_x = 0.2*data.min() + 0.75*data.max()
    T_x = [0,175, 300]
    #T_y = generateTextY(9,1.2*fit[2])
    T_y = generateTextY(10,max(y))
    ax1.text(T_x[1],T_y[1],"Data Stats:")
    ax1.text(T_x[1],T_y[2],"mean: "+format(data.mean(),"7.2f"))
    ax1.text(T_x[1],T_y[3],"std: "+format(data.std(),"7.4f"))
    ax1.text(T_x[1],T_y[4],"count: "+str(data.count()))
    ax1.set_ylim(0,1.1*ymax)
    ax1.set_xlim(0,1.1*xmax)


def ShowBestCTRArrayHist(df, SlaveID, ChipID):
    PCB_data, column = getChipData(df,SlaveID, ChipID)
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
    EResT = EResT["ERres"+column[-1]]
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
    #ax1 = fig1.add_subplot(111,xlim=(2,10), ylim=(0,3400))
    ax1 = fig1.add_subplot(111,xlim=(2,10), ylim=(0,20))
    ax1.set_xlabel("Energy Resolution",fontsize=15)
    ax1.set_ylabel("Counts", fontsize = 15)
    ax1.set_title("Energy Resolution of Active Channel Pairs",fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=13)

    jazz = data.ERres
    cabbage = data.ELres
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
    ax1.text(T_x,T_y[1],"mean: "+format(jazz.mean(),"4.1f"))
    ax1.text(T_x,T_y[2],"std: "+format(jazz.std(),"4.1f"))
    ax1.text(T_x,T_y[3],"count: "+str(len(jazz)))
    # ax1.text(T_x,T_y[4],"Fit Paramters:")
    # ax1.text(T_x,T_y[5],"mean: "+format(stats.skewnorm.mean(a,loc,scale),"5.2f"))
    # ax1.text(T_x,T_y[6],"std: "+format(stats.skewnorm.std(a,loc,scale),"5.2f"))
    # ax1.text(T_x,T_y[7],"skewness: "+format(skew,"5.2f"))
    ax1.text(T_x-0.002,T_y[3]-1.75*(T_y[0] - T_y[1]),"ERes Mean: "+format(stats.skewnorm.mean(a,loc,scale),"7.3f"), fontsize='x-large')
    ax1.set_xlim(4.5,10)
    plt.show()

def Scatter(df,xName,yName,xUnit = '',yUnit = '',savefig = False):
        data = df[0]
        channel_name = df[1]
        # arrays are better than series...
        x = data[xName].to_numpy()
        y = data[yName].to_numpy()
        # generate plots
        fig,ax = plt.subplots()
        ax.scatter(x,y)
        ymin,ymax = plt.ylim()
        xmin,xmax = plt.xlim()
        y_place = generateTextY(10,ymax)
        x_place = generateTextX(0,xmin)
        # getting staistics and adding them to the figure:
        yMean,xMean = str(np.round(np.mean(y),2)),str(np.round(np.mean(x),2))
        yStdev,xStdev = str(np.round(np.std(y),2)),str(np.round(np.std(x),2))
        counts = str(len(y))
        ax.text(xmin, y_place[0],r'$\bf{Statistics:}$' + '\n')
        ax.text(xmin, y_place[1],'%s Mean: ' % xName + "%s " % xMean + "%s" % xUnit)
        ax.text(xmin, y_place[2],'%s std: ' % xName + "%s " % xStdev + "%s" % xUnit)
        ax.text(xmin, y_place[3],'%s Mean: ' % yName + "%s " % yMean + "%s" % yUnit)
        ax.text(xmin, y_place[4],'%s std: ' % yName + "%s " % yStdev + "%s" % yUnit)
        ax.text(xmin, y_place[5],'%s Counts: ' % counts)
        # setting x and y labels based on if units were given or not
        if len(xUnit) > 0:
            plt.xlabel(xName + " (%s)" % xUnit)
        else:
            plt.xlabel(xName)
        if len(yUnit) > 0:
            plt.ylabel(yName + " (%s)" % yUnit)
        else:
            plt.ylabel(yName)
        # other plot formatting
        plt.title(yName + ' vs ' + xName)
        # # save plot if True ; default is False
        # if savefig == True:
        #     plt.savefig(‘foo.png’)
        #     plt.close()
        return 0


def ChipScatter(ChipData,column,OneToOne = False, weighted = False):
  if weighted == OneToOne and OneToOne == True:
      raise ValueError("Hey idiot you can't do a weighted AND most populated heatmap")

  GeoID = "Geo"+ChipData[1]
  ChipData = ChipData[0]
  #WLOG we can take the minumum ChannelID and get the slave, same for the chip
  SlaveID = ChipData[GeoID].min() // 10000
  ChipID =  ChipData[GeoID].min() // 100 - 100*SlaveID
  #If One_to_One = True then we only plot most populated
  fig = plt.figure()
  if OneToOne == True:
      #Sort Data from highest Count to lowest for each channel and save only most channel data with most counts
      Sorted_Data = ChipData.sort_values(by = "PPCount", ascending = False)
      groups = Sorted_Data.groupby([GeoID])[column].first().reset_index()
      MEAN = np.mean(groups[GeoID] % 64)
      MIN = min(groups[column])
      MAX = max(groups[column])
      hops = max(groups[GeoID] % 64)

      #Plotting

      ax1 = fig.add_subplot(111, ylim=(MIN,MAX),xlim=(MEAN-40,MEAN+40), xlabel="GeoChannelID", ylabel=column, title="Most Populated CP Scatter plot of {}".format(column))
      ax1.scatter(groups[GeoID] % 64, groups[column], marker="_")
      #Stat Box
      T_x = generateTextX(10, hops) #0.05
      T_y = generateTextY(10,MAX)
      ax1.text(T_x[0],T_y[0],"Data Stats:")
      ax1.text(T_x[0],T_y[1],"mean: "+format(groups[column].mean(),"7.1f"))
      ax1.text(T_x[0],T_y[2],"std: "+format(groups[column].std(),"7.2f"))
      ax1.text(T_x[0],T_y[3],"count: "+str(len(groups)))

  elif weighted == True:
      ScatterData = pd.DataFrame(columns=[GeoID, "weighted_avg"])
      data = ChipData
      data["weighted"] = data["PPCount"] * data[column]
      groups = data.groupby(GeoID)
      for group in groups:
          if group[1][GeoID].count() != 0:
              df = group[1].reset_index(drop=True)
              weighted_avg = df["weighted"].sum() / df["PPCount"].sum()
              new_CID = pd.DataFrame({GeoID:[group[0]], "weighted_avg":[weighted_avg]})
              ScatterData = pd.concat([ScatterData, new_CID], ignore_index=True)
      groups = ChipData.sort_values(by = column)
      ScatterDataMod64 = ScatterData[GeoID]%64
      # print (ScatterDataMod64)
      # print (ScatterData[GeoID])
      # print (ScatterData["weighted_avg"])
      #Defining Mean for x-limit and Max and min of ERes for y-limit
      MEAN = np.mean(groups[GeoID] % 64)
      MIN = min(groups[column])
      MAX = max(groups[column])
      hops = max(groups[GeoID] % 64)
      ax1 = fig.add_subplot(111, ylim=(MIN,MAX),xlim=(MEAN-40,MEAN+40), xlabel="GeoChannelID", ylabel=column, title="Weighted Average Scatter plot of {} \n for Slave {} Chip {}".format(column, int(SlaveID), int(ChipID)))
      ax1.scatter(ScatterDataMod64, ScatterData["weighted_avg"], marker="_")
      print (ScatterDataMod64)
      #Stat Box
      T_x = generateTextX(10, hops) #0.05
      T_y = generateTextY(10,MAX)
      ax1.text(T_x[0],T_y[0],"Data Stats:")
      ax1.text(T_x[0],T_y[1],"mean: "+format(groups[column].mean(),"7.1f"))
      ax1.text(T_x[0],T_y[2],"std: "+format(groups[column].std(),"7.2f"))
      ax1.text(T_x[0],T_y[3],"count: "+str(len(ScatterData)))

  else:
      #For all data in Chip
      groups = ChipData.sort_values(by = column)
      #Defining Mean for x-limit and Max and min of ERes for y-limit
      MEAN = np.mean(groups[GeoID])
      MIN = min(groups[column])
      MAX = max(groups[column])
      #Plotting
      ax1 = fig.add_subplot(111, ylim=(MIN,MAX),xlim=(MEAN-40,MEAN+40), xlabel="GeoChannelID", ylabel=column, title="Scatter plot of {}".format(column))#ylim=(4,10)
      ax1.scatter(groups[GeoID] % 64,groups[column], marker="_")
      #Stat Box
      T_x = generateTextX(10, hops) #0.05
      T_y = generateTextY(10,MAX)
      ax1.text(T_x[0],T_y[0],"Data Stats:")
      ax1.text(T_x[0],T_y[1],"mean: "+format(groups[column].mean(),"7.4f"))
      ax1.text(T_x[0],T_y[2],"std: "+format(groups[column].std(),"7.4f"))
      ax1.text(T_x[0],T_y[3],"count: "+str(len(groups)))



#The general heatmap class for plotting a single quantity as a colorbar for a single chipID. The squares on this heatmap will be the GeoChannelIDs
def ChipHeatmap(ChipData,column,weighted=False, OneToOne= False,units = ""):
    if weighted == OneToOne:
        raise ValueError("Hey idiot you can't do a weighted AND most populated heatmap")
    #extracting the two relevant columns: The ChannelID and the variable that will be the color axis
    GeoColumn = "Geo"+ChipData[1]
    Chipdata = ChipData[0][[GeoColumn,str(column), "PPCount"]]
    #WLOG we can take the minumum ChannelID and get the slave, same for the chip
    SlaveID = Chipdata[GeoColumn].min() // 10000
    ChipID =  Chipdata[GeoColumn].min() // 100 - 100*SlaveID
    #The range of all possible GeoChannelIDs in the dataset
    GeoColumnRange = range(int(10**4 * SlaveID + 100*ChipID), int(10**4 * SlaveID + 100*ChipID + 64))
    #If I wanted to plot the weighted average of the column on the heatmap
    HMapdata = pd.DataFrame(columns = [GeoColumn, str(column)])
    for i in GeoColumnRange:
        #Getting the data for a specific ChannelID
        IDSlice = Chipdata[Chipdata[GeoColumn] == i]
        #If the channelID doesn't appear in the list
        if IDSlice[GeoColumn].count() == 0:
            addedID = pd.DataFrame({GeoColumn:[i], str(column): [0]})
            #add it to the dataframe with a zero value
            HMapdata = pd.concat([HMapdata, addedID], ignore_index=True)
        else:
            #calculates the weighted average for all channel pairs sharing a GeochannelID
            if weighted:
                IDSlice["weighted_data"] = IDSlice[column]*IDSlice["PPCount"]
                total_values = IDSlice["weighted_data"].sum()
                total_weights = IDSlice["PPCount"].sum()
                addedID = pd.DataFrame({GeoColumn:[i], str(column): [total_values/total_weights]})
                HMapdata = pd.concat([HMapdata, addedID], ignore_index=True)
            #takes the value of the column from the most populated channel pair per GeochannelID
            elif OneToOne:
                most_pop = IDSlice.sort_values(by="PPCount").reset_index()
                most_pop_val = most_pop[str(column)][0]
                addedID = pd.DataFrame({GeoColumn:[i], str(column):[most_pop_val]})
                pd.concat([HMapdata, addedID], ignore_index=True)
    #getting the GeoChannelIDs in order
    HMapdata = HMapdata.sort_values(by=GeoColumn)
    #numpy goodness to get our data in to the right square format
    ColumnData = np.reshape(np.array(HMapdata[column]), (8,8))
    plt.figure()
    plt.xticks(np.arange(8))
    plt.yticks(np.arange(8))
    for i in range(8):
        for j in range(8):
            #putting the value in text of the column in question
            plt.text(j, 7-i, format(ColumnData[i][j],"3.0f"), ha="center", va="center", color="b", fontsize='x-large')
            #debug to show the index that should convert to GeoChannelID on the plane
            #plt.text(j, 7-i, str(j*(7-i)), ha="center", va="center", color="b", fontsize='x-large')
    cmap = plt.get_cmap('binary')
    cmap = truncate_colormap(cmap, 0.2, 0.8)
    #Showing the heatmap
    #vmin and vmax are asymetric because
    plt.imshow(ColumnData,cmap=cmap,vmin=(HMapdata[str(column)].mean() - 3*HMapdata[str(column)].std()), vmax=(HMapdata[str(column)].mean() + 3*HMapdata[str(column)].std()))
    plt.colorbar()
    plt.title("Heatmap of {} for Slave {}, Chip {}".format(str(column), int(SlaveID), int(ChipID)))


def FullScannerHeatmap(data):
    #Plotting the charge weighted average CTR of teh full scanner by channel. To change the variable of interest, simply replace the following block of code
    #Charge weighted CTR
    data["CWCTR"] = data["CTR"] * data["PPCount"]

    CWCTRL = data.groupby("GeoChannelIDL")["CWCTR"].sum().reset_index()
    CWCTRL.columns = ["GeoChannelIDL", "CWCTRsum"]

    SumL = data.groupby("GeoChannelIDL")["PPCount"].sum().reset_index()
    SumL.columns = ["GeoChannelIDL", "PPsum"]

    dataL = CWCTRL.merge(SumL, on="GeoChannelIDL")
    dataL["CWCTR"] = dataL["CWCTRsum"] / dataL["PPsum"]

    CWCTRR = data.groupby("GeoChannelIDR")["CWCTR"].sum().reset_index()
    CWCTRR.columns = ["GeoChannelIDR", "CWCTRsum"]

    SumR = data.groupby("GeoChannelIDR")["PPCount"].sum().reset_index()
    SumR.columns = ["GeoChannelIDR", "PPsum"]

    dataR = CWCTRR.merge(SumR, on="GeoChannelIDR")
    dataR["CWCTR"] = dataR["CWCTRsum"] / dataR["PPsum"]

    #properly formatting the channelID values so they increment by one every time for the purposes of graphing
    def XCol(ChannelID):
        Port = ChannelID // 10**6
        Slave = (ChannelID - 10**6*Port) // 10**4
        Chip = (ChannelID - Slave*10**4 - 10**6 * Port) // 100
        Channel = ChannelID % 100
        if Port == 1 and Chip < 8:
            XCol = 32*(2 - Slave) + 8*(Chip // 2 ) + (Channel%8)
            return XCol
        elif Port == 0 and Chip < 8:
            XCol = 32*(2 - Slave) + 8*(3 - Chip//2 ) + (Channel%8)
            return XCol
        elif Port == 1 and Chip >= 8:
            XCol = 32*(2 - Slave) + 8*((15-Chip) // 2 ) + (Channel%8)
            return XCol
        else:
            XCol = 32*(2 - Slave) + 8*(Chip // 2 -4) + (Channel%8)
            return XCol

    #getting the row of the channel ID to only plot one row at a time on the histograms
    def ZRow(ChannelID):
        Port = ChannelID // 10**6
        Slave = (ChannelID - 10**6*Port) // 10**4
        Chip = (ChannelID - Slave*10**4 - 10**6 * Port) // 100
        Channel = ChannelID % 100
        if Chip < 8:
            Row = 8*(Chip % 2) + (Channel // 8)
            return Row
        else:
            Row = 16 + 8*((Chip -8) % 2) + (Channel // 8)
            return Row

    #preparing the unique heatmap index
    dataL["Row"] = dataL["GeoChannelIDL"].apply(ZRow)
    dataL["Col"] = dataL["GeoChannelIDL"].apply(XCol)

    dataR["Row"] = dataR["GeoChannelIDR"].apply(ZRow)
    dataR["Col"] = dataR["GeoChannelIDR"].apply(XCol)

    dataL["HMapIndex"] = 96*dataL["Row"] + dataL["Col"]
    dataR["HMapIndex"] = 96*dataR["Row"] + dataR["Col"]
    #filling in dead pixels, 3072 = 32*96
    for i in range(3072):
        if i not in dataL.HMapIndex.unique():
            fillerL = pd.DataFrame({"GeoChannelIDL":[0], "CWCTRsum":[0], "PPsum":[0], "CWCTR": [0], "Row":[0], "Col":[0], "HMapIndex":[i]})
            dataL = pd.concat([dataL, fillerL], ignore_index=True)
        if i not in dataR.HMapIndex.unique():
            fillerR = pd.DataFrame({"GeoChannelIDR":[0], "CWCTRsum":[0], "PPsum":[0], "CWCTR": [0],"Row":[0], "Col":[0], "HMapIndex":[i]})
            dataR = pd.concat([dataR, fillerR], ignore_index=True)

    dataL = dataL.sort_values(by="HMapIndex").reset_index(drop=True)
    dataR = dataR.sort_values(by="HMapIndex").reset_index(drop=True)

    crescent_mapL = np.reshape(np.array(dataL.CWCTR), (32,96))
    crescent_mapR = np.reshape(np.array(dataR.CWCTR), (32,96))
    index_test = np.reshape(np.array(dataL.index), (32,96))

    # for i in range(32):
    #     for j in range(96):
    #         plt.text(j, i, crescent_mapL[i][j], ha="center", va="center", color="k", fontsize='small')

    fig = plt.figure()
    ax = fig.add_subplot(211)
    im = ax.imshow(crescent_mapL,vmin=100, vmax=500)
    # cax = plt.axes([0.85, 0.1, 0.025, 0.8])
    #fig.colorbar(im, label='CTR [ps]', orientation='horizontal') #cax = cax
    ax.set_title("Left Crescent")
    #plt.subplots_adjust(bottom=0.05, right=0.8, top=0.95,left=0.2)


    # fig1 = plt.figure()
    ax1 = fig.add_subplot(212)
    im1 = ax1.imshow(crescent_mapR,vmin=100, vmax=500)
    cax = plt.axes([0.85, 0.1, 0.025, 0.8])
    cbar = fig.colorbar(im1,cax = cax)
    cbar.set_label('CTR [ps]', size=35)
    ax1.set_title("Right Crescent")
    plt.subplots_adjust(bottom=0.05, right=0.8, top=0.95,left=0.2)
    plt.show()

    return 0
