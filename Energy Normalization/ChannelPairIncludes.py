import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from scipy import stats

# Helps old people see the code
size = "x-large"

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

def generateTextX(num,max, subplots=False):
    if 1000 > max > 80 and subplots==False:
        x = [max - 0.08*(i+1)*max for i in range(0,num+1)]
    elif max <= 80 and subplots == False:
        x = [max - 0.15*(i+1)*max for i in range(0,num+1)]
    elif max > 1000 and subplots==False:
        x = [max - 0.0007*(i+1)*max for i in range(0, num+1)]
    elif subplots == True:
        x = [0.75*max - 0.08*(i+1)*max for i in range(0,num+1)]
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

def Hist1D(data, units="", fit =[]):

    if fit == []:
        fit = [np.linspace(data.min(),data.max(),100), 3*[0], 100]
    fig = plt.figure()
    #ax1 = fig.add_subplot(111,xlim=(data.min(),data.max()), ylim=(0,1.2*fit[2]))
    ax1 = fig.add_subplot(111,xlim=(data.min(),data.max()), ylim=(0,25))
    ax1.set_xlabel(data.name+" in "+units,fontsize=size)
    ax1.set_ylabel("Counts", fontsize = size)
    ax1.set_title(data.name,fontsize=size)
    ax1.tick_params(axis='both', which='major', labelsize=13)
    y,x,_ = ax1.hist(data, bins=fit[0],color='b')
    xmax = max(x)
    ymax = max(y)
    #T_x = 0.2*data.min() + 0.75*data.max()
    T_x = generateTextX(4,max(x))
    #T_y = generateTextY(9,1.2*fit[2])
    T_y = generateTextY(8,max(y))
    x_p = np.linspace(data.min(),data.max(),400)
    ax1.text(T_x[1],T_y[1],"Data Stats:")
    ax1.text(T_x[1],T_y[2],"mean: "+format(data.mean(),"7.2f"))
    ax1.text(T_x[1],T_y[3],"std: "+format(data.std(),"7.4f"))
    ax1.text(T_x[1],T_y[4],"count: "+str(data.count()))
    ax1.plot(x_p,gauss(x_p,*fit[1]), color='k')
    ax1.set_ylim(0,1.1*ymax)
    ax1.set_xlim(0,1.1*xmax)
    # ax1.text(T_x,T_y[4],"Fit Paramters:")
    # ax1.text(T_x,T_y[5],"mu: "+format(fit[1][1],"7.4f"))
    # ax1.text(T_x,T_y[6],"sigma: "+format(fit[1][2],"7.4f"))
    # ax1.text(T_x,T_y[7],"A: "+format(fit[1][0],"5.2f"))
    # ax1.text(T_x,T_y[8],"FWHM: "+format(2.3548* fit[1][2],"7.4f"))
    # ax1.text(T_x-0.002,T_y[3]-1.75*(T_y[0] - T_y[1]),data.name+" Mean: "+format(fit[1][1],"6.3f")+" "+units, fontsize='x-large')

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