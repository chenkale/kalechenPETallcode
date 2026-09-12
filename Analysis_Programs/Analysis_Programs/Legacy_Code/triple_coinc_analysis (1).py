# Takes in triple coincidence data for double-ended readout setup with reference module.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
from scipy.optimize import curve_fit
from tqdm import tqdm
import statistics as stat

pd.set_option("display.max_columns",100)

class Channel_Pair:
    def __init__(self, ChannelPairData):
        self.df = ChannelPairData
        self.df["ChargeTot"] = self.df["ChargeR"] + self.df["ChargeL"]
        valuesR,binsR,presentersR = plt.hist(self.df.ChargeR, bins = np.linspace(0,29,80), ec='r',label='Right Channel', fill=False)
        valuesL,binsL,presentersL = plt.hist(self.df.ChargeL, bins = np.linspace(0,29,80), ec='b',label='Left Channel', fill=False)
        valuesT,binsT,presentersT = plt.hist(self.df.ChargeTot, bins = np.linspace(0,40,int((40*80)/29)), ec='g',label='Channel Sum', fill=False)
        bin_centers = [(binsR[q] + binsR[q+1])/2 for q in range(0,len(binsR)-1)]
        bin_centersT = [(binsT[q] + binsT[q+1])/2 for q in range(0,len(binsT)-1)]
        self.EHistR_y = valuesR.tolist()
        self.EHistR_x = bin_centers
        self.EHistL_y = valuesL.tolist()
        self.EHistL_x = bin_centers
        self.EHistT_y = valuesT.tolist()
        self.EHistT_x = bin_centersT
        self.PhotopeakR = 0
        self.PhotopeakL = 0
        self.PhotopeakT = 0
        self.EResR = 0
        self.EResL = 0
        self.EResT = 0
        self.PP_Count = 0
        self.EFitParam = 0
        self.CTRFitParam = 0
        self.DOI_mean = 0
        self.DOI_FWHM = 0
        plt.clf()
        plt.close()

    def findPhotopeak(self):
        #find the peak within 3 bins starting from the RHS of the graph
        for i in range(1, len(self.EHistR_y)-4):
            if self.EHistR_y[-1*i] > self.EHistR_y[-1*(i+1)] and self.EHistR_y[-1*i] > self.EHistR_y[-1*(i+2)] and self.EHistR_y[-1*i] > self.EHistR_y[-1*(i+3)] and self.EHistR_y[-1*(i+1)] != 0 and self.EHistR_y[-1*i] >= 30:
                self.PhotopeakR = self.EHistR_y.index(self.EHistR_y[-1*i],-1*i-1)
                break
        for i in range(1, len(self.EHistL_y)-4):
            if self.EHistL_y[-1*i] > self.EHistL_y[-1*(i+1)] and self.EHistL_y[-1*i] > self.EHistL_y[-1*(i+2)] and self.EHistL_y[-1*i] > self.EHistL_y[-1*(i+3)] and self.EHistL_y[-1*(i+1)] != 0 and self.EHistL_y[-1*i] >= 30:
                self.PhotopeakL = self.EHistL_y.index(self.EHistL_y[-1*i],-1*i-1)
                break
        for i in range(1, len(self.EHistT_y)-4):
            if self.EHistT_y[-1*i] > self.EHistT_y[-1*(i+1)] and self.EHistT_y[-1*i] > self.EHistT_y[-1*(i+2)] and self.EHistT_y[-1*i] > self.EHistT_y[-1*(i+3)] and self.EHistT_y[-1*(i+1)] != 0 and self.EHistT_y[-1*i] >= 30:
                self.PhotopeakT = self.EHistT_y.index(self.EHistT_y[-1*i],-1*i-1)
                break

    def getQINTs(self):
        # numpy requires the arguments y before x
        INTR = np.trapz(self.EHistR_y,self.EHistR_x)
        INTL = np.trapz(self.EHistL_y,self.EHistL_x)
        return [INTR, INTL]

    def PPFit(self):
        ExR_p = [0,1,1]
        ExL_p = [0,1,1]
        ExT_p = [0,1,1]

        fitR_x = self.EHistR_x[self.PhotopeakR-7:]
        fitL_x = self.EHistL_x[self.PhotopeakL-7:]
        fitT_x = self.EHistT_x[self.PhotopeakT-7:]
        fitR_y = self.EHistR_y[self.PhotopeakR-7:]
        fitL_y = self.EHistL_y[self.PhotopeakL-7:]
        fitT_y = self.EHistT_y[self.PhotopeakT-7:]


        # plt.scatter(fitT_x, fitT_y)
        # plt.show()

        #fitting said data and recording photpeak energy cuts
        if len(fitR_y) != 0 and len(fitR_x) != 0:
            try:
                ExR_p, ExR_co = curve_fit(gauss, fitR_x, fitR_y, p0=[max(fitR_y), fitR_x[fitR_y.index(max(fitR_y))], 0.5], bounds=[[max(fitR_y)-50, fitR_x[fitR_y.index(max(fitR_y))]-2, 0.001],[max(fitR_y)+50, fitR_x[fitR_y.index(max(fitR_y))]+2, 1]])
            except RuntimeError:
                print("Error - curve_fit failed")

        if len(fitL_x) != 0 and len(fitL_y) != 0:
            try:
                ExL_p, ExL_co = curve_fit(gauss, fitL_x, fitL_y, p0=[max(fitL_y), fitL_x[fitL_y.index(max(fitL_y))], 0.5], bounds=[[max(fitL_y)-50, fitL_x[fitL_y.index(max(fitL_y))]-2, 0.001],[max(fitL_y)+50, fitL_x[fitL_y.index(max(fitL_y))]+2, 1]])
            except RuntimeError:
                print("Error - curve_fit failed")

        if len(fitT_x) != 0 and len(fitT_y) != 0:
            try:
                ExT_p, ExT_co = curve_fit(gauss, fitT_x, fitT_y, p0=[max(fitT_y), fitT_x[fitT_y.index(max(fitT_y))], 0.5], bounds=[[max(fitT_y)-50, fitT_x[fitT_y.index(max(fitT_y))]-2, 0.001],[max(fitT_y)+50, fitT_x[fitT_y.index(max(fitT_y))]+2, 1]])
            except RuntimeError:
                print("Error - curve_fit failed")

        return [ExR_p, ExL_p, ExT_p]

    def SetEnergyParameters(self, run_param):
        self.EFitParam = run_param
        self.EResR = abs(100*(2.3548*abs(self.EFitParam[0][2]))/self.EFitParam[0][1])
        self.EResL = abs(100*(2.3548*abs(self.EFitParam[1][2]))/self.EFitParam[1][1])
        self.EResT = abs(100*(2.3548*abs(self.EFitParam[2][2]))/self.EFitParam[2][1])

    def CutOnEnergy(self,normed=True):
        if normed == True:
            R_cut = self.EFitParam[0][1] - 2.5*self.EFitParam[0][2]
            L_cut = self.EFitParam[0][1] - 2.5*(self.EFitParam[0][1]/self.EFitParam[1][1])*self.EFitParam[1][2]
            T_cut = self.EFitParam[2][1] - 2.5*self.EFitParam[2][2]
        elif normed == False:
            R_cut = self.EFitParam[0][1] - 2.5*self.EFitParam[0][2]
            L_cut = self.EFitParam[1][1] - 2.5*self.EFitParam[1][2]
            T_cut = self.EFitParam[2][1] - 2.5*self.EFitParam[2][2]
        print(R_cut)
        print(L_cut)
        print(self.df)
        #Changing wether the energy is cut on seperate L and R photopeaks or the summed photopeak
        self.df = self.df[self.df["ChargeR"] >= R_cut]
        self.df = self.df[self.df["ChargeL"] >= L_cut].reset_index(drop=True)

        # self.df = self.df[self.df["ChargeTot"] >= T_cut].reset_index(drop=True)
        self.df["Time_diff"] = 10**-3 * (self.df["TimeL"] - self.df["TimeR"])

    def EnergyCut(self):
        if self.EFitParam[0][1] > 5 and self.EFitParam[1][1] > 5: #and self.EFitParam[0][2] < 200 and self.EFitParam[1][2] < 200:
            return True
        else:
            return False

    def PPCut(self):
        if self.df["Time_diff"].count() < 100:
            return False
        else:
            self.PP_Count = self.df["Time_diff"].count()
            return True

    def SetCTRParameters(self, CTR_param):
        self.CTRFitParam = CTR_param

    def CTRFit(self):
        TRes5_UN = self.df["Time_diff"]

        x1_p = [0,1,1]
        low = -2
        high = low * -1
        num_bins = 150
        bins = np.linspace(low,high,num_bins)
        x = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),num_bins-1)
        y1 = TRes5_UN.value_counts(bins=bins,sort=False).tolist()
        peak = y1.index(max(y1))
        bin_centers = [(bins[i] + bins[i+1])/2 for i in range(0,len(bins)-1)]
        fit_x = bin_centers[peak-10:peak+11]
        fit_y = y1[peak-10:peak+11]
        try:
            x1_p, x1_co = curve_fit(gauss,fit_x,fit_y, p0=[max(fit_y),stat.mean(fit_x),0.5*stat.pstdev(fit_x)], bounds=[[max(fit_y)-10,stat.mean(fit_x)-0.1,0],[max(fit_y)+10,stat.mean(fit_x)+0.1,0.2]])
        except RuntimeError:
            print("Error - curve_fit failed")

        return x1_p


    def PlotEnergySpectrum(self, display=False):
        global LChannel
        global RChannel
        # plotting every energy spectrum
        x_f = np.linspace(0,29,300)
        plt.hist(self.df.ChargeR, bins = np.linspace(0,29,80), ec='r',label='Right Channel', fill=False)
        plt.hist(self.df.ChargeL, bins = np.linspace(0,29,80), ec='b',label='Left Channel', fill=False)
        # plt.hist(self.df.ChargeTot, bins = np.linspace(0,40,int((80*40)/29)), ec='g',label='Channel Sum', fill=False)
        plt.plot(x_f, gauss(x_f,*self.EFitParam[0]), color='k')
        # Displays the data shifted but with the inital energy paramters plotted
        plt.plot(x_f, gauss(x_f,*self.EFitParam[1]), color='k')
        # # Displats the data and fit funtion shifted, ONLY USE FOR DOI RUN IN CENTER OF CRYSTAL
        #plt.plot(x_f, gauss(x_f, self.EFitParam[1][0], self.EFitParam[0][1], self.EFitParam[1][2]), color='k')
        # plt.plot(x_f, gauss(x_f, *self.EFitParam[2]), color='k')
        t_y = generateTextY(10,1200)#400
        t_x = 1
        t_xT = 12
        plt.text(t_x,t_y[0], "Right Channel:")
        plt.text(t_x,t_y[1], "mean: "+format(self.EFitParam[0][1], "5.2f"))
        plt.text(t_x,t_y[2], "std: "+format(abs(self.EFitParam[0][2]), "5.2f"))
        plt.text(t_x,t_y[3], "count: "+format(ChannelPair.df.ChargeR.count(), "5"))
        plt.text(t_x,t_y[4], "ERes: "+format(self.EResR, "5.2f")+"%")
        plt.text(t_x,t_y[5], "Left Channel:")
        plt.text(t_x,t_y[6], "mean: "+format(self.EFitParam[1][1], "5.2f"))
        plt.text(t_x,t_y[7], "std: "+format(abs(self.EFitParam[1][2]), "5.2f"))
        plt.text(t_x,t_y[8], "count: "+format(ChannelPair.df.ChargeL.count(), "5"))
        plt.text(t_x,t_y[9], "ERes: "+format(self.EResL, "5.2f")+"%")
        plt.text(t_xT,t_y[0], "Channel Sum:")
        plt.text(t_xT,t_y[1], "mean: "+format(self.EFitParam[2][1], "5.2f"))
        plt.text(t_xT,t_y[2], "std: "+format(abs(self.EFitParam[2][2]), "5.2f"))
        plt.text(t_xT,t_y[3], "count: "+format(ChannelPair.df.ChargeTot.count(), "5"))
        plt.text(t_xT,t_y[4], "ERes: "+format(self.EResT, "5.2f")+"%")
        plt.ylim(0,1200)
        plt.title("Energy Spectrum, Pair {}, {}, {} mm DOI".format(LChannel, RChannel, Location))
        plt.legend()
        plt.savefig("Energy Spectrum, Pair {}, {}, {} mm DOI".format(LChannel,RChannel, Location))
        if display == True:
            plt.show()
            plt.clf()
            plt.close()
        else:
            plt.clf()
            plt.close()

    def PlotCTR(self, display=False):
        global LChannel
        global RChannel
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111,xlim=(-2,2), ylim=(0,20),xlabel="Time Difference in ns",ylabel="Counts",title="CTR: Channel Pair {}, {}".format(LChannel,RChannel))
        low = -2
        high = low * -1
        num_bins = 156
        bins = np.linspace(low,high,num_bins)
        x_p = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),400)
        T_x = -0.9#-1.9
        T_y = generateTextY(9,100)#100
        TRes5_UN = self.df["Time_diff"]
        ax1.hist(TRes5_UN, bins=bins,color='b')
        ax1.plot(x_p,gauss(x_p,*self.CTRFitParam), color='k')
        ax1.text(T_x,T_y[0],"Data Stats:")
        ax1.text(T_x,T_y[1],"mean: "+format(TRes5_UN.mean(),"5.2f"))
        ax1.text(T_x,T_y[2],"std: "+format(TRes5_UN.std(),"7.4f"))
        ax1.text(T_x,T_y[3],"count: "+str(len(TRes5_UN)))
        ax1.text(T_x,T_y[4],"Fit Paramters:")
        ax1.text(T_x,T_y[5],"mu: "+format(self.CTRFitParam[1],"5.2f"))
        ax1.text(T_x,T_y[6],"sigma: "+format(self.CTRFitParam[2],"7.4f"))
        ax1.text(T_x,T_y[7],"A: "+format(self.CTRFitParam[0],"5.2f"))
        ax1.text(T_x,T_y[8],"FWHM: "+format(2.3548* self.CTRFitParam[2],"7.4f"))
        # plt.savefig("CTR, Pair {}, {}".format(LChannel, RChannel))

        if int(LChannel) == 64 and int(RChannel) == 71:
            Hist = pd.DataFrame(columns = ['TimeDiff'])
            Hist["TimeDiff"] = self.df.Time_diff
            # Hist.to_csv(dir+"TimeDiff{}.csv".format(Location))

        if display == True:
            plt.show()
            plt.clf()
            plt.close()
        else:
            plt.clf()
            plt.close()

    def PlotCountsDiff(self, display=False):
        global LChannel
        global RChannel
        global Location
        #Define our Normalized Count diffferences
        self.df["NormCountDiff"] = (self.df.ChargeL - self.df.ChargeR)/(self.df.ChargeR + self.df.ChargeL)
        #self.df["NormCountDiff"] = self.df.ChargeL/self.df.ChargeR

        plt.figure()
        plt.hist(self.df.NormCountDiff, bins=50)
        plt.title("Normalized Count Differences \n Channel Pair {},{}  Location {}mm".format(LChannel, RChannel, Location))
        plt.ylim(0,400)
        plt.xlim(-.5,.5)
        T_y = generateTextY(4,400)
        T_x = -.45
        plt.text(T_x,T_y[0],"mean: "+format(self.df.NormCountDiff.mean(),"5.4f"))
        plt.text(T_x,T_y[1],"std: "+format(self.df.NormCountDiff.std(),"7.4f"))
        plt.text(T_x,T_y[2],"FWHM: "+format(2.3548*self.df.NormCountDiff.std(),"7.4f"))
        plt.text(T_x,T_y[3],"count: "+str(self.df.NormCountDiff.count()))
        plt.savefig("Count_Diff, Pair {}, {}, {} mm".format(LChannel, RChannel, Location))
        self.DOI_mean = self.df.NormCountDiff.mean()
        self.DOI_FWHM = 2.3548*self.df.NormCountDiff.std()
        if int(LChannel) == 64 and int(RChannel) == 71:
            Hist = pd.DataFrame(columns = ['SigDiff'])
            Hist["SigDiff"] = self.df.NormCountDiff
            # Hist.to_csv(dir+"SigDiff{}.csv".format(Location))

        if display == True:
            plt.show()
            plt.clf()
            plt.close()
        else:
            plt.clf()
            plt.close()

    def PlotChargeDist(self, display=False):
        global LChannel
        global RChannel
        global Location

        plt.figure()
        plt.scatter(self.df.ChargeL, self.df.ChargeR,s=2)
        plt.title("Left Charge vs Right Charge \n Channel Pair {},{}  Location {}mm".format(LChannel, RChannel, Location))
        plt.xlim(8,25)
        plt.ylim(8,25)
        T_y = generateTextY(5,25)
        T_x = 8.01
        T_x1 = 14.5
        plt.text(T_x, T_y[0], "Left Charge Statistics")
        plt.text(T_x,T_y[1],"mean: "+format(self.df.ChargeL.mean(),"5.2f"))
        plt.text(T_x,T_y[2],"std: "+format(self.df.ChargeL.std(),"7.4f"))
        plt.text(T_x,T_y[3],"FWHM: "+format(2.3548*self.df.ChargeL.std(),"7.4f"))
        plt.text(T_x,T_y[4],"count: "+str(self.df.ChargeL.count()))
        plt.text(T_x1, T_y[0], "Right Charge Statistics")
        plt.text(T_x1,T_y[1],"mean: "+format(self.df.ChargeR.mean(),"5.2f"))
        plt.text(T_x1,T_y[2],"std: "+format(self.df.ChargeR.std(),"7.4f"))
        plt.text(T_x1,T_y[3],"FWHM: "+format(2.3548*self.df.ChargeR.std(),"7.4f"))
        plt.text(T_x1,T_y[4],"count: "+str(self.df.ChargeR.count()))
        plt.savefig("Charge Scatter, Pair {}, {}, {} mm DOI".format(LChannel, RChannel, Location))

        if int(LChannel) == 64 and int(RChannel) == 71:
            Hist = pd.DataFrame(columns = ['ChargeL', "ChargeR"])
            Hist["ChargeL"] = self.df.ChargeL
            Hist["ChargeR"] = self.df.ChargeR
            Hist.to_csv(dir+"QDiff{}.csv".format(Location))

        if display == True:
            plt.show()
            plt.clf()
            plt.close()
        else:
            plt.clf()
            plt.close()

    def ChargeNormalization(self, dataframe):
        global LChannel
        global RChannel
        global Location
        if Location == 2:
            dataframe = dataframe.append({"LChannel":LChannel, "RChannel":RChannel, "Normalization":(self.EFitParam[0][1]/self.EFitParam[1][1])}, ignore_index=True)
            return dataframe
        else:
            return dataframe



def gauss(x,A,mu,sigma):
    y = A*np.exp(-(x - mu)**2 / (2 * sigma**2))
    return y

def generateTextY(num,max):
    y = [max - 0.04*(i+1)*max for i in range(0,num+1)]
    return y

def toGeo(x):
    y = 8*indices.get(x)[0] + indices.get(x)[1]
    return y


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

#Takes in trip coinc files and colects DOI information by channel ID
run_num = 39
Location = 2
dir = ""
df = pd.read_csv(dir+"run6_4_21__coinc.txt", sep="\t", header=None, low_memory=False)
df.columns = ["A", "B", "TimeL", "ChargeL", "AbsChannelIDL", "C", "D", "TimeRef",
              "ChargeRef", "AbsChannelIDRef", "E", "F", "TimeR", "ChargeR", "AbsChannelIDR"]

df["ChannelIDR"] = df["AbsChannelIDR"] % 128
df["ChannelIDL"] = df["AbsChannelIDL"] % 128
df["ChannelIDRef"] = df["AbsChannelIDRef"] % 128

# Find the geometric ID of each channel.
df["GeoChannelIDL"] = df["ChannelIDL"].apply(toGeo)
df["GeoChannelIDR"] = df["ChannelIDR"].apply(toGeo)
df["GeoChannelIDRef"] = df["ChannelIDRef"].apply(toGeo)

#make sure the channel IDs line up with something reasonable
pairs = df.groupby(["GeoChannelIDL", "GeoChannelIDR", "GeoChannelIDRef"]).size().reset_index(name="Counts").sort_values(by='Counts', ascending=False).reset_index()
# print(pairs)

#Normalizing the locations of the photopeaks for the two SiPMs attached to the DOI crystal. The dataset taken in the middle of the crystal determines the factor the left SiPM is scaled by
if Location == 10:
    NFactors = pd.DataFrame(columns=["LChannel", "RChannel", "Normalization"])
else:
    NFactors = pd.read_csv(dir+"NFactors.csv", sep=",", header=0, low_memory=False)

# ref_ch = df[df["GeoChannelIDRef"] == 12]
# plt.hist(ref_ch.ChargeRef, bins=np.linspace(0,29,100))
# plt.title("Charge Spectra of most active reference crystal")
# plt.show()
# Find the ASIC of each channel.
df["ChipIDL"] = df["AbsChannelIDL"] // 64 % 64
df["ChipIDR"] = df["AbsChannelIDR"] // 64 % 64
df["ChipIDRef"] = df["AbsChannelIDRef"] // 64 % 64

# Only consider events for which DR SiPM 1 and DR SiPM 2
# correspond to the same crystal.
df = df[((df["GeoChannelIDL"] % 8) + (df["GeoChannelIDR"] % 8) == 7) &
        (df["GeoChannelIDL"] // 8 == df["GeoChannelIDR"] // 8)]


# # Discard events which are not really triple coincidences
# # across three modules.
df = df[(df["ChipIDRef"] != df["ChipIDL"]) & (df["ChipIDRef"] != df["ChipIDL"])]

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('binary')
cmap = truncate_colormap(cmap, 0.2, 0.8)

#DEBUG GRAPHS SHOWING COUNTS BY CHANNEL ID
# On a heatmap
ref_count = df.GeoChannelIDRef.value_counts().reset_index()
ref_count.columns = ["Geo_Channel", "Counts"]
for i in range(0,128):
    run = ref_count[ref_count["Geo_Channel"] == i]
    if run.Counts.count() == 0:
        ref_count = ref_count.append({"Geo_Channel":i, "Counts":0},ignore_index=True)

ref_count = ref_count.sort_values(by="Geo_Channel")
CTRLT = ref_count[ref_count["Geo_Channel"] < 64]
CTRRT = ref_count[ref_count["Geo_Channel"] >= 64]
CTRLT = np.reshape(np.array(CTRLT.Counts), (8,8))
CTRRT = np.reshape(np.array(CTRRT.Counts), (8,8))


plt.subplot(121)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        plt.text(j, i, format(int(CTRLT[i][j]),"3"), ha="center", va="center", color="b",size='large')
plt.title("Number of Coincidences in Reference Array")
plt.imshow(CTRLT,cmap=cmap,vmin=0, vmax=500)

plt.subplot(122)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        plt.text(j, i, format(int(CTRRT[i][j]),"3"), ha="center", va="center", color="b",size='large')
plt.imshow(CTRRT,cmap=cmap,vmin=0, vmax=500)
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
cax = plt.axes([0.85, 0.1, 0.025, 0.8])
plt.colorbar(cax=cax)
plt.show()

# df = df[(df["GeoChannelIDRef"] <= 64) & (df.GeoChannelIDRef % 8 == 4)]
# On a heatmap
ref_count = df.GeoChannelIDL.value_counts().reset_index()
ref_count.columns = ["Geo_Channel", "Counts"]
for i in range(0,128):
    run = ref_count[ref_count["Geo_Channel"] == i]
    if run.Counts.count() == 0:
        ref_count = ref_count.append({"Geo_Channel":i, "Counts":0},ignore_index=True)

ref_count = ref_count.sort_values(by="Geo_Channel")
CTRLT = ref_count[ref_count["Geo_Channel"] < 64]
CTRRT = ref_count[ref_count["Geo_Channel"] >= 64]
CTRLT = np.reshape(np.array(CTRLT.Counts), (8,8))
CTRRT = np.reshape(np.array(CTRRT.Counts), (8,8))


plt.subplot(121)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        plt.text(j, i, format(int(CTRLT[i][j]),"3"), ha="center", va="center", color="b",size='large')
plt.title("Number of Coincidences in Left DOI Array")
plt.imshow(CTRLT,cmap=cmap,vmin=0, vmax=500)

plt.subplot(122)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        plt.text(j, i, format(int(CTRRT[i][j]),"3"), ha="center", va="center", color="b",size='large')
plt.imshow(CTRRT,cmap=cmap,vmin=0, vmax=500)
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
cax = plt.axes([0.85, 0.1, 0.025, 0.8])
plt.colorbar(cax=cax)
plt.show()

# On a heatmap
ref_count = df.GeoChannelIDR.value_counts().reset_index()
ref_count.columns = ["Geo_Channel", "Counts"]
for i in range(0,128):
    run = ref_count[ref_count["Geo_Channel"] == i]
    if run.Counts.count() == 0:
        ref_count = ref_count.append({"Geo_Channel":i, "Counts":0},ignore_index=True)

ref_count = ref_count.sort_values(by="Geo_Channel")
CTRLT = ref_count[ref_count["Geo_Channel"] < 64]
CTRRT = ref_count[ref_count["Geo_Channel"] >= 64]
CTRLT = np.reshape(np.array(CTRLT.Counts), (8,8))
CTRRT = np.reshape(np.array(CTRRT.Counts), (8,8))


plt.subplot(121)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        plt.text(j, i, format(int(CTRLT[i][j]),"3"), ha="center", va="center", color="b",size='large')
plt.title("Number of Coincidences in Right DOI Array")
plt.imshow(CTRLT,cmap=cmap,vmin=0, vmax=500)

plt.subplot(122)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        plt.text(j, i, format(int(CTRRT[i][j]),"3"), ha="center", va="center", color="b",size='large')
plt.imshow(CTRRT,cmap=cmap,vmin=0, vmax=500)
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
cax = plt.axes([0.85, 0.1, 0.025, 0.8])
plt.colorbar(cax=cax)
plt.show()

# Main Body Analysis code
pairs = df.groupby(["GeoChannelIDL", "GeoChannelIDR"]).size().reset_index(name="Counts").sort_values(by='Counts', ascending=False).reset_index()
DATA = pd.DataFrame(columns=["GeoChannelIDL", "GeoChannelIDR", "DOI_mean", "DOI_FWHM", "DOI_Loc"])
# for j in tqdm(range(pairs.GeoChannelIDL.count())):
for j in tqdm(range(10)):
     test = df[df["GeoChannelIDL"] == pairs.GeoChannelIDL[j]]
     test = test[test["GeoChannelIDR"] == pairs.GeoChannelIDR[j]].reset_index()

     # # test["ModGeoChannelIDRef"] = test["GeoChannelIDRef"] % 8
     # # plt.hist(test["ModGeoChannelIDRef"], bins=np.linspace(-0.5,7.5,9))
     # # plt.show()
     # test = test[(test["GeoChannelIDRef"] <= 64)]#
     LChannel = pairs.GeoChannelIDL[j]
     RChannel = pairs.GeoChannelIDR[j]
     ChannelPair = Channel_Pair(test)
     ChannelPair.findPhotopeak()
     ChannelPair.SetEnergyParameters(ChannelPair.PPFit())
     NFactors = ChannelPair.ChargeNormalization(NFactors)
     test_NFactor = NFactors[NFactors["LChannel"] == LChannel]
     test_NFactor = test_NFactor[test_NFactor["RChannel"] == RChannel].reset_index()
     # ChannelPair.df["ChargeL"] = ChannelPair.df["ChargeL"]*test_NFactor["Normalization"][0]
     ChannelPair.findPhotopeak()
     ChannelPair.SetEnergyParameters(ChannelPair.PPFit())
     ChannelPair.PlotEnergySpectrum(display=True)
     if ChannelPair.EnergyCut():
         ChannelPair.CutOnEnergy(normed=True)
         if ChannelPair.PPCut():
             ChannelPair.SetCTRParameters(ChannelPair.CTRFit())
             ChannelPair.PlotCTR(display=True)
     ChannelPair.PlotCountsDiff(display=True)
     ChannelPair.PlotChargeDist(display=True)
     DATA = DATA.append({"GeoChannelIDL":LChannel, "GeoChannelIDR":RChannel, "DOI_mean":ChannelPair.DOI_mean,"DOI_FWHM":ChannelPair.DOI_FWHM, "DOI_Loc":Location}, ignore_index=True)

DATA.to_csv(dir+"DOI{}v2_DATA.csv".format(Location))
NFactors.to_csv(dir+"NFactorsv2.csv")




#Currently Unused Code
# #Impose limit on referecne column and DR channel
# # df_filt = df[(df.ChannelIDRef % 8 == 4) & ((df.ChannelIDRef // 8).isin([0,1,2,3,4,5,6,7,8,9,10,11]) )]
# df_filt = df[(df.ChannelIDRef % 8 == 4)  ]
# df_filt = df_filt[ df_filt.ChannelIDL == 5 ]
#
# df_filt = df_filt[ (df_filt.ChargeL+df_filt.ChargeR > 26) & (df_filt.ChargeRef > 20) ]
#
# fig, ax = plt.subplots()
# vals, bins = ax.hist( df_filt.NormCountDiff, bins = np.linspace(-.3, +.15, 40) )
#
# bins = [x - (bins[1]-bins[0])/2 for x in bins[1:]]
# opt, cov = curve_fit(gauss, bins, vals, p0 = [df_filt.NormCountDiff.median(), df_filt.NormCountDiff.std(), max(vals)])
#
# ax.plot(np.linspace(min(bins), max(bins), 100 ), gauss(np.linspace(min(bins), max(bins), 100), opt[0], opt[1], opt[2]),color = "k")
#
# ax.set_xlabel("Normalized Charge differences (non-linearized)", fontsize = 14)
# ax.set_ylabel("Number of Detections", fontsize = 14)
#
# # xdata = [-12, -6, 6, 12]
# # ydata = [-.160055, -.150558, -.114187, -.157960]




# #Plots all the DOI information by Channel Pair from the DOI RUNS
# dir = "~/path/Geant4/Expirimental_Data/"
# df = pd.read_csv(dir+"DOI2_DATA.csv", sep=",", header=0, low_memory=False)
# df1 = pd.read_csv(dir+"DOI6_DATA.csv", sep=",", header=0, low_memory=False)
# df2 = pd.read_csv(dir+"DOI10_DATA.csv", sep=",", header=0, low_memory=False)
# df3 = pd.read_csv(dir+"DOI14_DATA.csv", sep=",", header=0, low_memory=False)
# df4 = pd.read_csv(dir+"DOI18_DATA.csv", sep=",", header=0, low_memory=False)
#
# for i in range(df.GeoChannelIDL.count()):
#     ChannelIDL = df["GeoChannelIDL"][i]
#     ChannelIDR = df["GeoChannelIDR"][i]
#     DOI2 = df[df["GeoChannelIDL"] == ChannelIDL]
#     DOI2 = DOI2[DOI2["GeoChannelIDR"] == ChannelIDR].reset_index()
#
#     DOI6 = df1[df1["GeoChannelIDL"] == ChannelIDL]
#     DOI6 = DOI6[DOI6["GeoChannelIDR"] == ChannelIDR].reset_index()
#
#     DOI10 = df2[df2["GeoChannelIDL"] == ChannelIDL]
#     DOI10 = DOI10[DOI10["GeoChannelIDR"] == ChannelIDR].reset_index()
#
#     DOI14 = df3[df3["GeoChannelIDL"] == ChannelIDL]
#     DOI14 = DOI14[DOI14["GeoChannelIDR"] == ChannelIDR].reset_index()
#
#     DOI18 = df4[df4["GeoChannelIDL"] == ChannelIDL]
#     DOI18 = DOI18[DOI18["GeoChannelIDR"] == ChannelIDR].reset_index()
#
#     dfs = [DOI2, DOI6, DOI10, DOI14, DOI18]
#     DOIs = [channel.DOI_mean[0] for channel in dfs]
#     Errors = [0.5*(channel.DOI_FWHM[0]) for channel in dfs]
#     locs = [channel.DOI_Loc[0] for channel in dfs]
#     plt.errorbar(locs, DOIs,Errors,[1.5,1.5,1.5,1.5,1.5], ecolor='k',capsize=3)
#     plt.title("Channel Signal Difference vs DOI \n {}, {}".format(ChannelIDL, ChannelIDR))
#     plt.xlabel("DOI in mm")
#     plt.ylabel("Normalized Signal Difference")
#     plt.show()

# #Signal Difference histograms for a single channel pair
# dir = "~/path/Geant4/Expirimental_Data/"
# df = pd.read_csv(dir+"SigDiff2.csv", sep=",", header=0, low_memory=False)
# df1 = pd.read_csv(dir+"SigDiff6.csv", sep=",", header=0, low_memory=False)
# df2 = pd.read_csv(dir+"SigDiff10.csv", sep=",", header=0, low_memory=False)
# df3 = pd.read_csv(dir+"SigDiff14.csv", sep=",", header=0, low_memory=False)
# df4 = pd.read_csv(dir+"SigDiff18.csv", sep=",", header=0, low_memory=False)
#
# data = [df.SigDiff.tolist(),df1.SigDiff.tolist(),df2.SigDiff.tolist(),df3.SigDiff.tolist(),df4.SigDiff.tolist()]
# colors = ['r', 'b', 'g', 'c','m']
# labels = ["2mm", "6mm", "10mm", "14mm", "18mm"]
# plt.hist(data[0], bins=np.linspace(0.85,1.25,80), ec=colors[0], label=labels[0], histtype='step')
# plt.hist(data[1], bins=np.linspace(0.85,1.25,80),ec=colors[1], label=labels[1], histtype='step')
# plt.hist(data[2], bins=np.linspace(0.85,1.25,80),ec=colors[2], label=labels[2], histtype='step')
# plt.hist(data[3], bins=np.linspace(0.85,1.25,80),ec=colors[3], label=labels[3], histtype='step')
# plt.hist(data[4], bins=np.linspace(0.85,1.25,80),ec=colors[4], label=labels[4], histtype='step')
# plt.title("Normalized Signal Difference at various DOIs")
# plt.xlabel("Normalized Signal Difference")
# plt.ylabel("Counts")
# plt.legend()
# plt.show()

# #Charge Difference histograms for a single channel pair
# dir = "~/path/Geant4/Expirimental_Data/"
# df = pd.read_csv(dir+"QDiff2.csv", sep=",", header=0, low_memory=False)
# df1 = pd.read_csv(dir+"QDiff6.csv", sep=",", header=0, low_memory=False)
# df2 = pd.read_csv(dir+"QDiff10.csv", sep=",", header=0, low_memory=False)
# df3 = pd.read_csv(dir+"QDiff14.csv", sep=",", header=0, low_memory=False)
# df4 = pd.read_csv(dir+"QDiff18.csv", sep=",", header=0, low_memory=False)
#
# dataR = [df.ChargeR.tolist(),df1.ChargeR.tolist(),df2.ChargeR.tolist(),df3.ChargeR.tolist(),df4.ChargeR.tolist()]
# dataL = [df.ChargeL.tolist(),df1.ChargeL.tolist(),df2.ChargeL.tolist(),df3.ChargeL.tolist(),df4.ChargeL.tolist()]
# colors = ['r', 'b', 'g', 'c','m']
# labels = ["2mm", "6mm", "10mm", "14mm", "18mm"]
# plt.scatter(dataL[0], dataR[0], color=colors[0], label=labels[0])
# plt.scatter(dataL[1],dataR[1], color=colors[1], label=labels[1])
# plt.scatter(dataL[2],dataR[2], color=colors[2], label=labels[2])
# plt.scatter(dataL[3],dataR[3], color=colors[3], label=labels[3])
# plt.scatter(dataL[4],dataR[4], color=colors[4], label=labels[4])
# plt.title("Left Charge vs Right Charge at various DOIs")
# plt.xlabel("Left Charge")
# plt.ylabel("Right Charge")
# plt.xlim(12,20)
# plt.ylim(12,18)
# plt.legend()
# plt.show()

# #Time Difference histograms for a single channel pair
# dir = "~/path/Geant4/Expirimental_Data/"
# df = pd.read_csv(dir+"TimeDiff2.csv", sep=",", header=0, low_memory=False)
# df1 = pd.read_csv(dir+"TimeDiff6.csv", sep=",", header=0, low_memory=False)
# df2 = pd.read_csv(dir+"TimeDiff10.csv", sep=",", header=0, low_memory=False)
# df3 = pd.read_csv(dir+"TimeDiff14.csv", sep=",", header=0, low_memory=False)
# df4 = pd.read_csv(dir+"TimeDiff18.csv", sep=",", header=0, low_memory=False)
#
# data = [df.TimeDiff.tolist(),df1.TimeDiff.tolist(),df2.TimeDiff.tolist(),df3.TimeDiff.tolist(),df4.TimeDiff.tolist()]
# colors = ['r', 'b', 'g', 'c','m']
# labels = ["2mm", "6mm", "10mm", "14mm", "18mm"]
# plt.hist(data[0], bins=np.linspace(0.6,1.6,25), ec=colors[0], label=labels[0], histtype='step')
# plt.hist(data[1], bins=np.linspace(0.6,1.6,25),ec=colors[1], label=labels[1], histtype='step')
# plt.hist(data[2], bins=np.linspace(0.6,1.6,25),ec=colors[2], label=labels[2], histtype='step')
# plt.hist(data[3], bins=np.linspace(0.6,1.6,25),ec=colors[3], label=labels[3], histtype='step')
# plt.hist(data[4], bins=np.linspace(0.6,1.6,25),ec=colors[4], label=labels[4], histtype='step')
# plt.title("Time Difference at various DOIs")
# plt.xlabel("Time diff in ns")
# plt.ylabel("Counts")
# plt.legend()
# plt.show()
