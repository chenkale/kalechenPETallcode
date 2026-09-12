#Will Matava, Kyle Klein

import glob
import os
from MiniPET_Vis import *
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import axes3d, Axes3D
import math
import matplotlib
import statistics as stat
#from tqdm import tqdm
#import time
# UNCOMMENT FOR FIRST PART OF PROGRAM
#matplotlib.use('Agg')


pd.set_option("display.max_columns",100)

class Channel_Pair:
    def __init__(self, ChannelPairData):
        self.df = ChannelPairData
        valuesR,binsR,presentersR = plt.hist(self.df.ChargeR, bins = np.linspace(0,45,160), ec='r',label='Reference PCB', fill=False)
        valuesL,binsL,presentersL = plt.hist(self.df.ChargeL, bins = np.linspace(0,45,160), ec='b',label='PCB of Interest', fill=False)
        bin_centers = [(binsR[q] + binsR[q+1])/2 for q in range(0,len(binsR)-1)]
        self.EHistR_y = valuesR.tolist()
        self.EHistR_x = bin_centers
        self.EHistL_y = valuesL.tolist()
        self.EHistL_x = bin_centers
        self.PhotopeakR = 0
        self.PhotopeakL = 0
        self.EResR = 0
        self.EResL = 0
        self.PP_Count = 0
        self.EFitParam = 0
        self.CTRFitParam = 0
        plt.clf()
        plt.close()

    def findPhotopeak(self):
        #find the peak within 3 bins starting from the RHS of the graph
        for i in range(1, len(self.EHistR_y)-4):
            if self.EHistR_y[-1*i] > self.EHistR_y[-1*(i+1)] and self.EHistR_y[-1*i] > self.EHistR_y[-1*(i+2)] and self.EHistR_y[-1*i] > self.EHistR_y[-1*(i+3)] and self.EHistR_y[-1*(i+1)] != 0 and self.EHistR_y[-1*i] >= 20:
                self.PhotopeakR = self.EHistR_y.index(self.EHistR_y[-1*i],-1*i-1)
                #print("   - Right photopeak estimate found!!!")
                break
        for i in range(1, len(self.EHistL_y)-4):
            if self.EHistL_y[-1*i] > self.EHistL_y[-1*(i+1)] and self.EHistL_y[-1*i] > self.EHistL_y[-1*(i+2)] and self.EHistL_y[-1*i] > self.EHistL_y[-1*(i+3)] and self.EHistL_y[-1*(i+1)] != 0 and self.EHistL_y[-1*i] >= 20:
                self.PhotopeakL = len(self.EHistL_y)-1*i
                #print("   - Left photopeak estimate found!!!")
                break

    def getQINTs(self):
        # numpy requires the arguments y before x
        INTR = np.trapz(self.EHistR_y,self.EHistR_x)
        INTL = np.trapz(self.EHistL_y,self.EHistL_x)
        return [INTR, INTL]

    def PPFit(self):
        ExR_p = [0,1,1]
        ExL_p = [0,1,1]

        fitR_x = self.EHistR_x[self.PhotopeakR-7:]
        fitL_x = self.EHistL_x[self.PhotopeakL-7:]
        fitR_y = self.EHistR_y[self.PhotopeakR-7:]
        fitL_y = self.EHistL_y[self.PhotopeakL-7:]
        # plt.scatter(fitL_x, fitL_y)
        # plt.show()

        #fitting said data and recording photpeak energy cuts
        if len(fitR_y) != 0 and len(fitR_x) != 0:
            try:
                ExR_p, ExR_co = curve_fit(gauss, fitR_x, fitR_y, p0=[max(fitR_y), fitR_x[fitR_y.index(max(fitR_y))], 0.5], bounds=[[max(fitR_y)-50, fitR_x[fitR_y.index(max(fitR_y))]-2, 0.001],[max(fitR_y)+50, fitR_x[fitR_y.index(max(fitR_y))]+2, 1]])
            except RuntimeError:
                print("Error - right photopeak curve fit failed!!!")

        if len(fitL_x) != 0 and len(fitL_y) != 0:
            try:
                ExL_p, ExL_co = curve_fit(gauss, fitL_x, fitL_y, p0=[max(fitL_y), fitL_x[fitL_y.index(max(fitL_y))], 0.5], bounds=[[max(fitL_y)-50, fitL_x[fitL_y.index(max(fitL_y))]-2, 0.001],[max(fitL_y)+50, fitL_x[fitL_y.index(max(fitL_y))]+2, 1]])
            except RuntimeError:
                print("Error - left photopeak curve fit failed!!!")

        return [ExR_p, ExL_p]

    def SetEnergyParameters(self, run_param):
        self.EFitParam = run_param
        self.EResR = abs(100*(2.3548*abs(self.EFitParam[0][2]))/self.EFitParam[0][1])
        self.EResL = abs(100*(2.3548*abs(self.EFitParam[1][2]))/self.EFitParam[1][1])

    def CutOnEnergy(self):
        R_cut = self.EFitParam[0][1] - 2.5*self.EFitParam[0][2]
        L_cut = self.EFitParam[1][1] - 2.5*self.EFitParam[1][2]
        # R_cut = self.EFitParam[0][1] - 1.0*self.EFitParam[0][2]
        # L_cut = self.EFitParam[1][1] - 1.0*self.EFitParam[1][2]
        self.df = self.df[self.df["ChargeR"] >= R_cut]
        self.df = self.df[self.df["ChargeL"] >= L_cut].reset_index(drop=True)
        self.df["Time_diff"] = 10**-3 * (self.df["TimeL"] - self.df["TimeR"])

    def EnergyCut(self):
        if self.EFitParam[0][1] > 10.5 and self.EFitParam[1][1] > 10.5: #and self.EFitParam[0][2] < 200 and self.EFitParam[1][2] < 200:
            return True
        else:
            return False

    def PPCut(self):
        if self.df["Time_diff"].count() < 50:
        #if self.df["Time_diff"].count() < 100:
        #if self.df["Time_diff"].count() < 300:
            return False
        else:
            self.PP_Count = self.df["Time_diff"].count()
            return True

    def SetCTRParameters(self, CTR_param):
        self.CTRFitParam = CTR_param

    def CTRFit(self):
        global LChannel
        global RChannel
        TRes5_UN = self.df["Time_diff"]

        # time_cutU = self.df[self.df["Time_diff"] < -2]
        # time_cutL = self.df[self.df["Time_diff"] > 2]
        # plt.hist(time_cutU.ChargeR, bins = np.linspace(0,max(time_cutU.ChargeR),int((160*max(time_cutU.ChargeR))/45)), ec='r', fill=False)
        # plt.hist(time_cutU.ChargeL, bins = np.linspace(0,max(time_cutU.ChargeL),int((160*max(time_cutU.ChargeL))/45)), ec='b', fill=False)
        # plt.title("Energy Spectrum, Pair {}, {}".format(toGeoChannelID(LChannel), toGeoChannelID(RChannel)))
        # plt.show()

        x1_p = [0,1,1]
        low = TRes5_UN.min()
        high = TRes5_UN.max()
        #changed to always keep the bin size constant to keep comparison possible
        #num_bins = 600
        num_bins = 200
        bins = np.linspace(low,high,num_bins)
        #x = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),num_bins-1) # WHAT IS THIS?? WHY IS IT HERE???
        y1 = TRes5_UN.value_counts(bins=bins,sort=False).tolist()
        peak = y1.index(max(y1))
        bin_centers = [(bins[i] + bins[i+1])/2 for i in range(0,len(bins)-1)]
        fit_x = bin_centers[max(0,peak-10):peak+11]
        fit_y = y1[max(0,peak-10):peak+11]
        if len(fit_x) != 0 and len(fit_y) != 0:
            try:
                x1_p, x1_co = curve_fit(gauss,fit_x,fit_y, p0=[max(fit_y),stat.mean(fit_x),0.5*stat.pstdev(fit_x)], bounds=[[max(fit_y)-10,stat.mean(fit_x)-0.1,0],[max(fit_y)+10,stat.mean(fit_x)+0.1,abs(stat.pstdev(fit_x))]])
            except RuntimeError:
                print("Error - CTR curve fit failed!!!")

        return x1_p


    def PlotEnergySpectrum(self, display=False):
        global LChannel
        global RChannel
        global PCB_num
        # plotting every energy spectrum
        x_max = max(self.df.ChargeR.tolist() + self.df.ChargeL.tolist())
        x_f = np.linspace(0,x_max,300)
        plt.hist(self.df.ChargeR, bins = np.linspace(0,max(self.df.ChargeR),int((160*max(self.df.ChargeR))/45)), ec='r',label='Channel ID {}'.format(toGeoChannelID(LChannel)), fill=False)
        plt.hist(self.df.ChargeL, bins = np.linspace(0,max(self.df.ChargeL),int((160*max(self.df.ChargeL))/45)), ec='b',label='Channel ID {}'.format(toGeoChannelID(RChannel)), fill=False)
        plt.plot(x_f, gauss(x_f,*self.EFitParam[0]), color='k')
        plt.plot(x_f, gauss(x_f,*self.EFitParam[1]), color='k')
        t_y = generateTextY(10,400)#400
        t_x = 1
        plt.text(t_x,t_y[0], 'Geo Channel ID {}'.format(toGeoChannelID(LChannel)))
        plt.text(t_x,t_y[1], "mean: "+format(self.EFitParam[0][1], "5.2f"))
        plt.text(t_x,t_y[2], "std: "+format(abs(self.EFitParam[0][2]), "5.2f"))
        plt.text(t_x,t_y[3], "count: "+format(ChannelPair.df.ChargeR.count(), "5"))
        plt.text(t_x,t_y[4], "ERes: "+format(self.EResR, "5.2f")+"%")
        plt.text(t_x,t_y[5], 'Geo Channel ID {}'.format(toGeoChannelID(RChannel)))
        plt.text(t_x,t_y[6], "mean: "+format(self.EFitParam[1][1], "5.2f"))
        plt.text(t_x,t_y[7], "std: "+format(abs(self.EFitParam[1][2]), "5.2f"))
        plt.text(t_x,t_y[8], "count: "+format(ChannelPair.df.ChargeL.count(), "5"))
        plt.text(t_x,t_y[9], "ERes: "+format(self.EResL, "5.2f")+"%")
        plt.ylim(0,400)
        #plt.title("Energy Spectrum, Pair {}, {}".format(toGeoChannelID(LChannel), toGeoChannelID(RChannel)))
        plt.title("Energy Spectrum, Pair {}, {}".format(LChannel, RChannel))
        plt.legend()
        #plt.savefig("Energy Spectrum, Pair {}, {}".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
        plt.savefig("Energy Spectrum, Pair {}, {}".format(int(LChannel), int(RChannel)))
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
        global PCB_num
        TRes5_UN = self.df["Time_diff"]
        fig1 = plt.figure()
        #ax1 = fig1.add_subplot(111,xlim=(0,9), ylim=(0,10),xlabel="Time Difference in ns",ylabel="Counts",title="CTR: Channel Pair {}, {}".format(toGeoChannelID(LChannel),toGeoChannelID(RChannel)))
        ax1 = fig1.add_subplot(111,xlim=(0,10), ylim=(0,20),xlabel="Time Difference in ns",ylabel="Counts",title="CTR: Channel Pair {}, {}".format(LChannel,RChannel))
        low = TRes5_UN.min()
        high = TRes5_UN.max()
        #num_bins = 600
        num_bins = 200
        bins = np.linspace(low,high,num_bins)
        x_p = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),400)
        T_x = 0.5#-1.9
        T_y = generateTextY(9,19.5)#100
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
        #plt.savefig("CTR, Pair {}, {}".format(toGeoChannelID(LChannel), toGeoChannelID(RChannel)))
        plt.savefig("CTR, Pair {}, {}".format(LChannel, RChannel))
        if display == True:
            plt.show()
            plt.clf()
            plt.close()
        else:
            plt.clf()
            plt.close()


#converts PETSys ID to geometric ID
def toGeo(x):
    y = 8*indices.get(x)[0] + indices.get(x)[1]
    return y


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

dir = ""

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

geo_channels = pd.DataFrame(columns = ["AbsChannelID","GeoChannelID"])
for i in range(128):
    slice = pd.DataFrame({"AbsChannelID": [i], "GeoChannelID": [toGeo(i)]})
    geo_channels = pd.concat([geo_channels, slice],ignore_index=True)

def toGeoChannelID(AbsChannelID):
    #old
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

#t = time.time()
chunks = pd.read_csv(dir+"OV_3.5_T1_10_Slave1Chip0.txt", sep="\t", header=None, names=["TimeL", "ChargeL", "ChannelIDL","TimeR", "ChargeR", "ChannelIDR"], usecols=[2,3,4,7,8,9],low_memory=False,chunksize=100000)
concat_df = pd.DataFrame()
#making a cut on energy that no photopeak is below, which just cuts comptons
for chunk in chunks:
    chunk = chunk.loc[chunk["ChargeL"] >= 10]
    chunk = chunk.loc[chunk["ChargeR"] >= 10].reset_index(drop=True)
    concat_df = pd.concat([concat_df, chunk],ignore_index=True)

del chunks

# #takes a little less than thirty seconds
# concat_df = pd.concat(chunks,ignore_index=True)

#making the almighty dataframe which will store everything we need, praise be
DATA = pd.DataFrame(columns=["ChannelIDL", "ChannelIDR", "PP_Count", "Channel_Eff","PeakL", "PeakR", "E_ResL", "E_ResR", "QINTL", "QINTR", "CTR", "E_Cut", "PP_Cut"])

pairs = concat_df.groupby(["ChannelIDL", "ChannelIDR"]).size().reset_index(name="Counts").sort_values(by='Counts', ascending=False).reset_index()
num_photons = 300
valid_channels = pairs[pairs["Counts"] >= num_photons].reset_index()
valid_channels.to_csv(dir+"reduced_List_mode.csv", sep="\t")
#print(time.time() - t)
#for j in tqdm(range(valid_channels.ChannelIDL.count())):
for j in range(valid_channels.ChannelIDL.count()):

    test = concat_df[concat_df["ChannelIDL"] == valid_channels.ChannelIDL[j]]
    test = test[test["ChannelIDR"] == valid_channels.ChannelIDR[j]].reset_index()
    RChannel = valid_channels.ChannelIDR[j]
    LChannel = valid_channels.ChannelIDL[j]

    ChannelPair = Channel_Pair(test)
    #print(" Looking for photopeaks in channel pair {} {}".format(LChannel,RChannel))
    ChannelPair.findPhotopeak() # Estimate the photopeak mean 
    ChannelPair.SetEnergyParameters(ChannelPair.PPFit()) # Do the photopeak fit
    QInts = ChannelPair.getQINTs()
    ChannelPair.PlotEnergySpectrum(display=False)

    #If the photopeak was found
    if ChannelPair.EnergyCut() == True:
        #print("     Photopeak is above threshold, cutting on events within 2.5 sigma of mean")
        ChannelPair.CutOnEnergy()
        if ChannelPair.PPCut() == True:
            #print("        Sufficient statistics (counts > 50) in photopeak")
            ChannelPair.SetCTRParameters(ChannelPair.CTRFit())  # Do the CTR fit
            ChannelPair.PlotCTR(display=False)
            line = pd.DataFrame({"ChannelIDL":[LChannel], "ChannelIDR":[RChannel], "PP_Count":[ChannelPair.PP_Count], "Channel_Eff": [ChannelPair.PP_Count/valid_channels.Counts[j]],"PeakL":[ChannelPair.EFitParam[1][1]], "PeakR":[ChannelPair.EFitParam[0][1]], "E_ResL":[ChannelPair.EResL], "E_ResR":[ChannelPair.EResR],"QINTL":[QInts[1]], "QINTR":[QInts[0]],"CTR":[abs(2.3548* ChannelPair.CTRFitParam[2])], "E_Cut":[1], "PP_Cut":[1]})
            DATA = pd.concat([DATA,line], ignore_index=True)
        elif ChannelPair.PPCut() == False:
            line = pd.DataFrame({"ChannelIDL":[LChannel], "ChannelIDR":[RChannel], "PP_Count":[ChannelPair.PP_Count], "Channel_Eff": [ChannelPair.PP_Count/valid_channels.Counts[j]],"PeakL":[ChannelPair.EFitParam[1][1]], "PeakR":[ChannelPair.EFitParam[0][1]], "E_ResL":[ChannelPair.EResL], "E_ResR":[ChannelPair.EResR],"QINTL":[QInts[1]], "QINTR":[QInts[0]],"CTR":[0], "E_Cut":[1], "PP_Cut":[0]})
            DATA = pd.concat([DATA,line], ignore_index=True)
    elif ChannelPair.EnergyCut() == False:
        line = pd.DataFrame({"ChannelIDL":[LChannel], "ChannelIDR":[RChannel], "PP_Count":[0], "Channel_Eff": [ChannelPair.PP_Count/valid_channels.Counts[j]],"PeakL":[ChannelPair.EFitParam[1][1]], "PeakR":[ChannelPair.EFitParam[0][1]], "E_ResL":[ChannelPair.EResL], "E_ResR":[ChannelPair.EResR],"QINTL":[QInts[1]], "QINTR":[QInts[0]],"CTR":[0], "E_Cut":[0], "PP_Cut":[0]})
        DATA = pd.concat([DATA,line], ignore_index=True)

DATA.to_csv(dir+"Output_OV_3.5_T1_10_Slave1Chip0.csv")
