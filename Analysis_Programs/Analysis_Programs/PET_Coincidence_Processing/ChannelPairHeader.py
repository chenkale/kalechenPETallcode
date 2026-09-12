#Will Matava, Kyle Klein
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import axes3d, Axes3D
import math
import matplotlib
import matplotlib.pyplot as plt
import statistics as stat
from tqdm import tqdm
import numpy as np
from numpy import genfromtxt
import time
from scipy.optimize import curve_fit
#plt.style.use('ChannelPairStyles.mplstyle')
# matplotlib.use('Agg')

def gauss(x,A,mu,sigma):
    # Simple Gaussian function to be used by the fitter
    y = A*np.exp(-(x - mu)**2 / (2 * sigma**2))
    return y

def generateTextY(num,max):
    # Function to help overlay useful text on histogram plots
    y = [max - 0.04*(i+1)*max for i in range(0,num+1)]
    return y

def HistogramCharge(ChannelPairData):
    # Takes in channel pair data frame, histograms ChargeL and ChargeR, and then returns values and bin centers
    ChargeData = ChannelPairData[:,[1,4]] # Only takes in columns for ChargeL and ChargeR
    # Can adjust the binning below to set a different range for the x-axis (energy range)
    valuesL, binsL = np.histogram(ChargeData[:,0], bins = np.linspace(0,45,160))
    valuesR, binsR = np.histogram(ChargeData[:,1], bins = np.linspace(0,45,160))
    # Binning will be same for both so just define bin centers once w.r.t. binsR
    bin_centers = [(binsR[q] + binsR[q+1])/2 for q in range(0,len(binsR)-1)]
    return valuesL, valuesR, bin_centers

def findPhotopeak(valuesL, valuesR, debug):
    # Find the photopeak within 3 bins starting from the RHS of the graph
    threshR = sum(valuesR) * 0.025
    for i in range(1, len(valuesR)-4):
        PhotopeakR = 0
        # If the -i_th bin is larger than ALL 3 bins to the left of it, set it as the photopeak
        if valuesR[-1*i] > valuesR[-1*(i+1)] and valuesR[-1*i] > valuesR[-1*(i+2)] and valuesR[-1*i] > valuesR[-1*(i+3)] and valuesR[-1*(i+1)] != 0 and valuesR[-1*i] >= threshR:
            PhotopeakR = len(valuesR)-1*i
            if debug:
                print("   - Right photopeak estimate found!!!")
            break
    # Repeat the algorithm for the left charge spectrum
    threshL = sum(valuesL) * 0.025
    for i in range(1, len(valuesL)-4):
        PhotopeakL = 0
        if valuesL[-1*i] > valuesL[-1*(i+1)] and valuesL[-1*i] > valuesL[-1*(i+2)] and valuesL[-1*i] > valuesL[-1*(i+3)] and valuesL[-1*(i+1)] != 0 and valuesL[-1*i] >= threshL:
            PhotopeakL = len(valuesL)-1*i
            if debug:
                print("   - Left photopeak estimate found!!!")
            break
    return PhotopeakL, PhotopeakR
    # def getQINTs(self):
    #     # numpy requires the arguments y before x
    #     INTR = np.trapz(valuesR,bin_centers)
    #     INTL = np.trapz(valuesL,bin_centers)
    #     return [INTR, INTL]


def fitPhotopeak(valuesL, valuesR, bin_centers, PhotopeakL, PhotopeakR):
    # function uses the photopeak location estimate to fit a Gaussian to the photopeaks
    # Arrays to hold the resultant fit parameters for the Gaussian: mean, sigma, and amplitude
    fitR_p = [0,1,1] # Initialized to amplitude=0, mean=1, sigma=1
    fitL_p = [0,1,1]
    # Set the range of x and y values to fit:
    #   from 7 bins to the left of the photopeak, up to the end of the spectra
    fitR_x = list(bin_centers[PhotopeakR-7:])
    fitL_x = list(bin_centers[PhotopeakL-7:])
    fitR_y = list(valuesR[PhotopeakR-7:])
    fitL_y = list(valuesL[PhotopeakL-7:])
        # plt.scatter(fitL_x, fitL_y)
        # plt.show()

    # Make sure fit ranges are non-empty, perform the fits, return parameters
    if len(fitR_y) != 0 and len(fitR_x) != 0:
        try:
            # curve_fit(function, x_data, y_data, p0=param_guesses, bounds=lower and upper param bounds)
            fitR_p, fitR_co = curve_fit(gauss, fitR_x, fitR_y, p0=[max(fitR_y), fitR_x[fitR_y.index(max(fitR_y))], 0.5], bounds=[[max(fitR_y)-50, fitR_x[fitR_y.index(max(fitR_y))]-2, 0],[max(fitR_y)+50, fitR_x[fitR_y.index(max(fitR_y))]+2, 1.5*abs(stat.pstdev(fitR_x))]])
        except RuntimeError:
            print("Error - right photopeak curve fit failed!!!")
    # Repeat the fit for the left charge spectrum
    if len(fitL_x) != 0 and len(fitL_y) != 0:
        try:
            fitL_p, fitL_co = curve_fit(gauss, fitL_x, fitL_y, p0=[max(fitL_y), fitL_x[fitL_y.index(max(fitL_y))], 0.5], bounds=[[max(fitL_y)-50, fitL_x[fitL_y.index(max(fitL_y))]-2, 0],[max(fitL_y)+50, fitL_x[fitL_y.index(max(fitL_y))]+2, 1.5*abs(stat.pstdev(fitL_x))]])
        except RuntimeError:
            print("Error - left photopeak curve fit failed!!!")
    # Return the fit parameter arrays (covariance arrays aren't needed)
    return [fitR_p, fitL_p] # First index denotes right(0) or left(1), second denotes params A(0), mu(1), sigma(2)

def PPThresholdCut(EFitParam):
    # Check that the photopeak locations are above 10.5 charge units
    if EFitParam[0][1] > 10.5 and EFitParam[1][1] > 10.5:
        return True
    else:
        return False

def CutOnEnergy(EFitParam, ChannelPairData):
    # Use the fitted mean and sigma to define an energy cut based on the photopeaks
    # Cut value is define to be 2.5 sigma below (to the left of) the mean
    R_cut = EFitParam[0][1] - 2.5*EFitParam[0][2]
    L_cut = EFitParam[1][1] - 2.5*EFitParam[1][2]
    # R_cut = EFitParam[0][1] - 1.0*EFitParam[0][2]
    # L_cut = EFitParam[1][1] - 1.0*EFitParam[1][2]
    R_cutu = EFitParam[0][1] + 2.5*EFitParam[0][2]
    L_cutu = EFitParam[1][1] + 2.5*EFitParam[1][2]
    # Only keep events with energies equal to or above the cut values
    ChannelPairData = ChannelPairData[ChannelPairData[:,4] >= R_cut]
    ChannelPairData = ChannelPairData[ChannelPairData[:,1] >= L_cut]
    ChannelPairData = ChannelPairData[ChannelPairData[:,4] <= R_cutu]
    ChannelPairData = ChannelPairData[ChannelPairData[:,1] <= L_cutu]
    return ChannelPairData

def PPOccupancyCut(PPCut_num, ChannelPairData, CTRFitParam):
    # # Check that the photopeaks have a sufficient number of events
    # if len(ChannelPairData) < PPCut_num:
    #     return False, 0
    # else:
    #     time_diff = np.subtract(ChannelPairData[:,0], ChannelPairData[:,3])
    #     time_diff.sort()
    #     STD = np.std(time_diff)
    #     for i in range(time_diff.size - 1):
    #         if (time_diff[i+1] - time_diff[i]) > 7*STD:
    #             break
    #         else:
    #             continue
    #     if i == time_diff.size - 2:
    #         #no outlier was found
    #         return True, time_diff
    #     else:
    #         time_diff = time_diff[i:]
    #         return True, time_diff
    # Check that the photopeaks have a sufficient number of events
    if len(ChannelPairData) < PPCut_num:
        return False, 0
    else:
        time_diff = np.subtract(ChannelPairData[:,0], ChannelPairData[:,3])
#        time_diff = np.sort(time_diff)
        Min = CTRFitParam[1] - 5*CTRFitParam[2]
        Max = CTRFitParam[1] + 5*CTRFitParam[2]

        time_diff.sort()
        STD = np.std(time_diff)
#        print(STD)
#        for i in range(time_diff.size - 1):
#            if (time_diff[i+1] - time_diff[i]) < 3*STD and (time_diff[i+2] - time_diff[i+1]) < 3*STD:
#                break
#            else:
#                continue
#        if i == time_diff.size - 2:
#            print("no left outlier was found")
#            return True, time_diff
#        for j in range(time_diff.size - 1):
#            if (time_diff[-1*j - 1] - time_diff[j]) < 3*STD and (time_diff[-1*j - 2] - time_diff[-1*j - 1]) < 3*STD:
#                break
#            else:
#                continue
#        if j == time_diff.size - 2:
#            time_diff = time_diff[i:]
##            print("no right outlier found")
#            return True, time_diff
#        else:
#            time_diff = time_diff[i: len(time_diff)- j]
##            print("left and right outlier found")
#            return True, time_diff
        if Max <= Min:
            print("Failed Input")
        else:
            time_diff = time_diff[time_diff >= Min]
            time_diff = time_diff[time_diff <= Max]
            return True, time_diff

def fitTimeDiff(time_diff):
    # Histogram the time differences for photopeak events and fit to Gaussian
    # Array to hold the resultant fit parameters for the Gaussian: mean, sigma, and amplitude
    fit_p = [0,1,1]
    # Set ranges for the data and histograms based on the max and min values
    low = time_diff.min()
    high = time_diff.max()
    # Calculate the bin numbers and sizes so they are constant to help with comparisons
    num_bins = int((high-low)/50)
    bins = np.linspace(low,high,num_bins)
    #sometimes the time differences are real close together, ya see
    if len(bins) < 20:
        bins = np.linspace(low, high, 20)
    # Create the time difference histogram
    values, bins = np.histogram(time_diff, bins=bins)
    # Prepare arguments for the fitter
    peak = np.argmax(values) # Estimate the peak to be the bin with the most counts
    bin_centers = bins + ((bins.max() - bins.min()) / (2*num_bins))
    bin_centers.resize((len(bin_centers) -1))
    # Set the x and y range for the fitter to be 10 bins around the peak
    fit_x = bin_centers[max(0,peak-10):peak+11]
    #fit_x = bin_centers
    # print(peak)
    # print(len(bin_centers))
    fit_y = values[max(0,peak-10):peak+11]
    #fit_y = values
    # Make sure fit ranges are non-empty, perform the fits, return parameters
    if len(fit_x) != 0 and len(fit_y) != 0:
        try:
            # curve_fit(function, x_data, y_data, p0=param_guesses, bounds=lower and upper param bounds)
            #print([[max(fit_y)-10,stat.mean(fit_x)-100,0],[max(fit_y)+10,stat.mean(fit_x)+100,1.5*abs(stat.pstdev(fit_x))]])
            #fit_p, fit_co = curve_fit(gauss, fit_x, fit_y, p0=[max(fit_y), stat.mean(fit_x), 0.5*stat.pstdev(fit_x)], bounds=[[max(fit_y)-10,stat.mean(fit_x)-100,0],[max(fit_y)+10,stat.mean(fit_x)+100,1.5*abs(stat.pstdev(fit_x))]])
            fit_p, fit_co = curve_fit(gauss, fit_x, fit_y, p0=[max(fit_y), stat.mean(fit_x), 0.5*stat.pstdev(fit_x)])
        except RuntimeError:
            print("Error - CTR curve fit failed!!!")
    return fit_p

def PlotEnergySpectrum(ChannelPairData, EFitParam, LChannel, RChannel, display=False):
    # plotting energy spectra for every channel pair and saving them all to .png files
    x_max = 60 # Hard coded value for the maximum x-range
    x_f = np.linspace(0,x_max,400) # Sets the binning: 300 bins from 0 to x_max
    # Plot the two histograms and their fitted functions
    fig = plt.figure(figsize=(14,9))
    valuesL, binsL, paramsL = plt.hist(ChannelPairData[:,4], bins = np.linspace(0,max(ChannelPairData[:,4]),int((160*max(ChannelPairData[:,4]))/45)), ec='r',label='Channel ID {}'.format(toGeoChannelID(LChannel)), fill=False)
    valuesR, binsR, paramsR = plt.hist(ChannelPairData[:,1], bins = np.linspace(0,max(ChannelPairData[:,1]),int((160*max(ChannelPairData[:,1]))/45)), ec='b',label='Channel ID {}'.format(toGeoChannelID(RChannel)), fill=False)
    plt.plot(x_f, gauss(x_f,*EFitParam[0]), color='k')
    plt.plot(x_f, gauss(x_f,*EFitParam[1]), color='k')
    # Code to help generate text overlays for plot info
    y_max = 1.05*max(max(valuesL), max(valuesR))
    t_y = generateTextY(10,y_max)
    t_x = 1
    plt.text(t_x,t_y[0], 'Geo Channel ID {}'.format(int(toGeoChannelID(LChannel))))
    plt.text(t_x,t_y[1], "mean: "+format(EFitParam[0][1], "5.2f"))
    plt.text(t_x,t_y[2], "std: "+format(abs(EFitParam[0][2]), "5.2f"))
    plt.text(t_x,t_y[3], "ERes: "+format(100*(2.355*abs(EFitParam[0][2]) / EFitParam[0][1]), "5.2f")+"%")
    plt.text(t_x,t_y[4], 'Geo Channel ID {}'.format(int(toGeoChannelID(RChannel))))
    plt.text(t_x,t_y[5], "mean: "+format(EFitParam[1][1], "5.2f"))
    plt.text(t_x,t_y[6], "std: "+format(abs(EFitParam[1][2]), "5.2f"))
    plt.text(t_x,t_y[7], "ERes: "+format(100*(2.355*abs(EFitParam[1][2]) / EFitParam[1][1]), "5.2f")+"%")
    plt.text(t_x,t_y[8], "count: "+format(len(ChannelPairData), "5"))
    # Set the y-axis limit
    plt.ylim(0,y_max)
    # Set the title and axis labels, draw the legend, and then save the figure to a .png
    plt.title("Energy Spectrum - Pair {}_{}".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
    plt.xlabel("Energy in DAQ units")
    plt.ylabel("Counts")
    plt.legend()
    plt.savefig("Energy Spectrum - Pair {}_{}".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
    if display:
        plt.show()
        plt.clf()
        plt.close()
    else:
        plt.clf()
        plt.close()

def PlotCTR(time_diff, CTRFitParam, LChannel, RChannel, display=False):
    # plotting time difference spectra for every channel pair and saving them all to .png files
    fig1 = plt.figure(figsize=(14,9))
    low = time_diff.min()
    high = time_diff.max()
    # Plot the histogram and set the title and axis labels
    ax1 = fig1.add_subplot(111,xlim=(low,high),xlabel="Time Difference in ps",ylabel="Counts",title="CTR -Channel Pair {}_{}".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
    # Calculate the bin numbers and sizes so they are constant to help with comparisons
    num_bins = int((high-low)/50)
    bins = np.linspace(low,high,num_bins)
    x_p = np.linspace(low+((high-low)/(2*num_bins)),high-((high-low)/(2*num_bins)),400)
    # Code to help generate text overlays for plot info
    T_x = 0.97*low + 0.03*high

    values, bins, params = ax1.hist(time_diff, bins=bins, color='b')
    ax1.plot(x_p,gauss(x_p,*CTRFitParam), color='k')
    ax1.set_ylim(0,1.05*max(values))
    T_y = generateTextY(9,1.05*max(values))
    ax1.text(T_x,T_y[0],"Data Stats:")
    ax1.text(T_x,T_y[1],"mean: "+format(time_diff.mean(),"5.2f"))
    ax1.text(T_x,T_y[2],"std: "+format(time_diff.std(),"7.4f"))
    ax1.text(T_x,T_y[3],"count: "+str(len(time_diff)))
    ax1.text(T_x,T_y[4],"Fit Paramters:")
    ax1.text(T_x,T_y[5],"mu: "+format(CTRFitParam[1],"5.2f"))
    ax1.text(T_x,T_y[6],"sigma: "+format(CTRFitParam[2],"7.4f"))
    ax1.text(T_x,T_y[7],"A: "+format(CTRFitParam[0],"5.2f"))
    ax1.text(T_x,T_y[8],"FWHM: "+format(2.355*CTRFitParam[2],"7.4f"))
    # Save the figure to a .png
    plt.savefig("CTR - Pair {}_{}.png".format(int(toGeoChannelID(LChannel)), int(toGeoChannelID(RChannel))))
    if display:
        plt.show()
        plt.clf()
        plt.close()
    else:
        plt.clf()
        plt.close()
def toGeo(x):
    #converts PETSys ID to geometric ID
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

geo_channels = []
for i in range(128):
    geo_channels.append([i,toGeo(i)])
geo_channels = np.asarray(geo_channels)

def toGeoChannelID(AbsChannelID):
    # Convert PETSys absolute channel IDs to geomteric IDs
    portID = AbsChannelID // 131072
    slaveID = (AbsChannelID - 131072*portID) // 4096
    chipID = (AbsChannelID - slaveID*4096 - 131072*portID) // 64
    channelID = AbsChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    AbsPCB_ChanID = geo_channels[geo_channels[:,0] == PCB_ChanID][0][1]

    #General formula can be found in above function "to AbsChannelID"
    GeoChannelID = 10**6 * portID + 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64
    return GeoChannelID
