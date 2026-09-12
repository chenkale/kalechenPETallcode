import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

# only part you need to change
coincidence_file = "Positronioum_Na-22button_30min_HWtrigOFF_SWtrigOFF_coinc.dat"
geometry_file = "TPPT-Geo.csv"


# we need Reconstruction (reco) ID to reference geometry tables
def toRecoChannelID(AbsChannelID):
    portID = AbsChannelID // 131072 # new
    slaveID = (AbsChannelID - 131072*portID) // 4096
    chipID = (AbsChannelID - 4096*slaveID - 131072*portID) // 64
    channelID = AbsChannelID % 64
    RecoChannelID = 3072*portID + 1024*slaveID + 64*chipID + channelID
    return RecoChannelID
    
def plot_active_pixels(coincidence_df,geometry_df): # functionalized so we can do this for MiniPET or TPPT
    
    litRight = np.unique(coincidence_df.ChannelIDR)
    litLeft = np.unique(coincidence_df.ChannelIDL)
    
    fig = plt.figure(figsize = (20,20))
    ax = plt.axes(projection='3d')

    print("Finding Active Channels:")
    inactive_chans = []
    for index,i,j,k in tqdm(zip(geometry_df.index,geometry_df.x,geometry_df.y,geometry_df.z)):
        if index in litRight or index in litLeft:
            ax.scatter(i,j,k,color = "blue")
        else:
            ax.scatter(i,j,k,color = "red")
            inactive_chans.append(index)
            
    
            
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.zaxis.set_tick_params(labelsize=20)
    ax.set_xlabel("X-axis [mm]",fontsize = 20,labelpad=20)
    ax.set_ylabel("Y-axis [mm]",fontsize = 20,labelpad=20)
    ax.set_zlabel("Z-axis [mm]",fontsize = 20,labelpad=20)
    fig.suptitle("3D Projection of Active Channels",fontsize = 30,y=0.92)
    plt.tight_layout()
    
    legend_handles = [Line2D([], [], marker='.', color='blue', linestyle='None',markersize = 17),
          Line2D([], [], marker='.', color='red', linestyle='None',markersize = 17)]
    fig.legend(legend_handles,["Active Channels","Inactive Channels"],fontsize = 20,bbox_to_anchor = (1,0.8))
    
    print("\n########################################################")
    print("There are {} inactive channels out of {} total Channels".format(len(inactive_chans),len(geometry_df)))
    print("########################################################\n")
    
    plt.savefig("Active_Channel_Map.pdf") # only option to save because its a big plot and needs a viewer
    
    return inactive_chans
    
# I only sample the first 100000 data points, if for whatever reason you want more change nrows in read_csv() below or remove the parameter altogether to read in all of the data
coincidence_df = pd.read_csv(coincidence_file,sep='\t',header = None,usecols=[2,3,4,7,8,9],nrows = 100000)
coincidence_df.columns = ['TimeL', 'ChargeL', 'ChannelIDL', 'TimeR', 'ChargeR', 'ChannelIDR']
coincidence_df["ChannelIDL"] = coincidence_df["ChannelIDL"].apply(toRecoChannelID)
coincidence_df["ChannelIDR"] = coincidence_df["ChannelIDR"].apply(toRecoChannelID)

# read in geometry LUT
geometry_df = pd.read_csv(geometry_file,usecols = [0,1,2],header = None)
geometry_df.columns = ['x','y','z']

# call function to plot
inactive_channels = plot_active_pixels(coincidence_df,geometry_df)

print("The inactive channel IDs are:",inactive_channels)
