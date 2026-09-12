import pandas as pd
from tqdm import tqdm

dir = "/home/kilroy101/path/Geant4/Expirimental_Data/"
chunks = pd.read_csv(dir+"Center_Source_8-8_coinc.dat", sep="\t", header=None, low_memory=False, usecols=[2,4,7,9],chunksize=1000000)

#fuck you Nuno lazy prick
def toRecoChannelID(AbsChannelID): #updated to account for port ID
    portID = AbsChannelID // 131072 # new
    slaveID = (AbsChannelID - 131072*portID) // 4096
    chipID = (AbsChannelID - 4096*slaveID - 131072*portID) // 64
    channelID = AbsChannelID % 64
    RecoChannelID = 3072*portID + 1024*slaveID + 64*chipID + channelID
    return RecoChannelID

df = pd.DataFrame()
for chunk in tqdm(chunks):
    chunk.columns = ["TimeL", "ChannelIDL","TimeR","ChannelIDR"]
    chunk.ChannelIDL = chunk.ChannelIDL.apply(toRecoChannelID)
    chunk.ChannelIDR = chunk.ChannelIDR.apply(toRecoChannelID)
    chunk["d1"] = chunk["ChannelIDL"]
    chunk["d2"] = chunk["ChannelIDR"]
    chunk["deltat"] = chunk["TimeL"] - chunk["TimeR"]
    chunk["t1"] = 10**-6 * chunk["TimeL"]
    chunk = chunk.drop(columns=["ChannelIDL", "ChannelIDR", "TimeR", "TimeL"])
    df = pd.concat([df,chunk],ignore_index=True)


df.to_csv("~/path/Geant4/Expirimental_Data/Scanner_Commissioning/Center_Source_8-8_reco.lm", header=False, index=False)
