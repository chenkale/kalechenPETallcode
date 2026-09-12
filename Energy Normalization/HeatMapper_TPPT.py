from ChannelPairHeader_TPPT import *
import pandas as pd
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import axes3d, Axes3D
import math
import matplotlib
import statistics as stat
pd.set_option("display.max_rows",1000)

# Code to create heat maps showing the number of photo peak coincidence events seen in each pixel of each array before and after
# normalization is applied to each LOR (channel pair). Takes the output .csv file produced from ChannelPairBuilder.py as the input

# --- Set Debug flag to true to turn on various "print" commands ---
Debug = True

# ===================================================================================================================
# ======================================== SECTION 1: Read and organize data ========================================
# ===================================================================================================================

# --- Choose while file to read in and strip the core info out for saving images to output filenames ---

# Your relative working directory
datadir = "NewOutputs/TPPT/"

# New Outputs
## MiniPET
file_name = "Output_Combined_NormalizationRun1_circle_calibration_30min_ccw_Run1_TrigOn_coinc"

# Old Outputs
#file_name = "Output_NormalizationRun2_210mmSep_21mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun3_210mmSep_23mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun4_210mmSep_25mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun5_210mmSep_27mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun6_210mmSep_29mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun7_210mmSep_31mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun8_210mmSep_33mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun9_210mmSep_35mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun10_210mmSep_37mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun11_210mmSep_39mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun12_210mmSep_41mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun13_210mmSep_43mmSourceHeight_coinc"
#file_name = "Output_NormalizationRun14_210mmSep_45mmSourceHeight_coinc"
vmax = 6000 # Set the vmax value for the color scale to make sense depending on the input data
           # Use this value for any of the above files with lower pp counts per pixel

#file_name = "Output_SourcePlaneNormalization_210mmSep"
#vmax = 5000

base_name = file_name.lstrip("Output_")
base_name = base_name.rstrip("_coinc")
print(base_name)

# --- Read the data into a pandas dataframe, keeping only the columns of interest (first three) related to the ChannelIDs and number of PP events ---
df = pd.read_csv(datadir+"{}.csv".format(file_name),
                     sep="\t",
                     low_memory=False,
                     header=None,
                     names=["ChannelIDL", "ChannelIDR", "PP_num_evts"],
                     usecols=[0, 1, 2])
#print(df)

# --- Append the corresponding GeoChannelIDs to the dataframe --- 
Geo_L = []
Geo_R = []
# TODO: Check if this is outputting the correct numbers
for i in range(df.ChannelIDL.count()):
    Geo_L.append(toGeoChannelID(df["ChannelIDL"][i])) # Port 8 was -864 for MiniPET (check to see if math is needed here, maybe -1000000?
    Geo_R.append(toGeoChannelID(df["ChannelIDR"][i])) # Port 1 in FEB/D, "Right" module when facing in beam direction (downstream) -> PCB #4
df["GeoChannelIDL"] = Geo_L
df["GeoChannelIDR"] = Geo_R
df["ChanPair"] = df["GeoChannelIDL"].astype(str) + "-" + df["GeoChannelIDR"].astype(str)
if(Debug):
    print("1")
    print(df.sort_values(by="GeoChannelIDL"))

# --- Reading in normalization constants .csv file into a second dataframe ---
df2 = pd.read_csv(datadir+"LOR_Geo_Norm_Consts.csv", sep="\t", low_memory=False, header=None, names=["GeoChannelIDL", "GeoChannelIDR", "Norm_Const"])
df2["ChanPair"] = df2["GeoChannelIDL"].astype(str) + "-" + df2["GeoChannelIDR"].astype(str)
if(Debug):
    print("2")
    print(df2)

# --- Merge the dataframes (omitting rows from df2 with no matching channel pair) ---
merged_df = df.merge(df2, how='right', on ='ChanPair')
merged_df["Scaled_PP_num_evts"] = merged_df["PP_num_evts"]/merged_df["Norm_Const"]
if(Debug):
    print("3")
    print(merged_df.sort_values(by="PP_num_evts"))

# --- Output merged df to sorted debug file for more easily checking results ---
#merged_df = merged_df.astype({'ChanIDL':'int32','ChanIDR':'int32'})
#np.savetxt("HeatMapper_Debug.txt", merged_df, fmt = (('%03i     %03i  %-i  %03i  %03i  %s  %03i  %03i  %1.6f  %1.6f')), delimiter = "\t")

# ==================================================================================================================
# ==================== SECTION 2: Generate PRE-corrected heat maps of PP occupancy per channel ====================
# ==================================================================================================================

# --- Setup heat map visualization parameters ---
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('jet') # cmap options: binary (default), *spring*, *summer*, autumn, winter, *cool*, Wistia, *hot*, copper, *hsv*, *rainbow*, *jet*
cmap = truncate_colormap(cmap, 0.2, 0.8)

# --- Generating left array heat map elements by looping over ChannelIDL and summing PP events across all LORs ---
Evts_L = pd.DataFrame(columns=["ChannelID_L", "PP_Evts_L"])
for i in range(0,128):
    run = df[df["GeoChannelIDL"] == i]
    #print(run)
    if run.ChannelIDL.count() == 0:
        Evts_L = pd.concat([Evts_L, pd.DataFrame([{"ChannelID_L": i, "PP_Evts_L": 0.}])], ignore_index=True)
    else:
        SumEvts_L = run.PP_num_evts.sum()
        Evts_L = pd.concat([Evts_L, pd.DataFrame([{"ChannelID_L": i, "PP_Evts_L": SumEvts_L}])], ignore_index=True)
if(Debug):
    print("4")
    print(Evts_L)

# --- Separating elements into "Top" and "Bottom" arrays based on geo ChannelIDs ---
Top_L = Evts_L[Evts_L["ChannelID_L"] < 64]
Bot_L = Evts_L[Evts_L["ChannelID_L"] >= 64]

# --- Setting any empty/zero pixels to "NaN" so that they appear blank on the heat maps (rather than the color associated with min_value) ---
Top_L[Top_L == 0] = np.nan
Bot_L[Bot_L == 0] = np.nan

# --- Reshape data into an 8x8 array for creating the heat map ---
Top_L = np.reshape(np.array(Top_L.PP_Evts_L), (8,8))
Bot_L = np.reshape(np.array(Bot_L.PP_Evts_L), (8,8))
if(Debug):
    print("5")
    print(Top_L)
    print(Bot_L)

# --- Now going through the same steps for the right array data ---
Evts_R = pd.DataFrame(columns=["ChannelID_R", "PP_Evts_R"])
for i in range(0,128):
    run = df[df["GeoChannelIDR"] == i]
    if run.ChannelIDR.count() == 0:
        Evts_R = pd.concat([Evts_R, pd.DataFrame([{"ChannelID_R": i, "PP_Evts_R": 0.}])], ignore_index=True)
    else:
        SumEvts_R = run.PP_num_evts.sum()
        Evts_R = pd.concat([Evts_R, pd.DataFrame([{"ChannelID_R": i, "PP_Evts_R": SumEvts_R}])], ignore_index=True)
if(Debug):
    print("6")
    print(Evts_R)
      
Top_R = Evts_R[Evts_R["ChannelID_R"] < 64]
Bot_R = Evts_R[Evts_R["ChannelID_R"] >= 64]

Top_R[Top_R == 0] = np.nan
Bot_R[Bot_R == 0] = np.nan

Top_R = np.reshape(np.array(Top_R.PP_Evts_R), (8,8))
Bot_R = np.reshape(np.array(Bot_R.PP_Evts_R), (8,8))
if(Debug):
    print("7")
    print(Top_R)
    print(Bot_R)

# --- Putting all the heat maps (left array top and bottom + right array top and bottom) on one 2x2 figure ---
plt.figure(figsize=(10,8))
plt.subplot(221) # Going into subfigure 1 of 2x2 figure array
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))

# --- Loop through each 8x8 array to add text into each pixel corresponding to the data value of that pixel (or set to zero if NaN) ---
for i in range(8):
    for j in range(8):
        if (Top_L[j][i] > 0):
            plt.text(i, j, format(int(Top_L[j][i]),"3"), ha="center", va="center", color="b")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="b")
plt.title("'Left' Array - Top")
plt.imshow(Top_L,cmap=cmap,vmin=0, vmax=vmax)

# --- Repeat the above but for subfigure 3 of 2x2 array (corresponds to subfigure below subfigure 1) ---
plt.subplot(223)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        if (Bot_L[j][i] > 0):
            plt.text(i, j, format(int(Bot_L[j][i]),"3"), ha="center", va="center", color="b")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="b")
plt.title("'Left' Array - Bottom")
plt.imshow(Bot_L,cmap=cmap,vmin=0, vmax=vmax)

# --- Repeat the above but for subfigure 2 of 2x2 array (corresponds to subfigure to the right of subfigure 1) ---
plt.subplot(222)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        if (Top_R[j][i] > 0):
            plt.text(i, j, format(int(Top_R[j][i]),"3"), ha="center", va="center", color="r")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="r")
plt.title("'Right' Array - Top")
plt.imshow(Top_R, cmap=cmap, interpolation = 'none', vmin=0, vmax=vmax)

# --- Repeat the above but for subfigure 4 of 2x2 array (corresponds to subfigure below subfigure 2) ---
plt.subplot(224)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        if (Bot_R[j][i] > 0):
            plt.text(i, j, format(int(Bot_R[j][i]),"3"), ha="center", va="center", color="r")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="r")
plt.title("'Right' Array - Bottom")
plt.imshow(Bot_R, cmap=cmap, interpolation = 'none', vmin=0, vmax=vmax)

# --- Adjust various plotting parameters (margins, axis, etc.) and save and display the figure ---
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
cax = plt.axes([0.85, 0.1, 0.025, 0.8])
plt.colorbar(cax=cax)
plt.savefig("Heat_Map_{}.eps".format(base_name), format='eps')
plt.show()
plt.clf()
plt.close()

# ==================================================================================================================
# ==================== SECTION 3: Generate POST-corrected heat maps of PP occupancy per channel ====================
# ==================================================================================================================

# --- Need to re-scale contributions from each LOR by the norm factor merged into the dataframe earlier ---
# --- Generating left array heat map elements by looping over ChannelIDL and summing SCALED PP events across all LORs ---
# --- Procedure is the exact same but just sums of the column named "Scaled_PP_num_evts" in the merged dataframe ---
NormEvts_L = pd.DataFrame(columns=["ChannelID_L", "Scaled_PP_Evts_L"])
for i in range(0,128):
    run = merged_df[merged_df["GeoChannelIDL_x"] == i]
    if(Debug):
        print("8")
        print(run)
    if run.ChannelIDL.count() == 0:
        NormEvts_L = pd.concat([
            NormEvts_L,
            pd.DataFrame([{"ChannelID_L": i, "Scaled_PP_Evts_L": 0.}])
        ], ignore_index=True)
    else:
        NormSumEvts_L = run.Scaled_PP_num_evts.sum()
        NormEvts_L = pd.concat([
            NormEvts_L,
            pd.DataFrame([{"ChannelID_L": i, "Scaled_PP_Evts_L": NormSumEvts_L}])
        ], ignore_index=True)
if(Debug):
    print("9")
    print(NormEvts_L)
        
NormTop_L = NormEvts_L[NormEvts_L["ChannelID_L"] < 64]
NormBot_L = NormEvts_L[NormEvts_L["ChannelID_L"] >= 64]

NormTop_L[NormTop_L == 0] = np.nan
NormBot_L[NormBot_L == 0] = np.nan

NormTop_L = np.reshape(np.array(NormTop_L.Scaled_PP_Evts_L), (8,8))
NormBot_L = np.reshape(np.array(NormBot_L.Scaled_PP_Evts_L), (8,8))
if(Debug):
    print("10")
    print(NormTop_L)
    print(NormBot_L)

# --- Following the same procedure for the right array ---
NormEvts_R = pd.DataFrame(columns=["ChannelID_R", "Scaled_PP_Evts_R"])
for i in range(0,128):
    run = merged_df[merged_df["GeoChannelIDR_x"] == i]
    if run.ChannelIDR.count() == 0:
        NormEvts_R = pd.concat([
            NormEvts_R,
            pd.DataFrame([{"ChannelID_R": i, "Scaled_PP_Evts_R": 0.}])
        ], ignore_index=True)
    else:
        NormSumEvts_R = run.Scaled_PP_num_evts.sum()
        NormEvts_R = pd.concat([
            NormEvts_R,
            pd.DataFrame([{"ChannelID_R": i, "Scaled_PP_Evts_R": NormSumEvts_R}])
        ], ignore_index=True)
if(Debug):
    print("11")
    print(NormEvts_R)

        
NormTop_R = NormEvts_R[NormEvts_R["ChannelID_R"] < 64]
NormBot_R = NormEvts_R[NormEvts_R["ChannelID_R"] >= 64]

NormTop_R[NormTop_R == 0] = np.nan
NormBot_R[NormBot_R == 0] = np.nan

NormTop_R = np.reshape(np.array(NormTop_R.Scaled_PP_Evts_R), (8,8))
NormBot_R = np.reshape(np.array(NormBot_R.Scaled_PP_Evts_R), (8,8))
if(Debug):
    print("12")
    print(NormTop_R)
    print(NormBot_R)

# --- Putting both heat maps on one 4x4 figure ---
plt.figure(figsize=(10,8))
plt.subplot(221)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        #if (NormTop_L[j][i] > 0):
        if not (math.isnan(NormTop_L[j][i])) and not (math.isinf(NormTop_L[j][i])):
            plt.text(i, j, format(int(NormTop_L[j][i]),"3"), ha="center", va="center", color="b")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="b")
plt.title("'Left' Array - Top")
plt.imshow(NormTop_L,cmap=cmap,vmin=0, vmax=vmax)

plt.subplot(223)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        if not (math.isnan(NormBot_L[j][i])) and not (math.isinf(NormBot_L[j][i])):
            plt.text(i, j, format(int(NormBot_L[j][i]),"3"), ha="center", va="center", color="b")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="b")
plt.title("'Left' Array - Bottom")
plt.imshow(NormBot_L,cmap=cmap,vmin=0, vmax=vmax)

# plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
# cax = plt.axes([0.85, 0.1, 0.025, 0.8])
# plt.colorbar(cax=cax)
# plt.show()

#plt.figure(figsize=(5,8))
plt.subplot(222)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        if not (math.isnan(NormTop_R[j][i])) and not (math.isinf(NormTop_R[j][i])):
            plt.text(i, j, format(int(NormTop_R[j][i]),"3"), ha="center", va="center", color="r")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="r")
plt.title("'Right' Array - Top")
plt.imshow(NormTop_R, cmap=cmap, interpolation = 'none', vmin=0, vmax=vmax)

plt.subplot(224)
plt.xticks(np.arange(8))
plt.yticks(np.arange(8))
for i in range(8):
    for j in range(8):
        if not (math.isnan(NormBot_R[j][i])) and not (math.isinf(NormBot_R[j][i])):        
            plt.text(i, j, format(int(NormBot_R[j][i]),"3"), ha="center", va="center", color="r")
        else:
            plt.text(i, j, 0, ha="center", va="center", color="r")
plt.title("'Right' Array - Bottom")
plt.imshow(NormBot_R, cmap=cmap, interpolation = 'none', vmin=0, vmax=vmax)

plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,left=0.2)
cax = plt.axes([0.85, 0.1, 0.025, 0.8])
plt.colorbar(cax=cax)
plt.savefig("Scaled_Heat_Map_{}.eps".format(base_name), format='eps')
plt.show()
plt.clf()
plt.close()
