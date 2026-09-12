import numpy as np
import pandas as pd
from tqdm import tqdm

def filter_data(data_file, stats_dir, output_file, sigma=3):
    # Load the data file into a DataFrame
    data_df = pd.read_csv(data_file, sep='\t', names=["TimeLeft", "EnergyLeft", "ChannelIDLeft", "TimeRight", "EnergyRight", "ChannelIDRight"], usecols=(2,3,4,7,8,9), dtype=np.float32)
    
    # Load left and right statistics files containing mean energy and standard deviation per channel
    left_stats = pd.read_csv(f"{stats_dir}/averaged_data_left.tsv", sep='\t')
    right_stats = pd.read_csv(f"{stats_dir}/averaged_data_right.tsv", sep='\t')
    
    # Compute global mean values and standard deviations for missing channel IDs
    global_Lmean, global_Lstd = left_stats["LenergyMean"].mean(), left_stats["LenergyStd"].mean()
    global_Rmean, global_Rstd = right_stats["RenergyMean"].mean(), right_stats["RenergyStd"].mean()
    
    # Create lookup dictionaries for quick access
    left_stats_dict = left_stats.set_index("LChannel")[['LenergyMean', 'LenergyStd']].to_dict(orient='index')
    right_stats_dict = right_stats.set_index("RChannel")[['RenergyMean', 'RenergyStd']].to_dict(orient='index')
    
    # Define filtering function
    def filter_row(row):
        # Retrieve designated mean/std from dictionaries or fallback to global values
        left_data = left_stats_dict.get(row["ChannelIDLeft"], {"LenergyMean": global_Lmean, "LenergyStd": global_Lstd})
        right_data = right_stats_dict.get(row["ChannelIDRight"], {"RenergyMean": global_Rmean, "RenergyStd": global_Rstd})
        
        left_min = left_data["LenergyMean"] - sigma * left_data["LenergyStd"]
        right_min = right_data["RenergyMean"] - sigma * right_data["RenergyStd"]
        
        return (row["EnergyLeft"] > left_min) and (row["EnergyRight"] > right_min)
    
    # Apply filtering while maintaining order
    filtered_df = data_df[data_df.apply(filter_row, axis=1)]
    
    # Save filtered data
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    return filtered_df

# Example usage
data_file = "/home/michaelgajda/mda/postSpillDirTPPT/trip2/postSpillData_Run5_TPPT_Acrylic1_NoCollimator_1shot_15min_HWTrigOn_coinc.dat"
stats_dir = "/home/michaelgajda/mda/postSpillDirTPPT/channel_pair_energy_map"
output_file = "/home/michaelgajda/mda/postSpillDirTPPT/two_point_five_sigma_data/Filtered_Run12.dat"

# Apply the filtering process
filtered_df = filter_data(data_file, stats_dir, output_file)
print(filtered_df.head())
