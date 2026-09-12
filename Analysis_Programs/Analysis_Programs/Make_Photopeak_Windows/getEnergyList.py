import pandas as pd
import glob
import os

def read_and_process_csvs(directory, output_file, stats_output_file):
    # Define column names
    columns = ["LChannel", "RChannel", "counts", "Left Mean", "Right Mean", 
               "left sigma", "right sigma", "flag1", "flag2", "flag3"]
    
    # Read all CSV files
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    dataframes = []
    
    for file in csv_files:
        df = pd.read_csv(file, sep='\t', names=columns, engine='python')
        dataframes.append(df)
    
    # Concatenate all dataframes
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Group by LChannel and calculate:
    # - mean of Left Mean
    # - mean of left sigma
    left_stats = combined_df.groupby("LChannel").agg({
        "Left Mean": "mean",
        "left sigma": "mean"
    }).reset_index()
    left_stats.columns = ["LChannel", "LenergyMean", "LenergyStd"]
    
    # Group by RChannel and calculate:
    # - mean of Right Mean
    # - mean of right sigma
    right_stats = combined_df.groupby("RChannel").agg({
        "Right Mean": "mean",
        "right sigma": "mean"
    }).reset_index()
    right_stats.columns = ["RChannel", "RenergyMean", "RenergyStd"]
    
    # Compute the mean standard deviation excluding NaNs
    mean_LenergyStd = left_stats["LenergyStd"].mean(skipna=True)
    mean_RenergyStd = right_stats["RenergyStd"].mean(skipna=True)
    
    # Replace NaNs with the mean standard deviation
    left_stats["LenergyStd"].fillna(mean_LenergyStd, inplace=True)
    right_stats["RenergyStd"].fillna(mean_RenergyStd, inplace=True)
    
    # Save the results
    left_stats.to_csv(output_file.replace(".tsv", "_left.tsv"), sep='\t', index=False)
    right_stats.to_csv(output_file.replace(".tsv", "_right.tsv"), sep='\t', index=False)
    
    return left_stats, right_stats

# Example usage
directory = "/home/michaelgajda/mda/postSpillDirTPPT/channel_pair_energy_map"
averaged_output_file = "/home/michaelgajda/mda/postSpillDirTPPT/channel_pair_energy_map/averaged_data.tsv"
stats_output_file = "/home/michaelgajda/mda/postSpillDirTPPT/channel_pair_energy_map/statistics.csv"

left_stats, right_stats = read_and_process_csvs(directory, averaged_output_file, stats_output_file)

# To check results
print(left_stats.head())
print(len(left_stats))
print(right_stats.head())
print(len(right_stats) + len(left_stats))
