import pandas as pd

df = pd.read_csv("map_trigger.tsv", sep="\t", header=None)

zone1 = [4,5,6,7,20,21,22,23, 28,29,30,31]
zone2 = [36,37,38,39,44,45,46,47,60,61,62,63]

for region1 in zone1:
    for region2 in zone2:
        df = pd.concat([df, pd.DataFrame({0: [region1], 1:[region2], 2:["C"]})], ignore_index=True)

df.to_csv("map_trigger.tsv", sep='\t', header=False, index=False)