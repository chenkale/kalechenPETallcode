import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from tqdm import tqdm
from scipy.optimize import curve_fit

name = 'pointsrpgd'
data_dir = f'/home/kale-chen/Documents/PET/TPPT2026/data/{name}.bin'
num_cols, data_size = 4, 2
# Define time ranges to splice on
timeranges = [(0, 5), (70, 80), (145, 155), (220, 230), (295, 305), (370, 380), (445, 455), (520, 530)]
timeranges = [(5, 55), (65, 115), (125, 175), (185, 235), (245, 295), (305, 355), (365, 415), (425, 475)]
num_rows = os.path.getsize(data_dir) // (num_cols * data_size)
data = np.memmap(data_dir, dtype=np.int16, mode='r', shape=(num_rows, num_cols))
print(num_rows,)
print(data[0:5])

for i in range(len(timeranges)):
    timerange = timeranges[i]
    totalevents = 0
    if os.path.exists(f'/home/kale-chen/Documents/PET/TPPT2026/data/{name}{i}.bin'):
        os.remove(f'/home/kale-chen/Documents/PET/TPPT2026/data/{name}{i}.bin')
    with open(f'/home/kale-chen/Documents/PET/TPPT2026/data/{name}{i}.bin', 'ab') as f:
        for j in range(0, num_rows, 1000000):
            chunk = data[j:j+1000000]
            if num_cols == 4:
                mask = (chunk[:, 3] >= timerange[0]) & (chunk[:, 3] < timerange[1])
                chunk = chunk[mask]
            else:
                mask = (chunk[:, 5] >= timerange[0]) & (chunk[:, 5] < timerange[1])
                chunk = chunk[mask]
            totalevents += len(chunk)
            chunk.tofile(f)
    print(totalevents)