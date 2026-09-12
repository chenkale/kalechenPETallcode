import numpy as np
import os
from tqdm import tqdm

name = 'pointsrpg'
#maskname = 'mcirc_geommask'
maskname = 'pointsrp_intertimesmask'
data_dir = f'/home/kale-chen/Documents/PET/TPPT2026/data/'
data_path = os.path.join(data_dir, f'{name}.bin')
num_cols = 4
mask_path = os.path.join(data_dir, f'{maskname}.bin')

num_rows = os.path.getsize(data_path) // (num_cols * 2)
num_rows_mask = os.path.getsize(mask_path) // (1)
print(num_rows, num_rows_mask)
if num_rows != num_rows_mask:
    raise ValueError('Number of rows in data and mask do not match')

data = np.memmap(data_path, dtype = np.int16, mode = 'r', shape = (num_rows, num_cols))
mask = np.memmap(mask_path, dtype = np.bool, mode = 'r', shape = (num_rows))

chunksize = 1000000
event_count = 0
if os.path.exists(f'{data_path.replace('.bin', '_masked.bin')}'):
    os.remove(f'{data_path.replace('.bin', '_masked.bin')}')
for i in tqdm(range(0, num_rows, chunksize)):
    maskchunk = mask[i:i+chunksize]
    chunk = data[i:i+chunksize]
    chunk = chunk[maskchunk]
    with open(f'{data_path.replace('.bin', '_masked.bin')}', 'ab') as f:
        chunk.tofile(f)
    event_count += len(chunk)
print(event_count)  