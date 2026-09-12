import numpy as np
import os
from tqdm import tqdm

name = 'htof'
#maskname = 'mcirc_geommask'
maskname = f'{name}_pcut_mask'
data_dir = f'/home/kale-chen/Documents/PET/TPPT2026/data/'
data_path = os.path.join(data_dir, f'{name}.dat')
mask_path = os.path.join(data_dir, f'{maskname}.bin')

num_rows = os.path.getsize(mask_path) // (1)
mask = np.memmap(mask_path, dtype = np.bool, mode = 'r', shape = (num_rows))

i = 0
kept = 0
with open(f'{data_path.replace('.dat', 'p.dat')}', 'w') as f:
    for line in tqdm(open(data_path, 'r')):
        if mask[i]:
            f.write(line)
            kept += 1
        i += 1
print(f'Kept {kept} events')