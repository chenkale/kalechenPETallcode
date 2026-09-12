# Event cutter

import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

data_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\compressed_events_cut.bin'
out_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\compressed_events_cut.bin'

cut_large_time_diffs = True # Cut events with large time differences, unphysical
counts_threshold = 100 # Only keep LORs with this number of events or more

lor_pops = np.zeros(3072 * 3072, dtype = np.int32)

num_rows = os.path.getsize(data_dir) // 10 # 10 = 2 (size of int16) * 5 (entries per row)
data = np.memmap(data_dir, dtype=np.int16, mode='r', shape=(num_rows,5))
for i in tqdm(range(0, num_rows, 1000000)):
    chunk = data[i:i+1000000] # Selects ID information only
    if cut_large_time_diffs:
        chunk = chunk[(chunk[:, 4] > -2500) & (chunk[:, 4] < 2500)]
    ids = np.asarray(chunk[:, [1, 3]], dtype = np.int32)
    indices = ids[:, 0] * 3072 + ids[:, 1]
    np.add.at(lor_pops, indices, 1)

if counts_threshold > 0:
    mask = np.asarray([pop >= counts_threshold for pop in lor_pops])
    for i in tqdm(range(0, num_rows, 1000000)):
        chunk = data[i:i+1000000]
        if cut_large_time_diffs:
            chunk = chunk[(chunk[:, 4] > -2500) & (chunk[:, 4] < 2500)].astype(np.int32)
        ids = np.asarray(chunk[:,1] * 3072 + chunk[:,3], dtype = np.int32)
        # Apply mask: keep rows where the ID is included
        keep = mask[ids]  # mask[ids] gives True/False for each row
        chunk = chunk[keep].astype(np.int16)
        with open(out_dir, 'ab') as f:
            chunk.tofile(f)

for x in tqdm(range(len(lor_pops))):
    if lor_pops[x] < counts_threshold:
        lor_pops[x] = 0
lor_pops.tofile(r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\event_counts_cut1.bin')
