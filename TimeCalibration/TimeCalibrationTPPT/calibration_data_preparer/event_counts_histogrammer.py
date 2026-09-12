# Events counts histogrammer 

import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

data_dir = r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\compressed_events_cut.bin'

lor_pops = np.zeros(3072 * 3072, dtype = np.int32)

num_rows = os.path.getsize(data_dir) // 10 # 10 = 2 (size of int16) * 5 (entries per row)
data = np.memmap(data_dir, dtype=np.int16, mode='r', shape=(num_rows,5))
for i in tqdm(range(0, num_rows, 2000000)):
    chunk = data[i:i+2000000] # Selects ID information only
    ids = np.asarray(chunk[:, [1, 3]], dtype = np.int32)
    indices = ids[:, 0] * 3072 + ids[:, 1]
    np.add.at(lor_pops, indices, 1)
        
#lor_pops.tofile(r'C:\Users\burri\Documents\PET\TimeCalibration\CalibrationDataPreparer\event_counts.bin')

# Histogram of LOR populations
print(f'Total events: {sum(lor_pops)}')
plt.hist(lor_pops, bins = np.linspace(0, 1000, 1000), color = 'blue', alpha = 0.4)
plt.yscale('log')
plt.title('LOR Populations')
plt.xlabel('Pop')
plt.ylabel('Counts')
plt.grid()
plt.show()