# Master script for CASToR

# Castor Streamlined

'''
Paths need to be specified
'''
name = 'run3_pcut_tcut'

data_file = f'/home/kale-chen/Documents/PET/MDA_04112026/Data/{name}.bin'
#data_file = f'/home/kale-chen/Documents/PET/TimeCalibration/MiniPET/{name}.bin'
num_cols = 3
#data_file = '/home/kale-chen/Documents/PET/test_map.bin'
#name = 'test_map'
geometry_file = '/home/kale-chen/Documents/PET/TPPT_Scanner_map_adjusted.csv'
#geometry_file = '/home/kale-chen/Documents/PET/test_map.csv'
#geometry_file = '/home/kale-chen/Documents/PET/remap.csv'
#geometry_file = '/home/kale-chen/Downloads/TPPT_Scanner_map.csv'
config_dir = '/home/kale-chen/Documents/CASToR/castor_v3.2/config/scanner'
castor_dir = '/home/kale-chen/Documents/CASToR/castor_v3.2/build'
output_dir = '/home/kale-chen/Documents/CASToR/images/'
sens_file = '/home/kale-chen/Documents/CASToR/images/hcenter_pcut_a/hcenter_pcut_a_sensitivity.hdr'
atn_file = '/home/kale-chen/Documents/CASToR/images/mumap/mumap512.hdr'

normalization = False
normalization_file = '/home/kale-chen/Documents/PET/Energy Normalization/event_counts_norm.bin'

dim = (256, 256, 256)
vox = (0.5, 0.5, 0.5)

total_duration = 900
frame_duration = 30



'''
Code begins here
'''


from configparser import NoOptionError
import os
import subprocess
from tqdm import tqdm
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from tqdm import tqdm

start = time.time()

def geometry_to_lut():
    scanner = np.asarray(pd.read_csv(geometry_file, dtype='float32', delimiter=',', header=None))
    scanner.tofile(os.path.join(config_dir, 'TPPT.lut'))

def write_hscan():
    with open(os.path.join(config_dir, 'TPPT.hscan'), 'w') as f:
        f.write('scanner name: TPPT\n' +
        'modality:    PET\n' +
        'description: This scanner description is based on an actual scanner, however, this implementation is not supported nor validated by its manufacturer.\n' +
        'scanner radius:    ' + '167.5' + '\n' +
        'number of rings in scanner:    32\n' +
        'number of elements:    6144\n' +
        #'number of elements:    512\n' +
        #'number of elements:    1041\n' +
        'number of layers:    1\n' +
        'number of crystals in layer(s):   6144\n' +
        #'number of crystals in layer(s):   512\n' +
        #'number of crystals in layer(s):   1041\n' +
        'crystals size depth:    15\n' +
        'crystals size transaxial:    3.00\n' +
        'crystals size axial:    3.00\n' +
        'voxels number transaxial:    256\n' +
        'voxels number axial:    128\n' + 
        'field of view transaxial:    4000\n' +
        'field of view axial:    4000\n' +
        'min angle difference:    0' )

def make_cdf(name, num_columns):
    if os.path.exists(os.path.join(output_dir, name) + '.cdf'):
        os.remove(os.path.join(output_dir, name) + '.cdf')
    num_rows = os.path.getsize(data_file) // (num_columns * 2)
    data = np.memmap(data_file, dtype = np.int16, mode = 'r', shape = (num_rows, num_columns))
    
    # Handle normalization information
    if normalization:
        normalization_data = np.fromfile(normalization_file, dtype = np.float32).reshape(3072, 3072)
    for i in tqdm(range(0, num_rows, 1000000)):
        chunk = np.copy(data[i:i+1000000])
        #chunk[:, 0] = chunk[:, 0] + 256
        #chunk[:, 2] = chunk[:, 2] // 2 # Time diff to TOF
        if normalization:
            towrite = np.zeros((len(chunk), 5), dtype = np.float32) # For reordering
            norm = np.array(normalization_data[chunk[:, 0], chunk[:, 1]], dtype = np.float32)
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            towrite[:, 2] = norm
            towrite[:,[1,3,4]] = chunk[:,[2,0,1]]
            print(towrite[0:5])
        else:
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            # Structured array towrite has [time (uint32), tof (float32), idl (uint32), idr (uint32)]
            towrite = np.zeros((len(chunk)), dtype = [('time', 'uint32'), ('tof', 'float32'), ('idl', 'uint32'), ('idr', 'uint32')])
            towrite['tof'] = chunk[:, 2].astype(np.float32)
            towrite['idl'] = chunk[:, 0].astype(np.uint32)
            towrite['idr'] = chunk[:, 1].astype(np.uint32)
            #towrite[:,[1,2]] = chunk[:,[0,1]] # Reorder columns
            print(towrite[0:5])
        with open(os.path.join(output_dir, name) + '.cdf', 'ab') as f:
            towrite.tofile(f)
    if normalization:
        return os.path.getsize(os.path.join(output_dir, name) + '.cdf') // 20
    else:
        return os.path.getsize(os.path.join(output_dir, name) + '.cdf') // 16
    #return os.path.getsize(os.path.join(output_dir, name) + '.cdf') // 12

# make a cdh file by specifying the data file and hscan file
def make_cdh(name, num_conic, duration=10, start_time=0):
    
    with open(os.path.join(output_dir, f'{name}.cdh'), 'w') as f:
        f.write('Data filename: ' + name + '.cdf\n' +
        'Number of events: '+  str(num_conic) + '\n' +
        'Data mode: list-mode\n' +
        'Data type: PET\n' +
        'Start time (s):' + str(start_time) + '\n' +
        'Duration (s): ' + str(duration) + '\n' +
        'Scanner name: TPPT\n' +
        'TOF information flag: 1\n' +
        'TOF resolution (ps): 500\n' +
        'List TOF measurement range (ps): 20000.\n'
        'Attenuation correction flag: 0\n'
        'Normalization correction flag: 0\n'
        'Scatter correction flag: 0\n'
        'Random correction flag: 0')

# run the castor recon command with the relevent address and info
def castor_recon(coinc_file_name, post, psf):
    os.chdir(output_dir)
    start = time.time()
    subprocess.run(['castor-recon',
                    '-df', output_dir + coinc_file_name + '.cdh',
                    '-dout', coinc_file_name,
                    '-vb' , '2',
                    '-it', '3:16', 
                    '-conv', psf + '::psf', 
                    '-conv', post + '::post', 
                    '-vox', f'{vox[0]},{vox[1]},{vox[2]}', 
                    '-dim', f'{dim[0]},{dim[1]},{dim[2]}',
                    #'-sens', sens_file,
                    #'-atn', atn_file
                    # 'opt'
                    #'-frm', frames,
                    ])
    print('Time for reconstruction = ' + str(time.time() - start))

def make_frames(frame_duration, total_duration):
    frames = f'0:{frame_duration}'
    for i in range(1, total_duration // frame_duration):
        frames += ',' + str(i * frame_duration)
    frames += f':{total_duration}'
    return frames

frames = make_frames(frame_duration, total_duration)
print(frames)

geometry_to_lut()
write_hscan()
try:
    os.mkdir(output_dir)
except FileExistsError:
    print("The directory for storing figures for this run already exists")

num_coinc = make_cdf(name, num_cols)
make_cdh(name, num_coinc, duration = total_duration, start_time = 0)
#castor_recon(name, f'gaussian,2,2,3',f'gaussian,2,2,3')
castor_recon(name, f'gaussian,2,2,1.5',f'gaussian,2,2,1.5')