# Castor Streamlined

'''
Paths need to be specified
'''
#names = ['0932','0923','0916','1732','1723','1716','2532','2523','2516']
#names = ['hcenterpgd','mcenterpgd']
names = ['pr0', 'pR0', 'pR1', 'pR2', 'pR3', 'pR4', 'pR5', 'pR6', 'pR7']
#names = ['pra0', 'pRa0', 'pRa1', 'pRa2', 'pRa3', 'pRa4', 'pRa5', 'pRa6', 'pRa7']
#names = ['tofpair1', 'tofpair2', 'tofpair3', 'tofpair4', 'tofpair5']
#names = ['htofpgd0', 'htofpgd1', 'htofpgd2', 'htofpgd3', 'htofpgd4', 'htofpgd5']
#names = ['empty']
#names = ['hcenterpgd']
#names = ['hcircpgd0','hcircpgd1','hcircpgd2','hcircpgd3','hcircpgd4','hcircpgd5','hcircpgd6','hcircpgd7']
out_names = [''] * len(names) # Use to specify output name, if left blank, will use name
tag = '1ns'
#centers = [(-78, 73, 0),(3, 76, 0),(75, 73, 0),(-78, 0, 0),(4, 0, 0),(73, -2, 0),(-77, -78, 0),(4, -78, 0),(74, -76, 0)]
centers = [(-2, 11, 0),(-3, 96, 0),(-71, 68, 0),(-99, 0, 0),(-71, -69, 0),(-3, -97, 0),(65,-70, 0),(93, -2, 0),(66, 67, 0)]
#centers = [(40, 0, 0), (40, 0, 0), (40, 0, 0), (40, 0, 0), (40, 0, 0), (40, 0, 0)]
#centers = [(0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]
#centers = [(0, 120, 0), (-85, 85, 0), (-120, 0, 0), (-85, -85, 0), (0, -120, 0), (85, -85, 0), (120, 0, 0), (85, 85, 0)]
#centers = [(120, 50, 0)]
#data_dir = '/home/kale-chen/Documents/PET/Spatial Resolution/Data2/'
data_dir = '/home/kale-chen/Documents/PET/TPPT2026/data/'

# Flags and processing
autoTag = True
if autoTag:
    normalization = True if 'n' in tag else False
    tof = True if 't' in tag else False
    time_alignment = True if 'c' in tag else False
    photopeak_cut = True if 'p' in tag else False
    usesens = True if 's' in tag else False
    projector = 'distanceDriven' if 'd' in tag else 'joseph'
else:
    normalization = False
    tof = False
    time_alignment = False
    photopeak_cut = False
    usesens = False
geometry_file = '/home/kale-chen/Documents/PET/TPPT_Scanner_map_adjusted.csv'
#geometry_file = '/home/kale-chen/Documents/PET/TPPT_Scanner_map2.csv'
config_dir = '/home/kale-chen/Documents/CASToR/castor_v3.2/config/scanner'
castor_dir = '/home/kale-chen/Documents/CASToR/castor_v3.2/build'
output_dir = '/home/kale-chen/Documents/CASToR/images/'


# Included information 
normalization_file = '/home/kale-chen/Documents/PET/Sensitivity Normalization/normfactors.bin'
norm_cdh = '/home/kale-chen/Documents/PET/Sensitivity Normalization/normfactors.cdh'
time_offsets_file = '/home/kale-chen/Documents/PET/TPPT2026/data/tsvs/time_offset_calibrationit10.tsv'
tof_resolution = 520
sens_file = '/home/kale-chen/Documents/PET/Sensitivity Normalization/sensimgs/correctedimg.img'
#sens_file = '/home/kale-chen/Documents/PET/Sensitivity Normalization/sensimgs/sensjosephbipartite.img'
atn_file = '/home/kale-chen/Documents/CASToR/images/mumap/mumap512.hdr'
maskmodules = False
# castor-recon parameters
#dim = (64, 64, 64)
#vox = (0.5, 0.5, 0.5)
#dim = (120, 15, 100)
#vox = (1, 1, 1)
#dim = (20, 20, 100)
#vox = (1, 1, 1)
#dim = (10, 10, 50)
#vox = (2, 2, 2)
dim = (360, 360, 110)
#vox = (1, 1, 1)
#dim = (180,180, 55)
#vox = (2, 2, 2)
#dim = (100, 100, 50)
vox = (1, 1, 1)
if '2' in tag:
    dim = (10, 10, 50)
    vox = (2, 2, 2)
elif '1' in tag:
    dim = (32, 32, 100)
    vox = (1, 1, 1)

psf = 'gaussian,0,0,0'
post = 'gaussian,0,0,0'
#psf = 'gaussian,1.5,1.5,3'
#post = 'gaussian,1.5,1.5,3'

total_duration = 300
frame_duration = 30

num_iterations = 3

# Preprocessing
time_offsets_file = '/home/kale-chen/Documents/PET/TimeCalibration/CalibrationDataPreparer/tsvs/time_offset_calibrationit10.tsv'

# Low statistics studies
stride = 1

# For now define guess for number of columns
guess = 4

'''
Code begins here
'''

import os
import subprocess
from tqdm import tqdm
import numpy as np
import pandas as pd
import time
from tqdm import tqdm

if normalization:
    normalization_data = np.fromfile(normalization_file, dtype = np.float32).reshape(3072, 3072)
if time_alignment:
    time_offsets = np.genfromtxt(time_offsets_file, delimiter = '\t')[:, 4]
if maskmodules:
    modulemask = np.fromfile('/home/kale-chen/Documents/PET/mask.bin', dtype = np.bool)

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
        'voxels number transaxial:    2projector56\n' +
        'voxels number axial:    128\n' + 
        'field of view transaxial:    4000\n' +
        'field of view axial:    4000\n' +
        'min angle difference:    0' )

def guess_num_cols(data_file):
    global guess

    return guess
    guess = 2
    guess_found = False
    data = np.memmap(data_file, dtype = np.int16, mode = 'r', shape = (1000000, guess))
    while not guess_found:
        chunk = data[:(1000000 // guess * guess)]
        vals = np.unique(chunk[:, 0])
        if len(vals) < 3072:
            guess_found = True
    return guess
if maskmodules:
    modulestokeep = np.array([9, 10, 13, 14])
    channelstokeep = np.array([n * 128 + m for n in modulestokeep for m in range(128)])
    channelstokeep = np.zeros(6144, dtype = np.bool)
    channelstokeep = channelstokeep + modulemask
else:
    channelstokeep = np.zeros(6144, dtype = np.bool)
def process_chunk(chunk):
    if num_cols == 6:
        chunk = chunk[:, [1, 3, 4]]
    elif num_cols == 4:
        chunk = chunk[:, [0, 1, 2]]
    newchunk = np.zeros((len(chunk), 3), dtype = np.int16)
    if time_alignment:
        writehead = 0
        for row in chunk:
            newchunk[writehead] = np.array([row[0], row[1], (row[2] - time_offsets[row[0] + 3072] + time_offsets[row[1]])], dtype = np.int16)
            writehead += 1
    else:
        newchunk = chunk
    if maskmodules:
        #newchunk = newchunk[np.isin(newchunk[:, 0], channelstokeep) & np.isin(newchunk[:, 1], channelstokeep)]
        newchunk = newchunk[modulemask[newchunk[:, 0] + 3072] & modulemask[newchunk[:, 1]]]
    return newchunk

def make_cdf(data_file, out_name, num_cols):
    if os.path.exists(os.path.join(output_dir, out_name) + '.cdf'):
        os.remove(os.path.join(output_dir, out_name) + '.cdf')
    num_rows = os.path.getsize(data_file) // (num_cols * 2)
    data = np.memmap(data_file, dtype = np.int16, mode = 'r', shape = (num_rows, num_cols))
    eventcount = 0
    for i in tqdm(range(0, num_rows, 1000000)):
        chunk = np.copy(data[i:i+1000000:stride])
        chunk = process_chunk(chunk)
        if normalization and tof:
            # Structured array towrite has [time (uint32), norm (float32), tof (float32), idl (uint32), idr (uint32)]
            towrite = np.zeros((len(chunk)), dtype = [('time', 'uint32'), ('norm', 'float32'), ('tof', 'float32'), ('idl', 'uint32'), ('idr', 'uint32')]) # For reordering
            norm = np.array(normalization_data[chunk[:, 0], chunk[:, 1]], dtype = np.float32)
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            towrite['norm'] = norm
            towrite['tof'] = chunk[:, 2].astype(np.float32)
            towrite['idl'] = chunk[:, 0].astype(np.uint32)
            towrite['idr'] = chunk[:, 1].astype(np.uint32)
        elif normalization and not tof:
            # Structured array towrite has [time (uint32), norm (float32), idl (uint32), idr (uint32)]
            towrite = np.zeros((len(chunk)), dtype = [('time', 'uint32'), ('norm', 'float32'), ('idl', 'uint32'), ('idr', 'uint32')])
            norm = np.array(normalization_data[chunk[:, 0], chunk[:, 1]], dtype = np.float32)
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            towrite['norm'] = norm
            towrite['idl'] = chunk[:, 0].astype(np.uint32)
            towrite['idr'] = chunk[:, 1].astype(np.uint32)
        elif not normalization and tof:
            # Structured array towrite has [time (uint32), tof (float32), idl (uint32), idr (uint32)]
            towrite = np.zeros((len(chunk)), dtype = [('time', 'uint32'), ('tof', 'float32'), ('idl', 'uint32'), ('idr', 'uint32')])
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            towrite['tof'] = chunk[:, 2].astype(np.float32)
            towrite['idl'] = chunk[:, 0].astype(np.uint32)
            towrite['idr'] = chunk[:, 1].astype(np.uint32)
        else:
            # Structured array towrite has [time (uint32), idl (uint32), idr (uint32)]
            towrite = np.zeros((len(chunk)), dtype = [('time', 'uint32'), ('idl', 'uint32'), ('idr', 'uint32')])
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            towrite['idl'] = chunk[:, 0].astype(np.uint32)
            towrite['idr'] = chunk[:, 1].astype(np.uint32)
        with open(os.path.join(output_dir, out_name) + '.cdf', 'ab') as f:
            towrite.tofile(f)
        eventcount += len(chunk)
    return eventcount

# make a cdh file by specifying the data file and hscan file
def make_cdh(num_rows, out_name, duration=10000, start_time=0):
    
    with open(os.path.join(output_dir, out_name + '.cdh'), 'w') as f:
        f.write('Data filename: ' + out_name + '.cdf\n' +
        'Number of events: '+  str(num_rows) + '\n' +
        'Data mode: list-mode\n' +
        'Data type: PET\n' +
        'Start time (s):' + str(start_time) + '\n' +
        'Duration (s): ' + str(duration) + '\n' +
        'Scanner name: TPPT\n' +
        f'TOF information flag: {1 if tof else 0}\n' +
        f'TOF resolution (ps): {tof_resolution}\n' +
        'List TOF measurement range (ps): 20000.\n'
        'Attenuation correction flag: 0\n'
        f'Normalization correction flag: {1 if normalization else 0}\n'
        'Scatter correction flag: 0\n'
        'Random correction flag: 0')

def make_sens_img(out_name, center):
    if not os.path.exists(os.path.join(output_dir, out_name)):
        os.mkdir(os.path.join(output_dir, out_name))
    if os.path.exists(os.path.join(output_dir, out_name, f'{out_name}_sensitivity.img')):
        os.remove(os.path.join(output_dir, out_name, f'{out_name}_sensitivity.img'))
    sensimg = np.fromfile(sens_file, dtype = np.float32).reshape(360, 360, 110, order = 'F')
    #sensimg = np.flip(np.flip(sensimg, axis = 1), axis = 2)
    region = sensimg[180 + center[0]-dim[0]//2:180 + center[0] + dim[0]//2, 180 + center[1] - dim[1]//2:180 + center[1] + dim[1]//2, 55 + center[2] - dim[2]//2:55 + center[2] + dim[2]//2]
    #region = np.flip(np.flip(region, axis = 1), axis = 2)
    region.astype(np.float32).ravel(order='F').tofile(os.path.join(output_dir, out_name, f'{out_name}_sensitivity.img'))

def make_sens_hdr(out_name):
    with open(os.path.join(output_dir, out_name, out_name + '_sensitivity.hdr'), 'w') as f:
        f.write('!INTERFILE := \n' +    
        '!imaging modality := PET\n' +
        '!version of keys := CASToRv3.2.1\n' +
        'CASToR version := 3.2.1\n\n' +
        '!GENERAL DATA := \n' +
        '!originating system := TPPT\n' +
        '!data offset in bytes := 0\n' +
        f'!name of data file := {out_name}_sensitivity.img\n' +
        f'patient name := {out_name}_sensitivity\n\n' +
        '!GENERAL IMAGE DATA \n' +
        '!type of data := Dynamic\n' +
        f'!total number of images := {dim[2]}\n' +
        'imagedata byte order := LITTLEENDIAN\n' +
        '!number of frame groups :=1 \n' +
        'process status := \n\n' +
        '!STATIC STUDY (General) := \n' +
        'number of dimensions := 3\n' +
        f'!matrix size [1] := {dim[0]}\n' +
        f'!matrix size [2] := {dim[1]}\n' +
        f'!matrix size [3] := {dim[2]}\n' +
        '!number format := short float\n' +
        '!number of bytes per pixel := 4\n' +
        f'scaling factor (mm/pixel) [1] := {vox[0]}\n' +
        f'scaling factor (mm/pixel) [2] := {vox[1]}\n' +
        f'scaling factor (mm/pixel) [3] := {vox[2]}\n' +
        'first pixel offset (mm) [1] := 0\n' +
        'first pixel offset (mm) [2] := 0\n' +
        'first pixel offset (mm) [3] := 0\n' +
        'data rescale offset := 0\n' +
        'data rescale slope := 1\n' +
        'quantification units := 1\n' +
        f'!number of images in this frame group := {dim[2]}\n' +
        f'!image duration (sec) := {total_duration}\n' +
        '!image start time (sec) := 0\n' +
        'pause between frame groups (sec) := 0\n' +
        '!END OF INTERFILE := \n\n' +
        '!COPY OF INPUT HEADER 1\n' +
        f'Data filename: {out_name}.cdf\n' +
        f'Number of events: {num_rows}\n' +
        'Data mode: list-mode\n' +
        'Data type: PET\n' +
        f'Start time (s):{0}\n' +
        f'Duration (s): {total_duration}\n' +
        'Scanner name: TPPT\n' +
        f'TOF information flag: {1 if tof else 0}\n' +
        f'TOF resolution (ps): {tof_resolution}\n' +
        f'List TOF measurement range (ps): {20000}.\n' +
        'Attenuation correction flag: 0\n' +
        f'Normalization correction flag: {1 if normalization else 0}\n' +
        'Scatter correction flag: 0\n' +
        '!END OF COPY OF INPUT HEADER 1')

    return os.path.join(output_dir, out_name, out_name + '_sensitivity.hdr')

# run the castor recon command with the relevent address and info
def castor_recon(coinc_file_name, center, senshdr):
    os.chdir(output_dir)
    subprocess.run(['castor-recon',
                    '-df', os.path.join(output_dir, coinc_file_name + '.cdh'),
                    '-dout', coinc_file_name,
                    '-vb' , '2',
                    '-it', f'{num_iterations}:16', 
                    #'-conv', psf + '::psf', 
                    #'-conv', post + '::post', 
                    '-vox', f'{vox[0]},{vox[1]},{vox[2]}', 
                    '-dim', f'{dim[0]},{dim[1]},{dim[2]}',
                    '-th', '4',
                    *(['-sens', senshdr] if usesens else []),
                    #'-atn', atn_file
                    # 'opt'
                    #'-frm', frames,
                    *(["-norm", norm_cdh] if normalization else []),
                    '-off', f'{center[0]},{center[1]},{center[2]}',
                    '-proj', projector
                    ])

def make_frames(frame_duration, total_duration):
    frames = f'0:{frame_duration}'
    for i in range(1, total_duration // frame_duration):
        frames += ',' + str(i * frame_duration)
    frames += f':{total_duration}'
    return frames


#frames = make_frames(frame_duration, total_duration)
#print(frames)


geometry_to_lut()
write_hscan()
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass
start = time.time()
times = []
for name, out_name, center in zip(names, out_names, centers):
    realcenter = (center[0] * -1, center[1] * -1, center[2])
    data_file = os.path.join(data_dir, name + '.bin')
    if out_name == '':
        out_name = name + tag
    num_cols = guess_num_cols(data_file)
    num_rows = make_cdf(data_file, out_name, num_cols)
    make_cdh(num_rows, out_name, duration = total_duration, start_time = 0)
    make_sens_img(out_name, center)
    senshdr = make_sens_hdr(out_name)
    castor_recon(out_name, realcenter, senshdr = None if not usesens else senshdr)
    times.append(round(time.time() - start, 1))
print(f'Times = {times}')

