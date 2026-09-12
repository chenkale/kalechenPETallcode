# Castor Streamlined

'''
Paths need to be specified
'''

data_file = '/home/kale-chen/Documents/PET/MDA_04112026/R.dat'
name = 'run4'
geometry_file = '/home/kale-chen/Documents/PET/TPPT_Scanner_map_adjusted.csv'
config_dir = '/home/kale-chen/Documents/CASToR/castor_v3.2/config/scanner'
castor_dir = '/home/kale-chen/Documents/CASToR/castor_v3.2/build'
output_dir = '/home/kale-chen/Documents/CASToR/images/'


cut_on_photopeak = False
#cut_on_photopeak = False
# If cutting on photopeak, need to specify:
photopeakL_dir = '/home/kale-chen/Documents/PET/TimeCalibration/ToUpload/photopeakboundsL.txt'
photopeakR_dir = '/home/kale-chen/Documents/PET/TimeCalibration/ToUpload/photopeakboundsR.txt'


# Max threshold for reconstruction viewing
maxthresh = 0.10 # Scales intensity below maxthresh to 0-1, above maxthresh becomes 1

'''
Code begins here
'''


import os
import subprocess
from tqdm import tqdm
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt




start = time.time()

def geometry_to_lut(): # Write geometry lut
    # Read geometry file as np array
    scanner = np.asarray(pd.read_csv(geometry_file, dtype='float32', delimiter=',', header=None))
    scanner.tofile(os.path.join(config_dir, 'TPPT.lut')) # Write file to scanner config directory

def write_hscan(): # Write hscan file
    # write hscan file to scanner config directory
    with open(os.path.join(config_dir, 'TPPT.hscan'), 'w') as f:
        f.write('scanner name: TPPT\n' +
        'modality:    PET\n' +
        'description: This scanner description is based on an actual scanner, however, this implementation is not supported nor validated by its manufacturer.\n' +
        'scanner radius:    ' + '167.5' + '\n' +
        'number of rings in scanner:    32\n' +
        'number of elements:    6144\n' +
        'number of layers:    1\n' +
        'number of crystals in layer(s):   6144\n' +
        'crystals size depth:    15\n' +
        'crystals size transaxial:    3.00\n' +
        'crystals size axial:    3.00\n' +
        'voxels number transaxial:    256\n' +
        'voxels number axial:    128\n' + 
        'field of view transaxial:    4000\n' +
        'field of view axial:    4000\n' +
        'min angle difference:    0' )

def make_cdf(name):
    if cut_on_photopeak:
        photopeakL = np.genfromtxt(photopeakL_dir, dtype = np.float32)
        photopeakR = np.genfromtxt(photopeakR_dir, dtype = np.float32)
# To binary with dtype int16, keep only needed data: EnergyL, IDL, EnergyR, IDR, TimeDiff
# Performs time difference calculation, IDL and IDR mod conversions
# int16 is precise for IDL, IDR, and TimeDiff;
# EnergyL and EnergyR are multiplied by 100 to keep (essentially) 2 decimal precision, exact precision is not needed for photopeak cuts
    if os.path.exists(os.path.join(output_dir, name) + '.cdf'):
        pass
    else:
        chunk_count = 0
        for chunk in pd.read_csv(data_file, sep = '\t', chunksize = 10000000, usecols = [2, 3, 4, 7, 8, 9], header = None):
            chunk.columns = ['TimeL', 'EnergyL', 'IDL', 'TimeR', 'EnergyR', 'IDR']
            chunk['TimeDiff'] = chunk['TimeL'] - chunk['TimeR']     # Calculate time difference, smaller number than individual times
            chunk = chunk[['EnergyL', 'IDL', 'EnergyR', 'IDR', 'TimeDiff']]
            chunk['IDL'] = (chunk['IDL'] - 131072) % 3072           # Convert IDs to 0-3071 range for ease of handling, convert back when writing tsv
            chunk['IDR'] = chunk['IDR'] % 3072                      
            chunk['EnergyL'] = chunk['EnergyL'] * 100               # Multiply Energy by 100 to allow for storage as int16 with 2 decimal precision
            chunk['EnergyR'] = chunk['EnergyR'] * 100

            if cut_on_photopeak:
                energyL = chunk["EnergyL"].to_numpy(dtype=float)
                energyR = chunk["EnergyR"].to_numpy(dtype=float)
                lowerboundLs = photopeakL[chunk["IDL"], 0]
                upperboundLs = photopeakL[chunk["IDR"], 1]
                lowerboundRs = photopeakR[chunk["IDL"], 0]
                upperboundRs = photopeakR[chunk["IDR"], 1]
                mask = ((energyL > lowerboundLs) & (energyL < upperboundLs)
                & (energyR > lowerboundRs) & (energyR < upperboundRs))
                chunk = chunk[mask]
            chunk = chunk[['IDL', 'IDR', 'TimeDiff']]

            # Process into cdf format: (Time, TOF, IDL, IDR)
            chunk = np.asarray(chunk, dtype = np.int32) 
            chunk[:, 0] = chunk[:, 0] + 3072 # Convert to RecoID
            chunk[:, 2] = chunk[:, 2] // 2 # Time diff to TOF
            towrite = np.zeros((len(chunk), 4), dtype = np.int32) # For reordering
            towrite[:,[1,2,3]] = chunk[:,[2,0,1]] # Reorder columns
            with open(os.path.join(output_dir, name) + '.cdf', 'ab') as f:
                towrite.tofile(f)
            print(f'cdf: {time.time() - start} chunk {chunk_count}')
            chunk_count += 1
    
    print(f'Written cdf: {time.time() - start}')
    return os.path.getsize(os.path.join(output_dir, name) + '.cdf') // 16

# make a cdh file by specifying the data file and hscan file
def make_cdh(name, num_conic, tof=60):
    
    with open(os.path.join(output_dir, f'{name}.cdh'), 'w') as f:
        f.write('Data filename: ' + name + '.cdf\n' +
        'Number of events: '+  str(num_conic) + '\n' +
        'Data mode: list-mode\n' +
        'Data type: PET\n' +
        'Start time (s):0\n' +
        'Duration (s): ' + str(tof) + '\n' +
        'Scanner name: TPPT\n' +
        'TOF information flag: 1\n' +
        'TOF resolution (ps): 20000\n' +
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
                    '-it', '1:16', 
                    '-conv', psf + '::psf', 
                    '-conv', post + '::post', 
                    '-vox', '1,1,1', 
                    '-dim', '256,256,128'
                    # 'opt'
                    # '-frm', frames,
                    ])
    print('Time for reconstruction = ' + str(time.time() - start))


geometry_to_lut()
write_hscan()
try:
    os.mkdir(output_dir)
except FileExistsError:
    print("The directory for storing figures for this run already exists")

num_coinc = make_cdf(name)
make_cdh(name, num_coinc)
castor_recon(name, f'gaussian,1.5,1.5,1.75',f'gaussian,1.5,1.5,1.75')

# Reconstruction viewer

img_path = os.path.join(output_dir, name, f'{name}_it1.img')

maxthresh = 0.10 # Scales intensity below maxthresh to 0-1, above maxthresh becomes 1
dim = (128, 128, 64)
vox = (1.5, 1.5, 1.5)

# Define bounds for use in histogramming/heatmapping
xlower, xupper = -dim[0] * vox[0] / 2, dim[0] * vox[0] / 2
ylower, yupper = -dim[1] * vox[1] / 2, dim[1] * vox[1] / 2
zlower, zupper = -dim[2] * vox[2] / 2, dim[2] * vox[2] / 2

def applythresh(img2d):
    img2d = img2d * (1 / maxthresh)
    for i in range(len(img2d)):
        for j in range(len(img2d[i])):
            if img2d[i][j] > 1:
                img2d[i][j] = 1
    return img2d

img = np.fromfile(img_path, dtype = np.float32).reshape(dim, order = 'F')

# Transverse (xy)
fig = plt.figure(0)
axxyhm = fig.add_subplot() # Heatmap
imgxy = applythresh(img.sum(axis = 2)) # Image with threshold cast
xyhm = axxyhm.imshow(imgxy.T, cmap='viridis', extent=(xlower, xupper, ylower, yupper), vmin = 0, vmax = 1) # Display heatmap
axxyhm.set_xlim(xlower, xupper) 
axxyhm.set_ylim(ylower, yupper)
axxyhm.set_xlabel("x (mm)")
axxyhm.set_ylabel("y (mm)")
axxyhm.set_title('Transverse Sum (collapse z)')

# Coronal (xz)
fig = plt.figure(1)
axxzhm = fig.add_subplot() # Heatmap
imgxz = applythresh(img.sum(axis = 1)) # Image with threshold cast
xzhm = axxzhm.imshow(imgxz.T, cmap='viridis', extent=(xlower, xupper, zlower, zupper), vmin = 0, vmax = 1) # Display heatmap
axxzhm.set_xlim(xlower, xupper) 
axxzhm.set_ylim(zlower, zupper)
axxzhm.set_xlabel("x (mm)")
axxzhm.set_ylabel("z (mm)")
axxyhm.set_title('Transverse Sum (collapse y)')

# Sagittal (yz)
fig = plt.figure(2)
axyzhm = fig.add_subplot() # Heatmap
imgyz = applythresh(img.sum(axis = 0)) # Image with threshold cast
yzhm = axyzhm.imshow(imgyz.T, cmap='viridis', extent=(ylower, yupper, zlower, zupper), vmin = 0, vmax = 1) # Display heatmap
axyzhm.set_xlim(ylower, yupper) 
axyzhm.set_ylim(zlower, zupper)
axyzhm.set_xlabel("y (mm)")
axyzhm.set_ylabel("z (mm)")
axxyhm.set_title('Transverse Sum (collapse x)')

plt.show()
