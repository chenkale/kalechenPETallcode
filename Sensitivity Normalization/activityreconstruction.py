
# Reconstruction viewer

from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit

name1 = 'hcenter_pcut_a'
name2 = 'mcenter_pcut_a'

img_path = f'/home/kale-chen/Documents/CASToR/images/'
iteration = 1
save = False
out_path = f'/home/kale-chen/Documents/PET/Spatial Resolution/images/'

dim = (10, 10, 50)
vox = (2, 2, 2)

img1 = np.fromfile(os.path.join(img_path, name1, f'{name1}_it{iteration}.img'), dtype = np.float32).reshape(dim, order = 'F')
img2 = np.fromfile(os.path.join(img_path, name2, f'{name2}_it{iteration}.img'), dtype = np.float32).reshape(dim, order = 'F')

# Look at coronal view (collapse y)
img1 = img1.sum(axis = 1)
img2 = img2.sum(axis = 1)

extent = [-dim[0] * vox[0] / 2, dim[0] * vox[0] / 2, -dim[2] * vox[2] / 2, dim[2] * vox[2] / 2]

fig, subplots = plt.subplots(1, 2, figsize = (10, 10))
subplots[0].imshow(img1.T, cmap = 'jet', extent = extent)
subplots[0].set_title(f'')
subplots[0].set_xlabel('x (mm)')
subplots[0].set_ylabel('y (mm)')
subplots[1].imshow(img2.T, cmap = 'jet', extent = extent)
subplots[1].set_title(f'')
subplots[1].set_xlabel('x (mm)')
subplots[1].set_ylabel('y (mm)')
plt.subplots_adjust(hspace = 0)
plt.show()
