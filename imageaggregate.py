from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec

names = ['0932',
         '0923',
         #'0916',
         '1732',
         '1723',
         '1716',
         '2532',
         '2523',
         '2516']
it = 1
#img_path = '/home/kale-chen/Documents/CASToR/images/'
img_path = '/home/kale-chen/Documents/PET/Spatial Resolution/Data2/'


dim = (201, 201, 1)
vox = (160 * 2 / 200, 160 * 2 / 200, 1)
#cmap = 'jet'
cmap = 'hot'

aggregate = np.zeros(dim)
for name in names:
    #img = np.fromfile(os.path.join(img_path, name, f'{name}_it{it}.img'), dtype = np.float32).reshape(dim, order = 'F')
    img = np.fromfile(os.path.join(img_path, f'{name}_fbp.bin'), dtype = np.float64).reshape(dim, order = 'F')
    aggregate += img
aggregate = np.flip(aggregate, axis=2)

aggregate = aggregate.sum(axis=2).T

plt.figure(0, figsize=(5,5))
plt.imshow(aggregate, cmap=cmap, aspect='auto', extent=(-dim[0]/2 * vox[0], dim[0]/2 * vox[0], -dim[1]/2 * vox[1], dim[1]/2 * vox[1]))
plt.xticks(np.arange(-90, 90.00001, 30))
plt.yticks(np.arange(-90, 90.00001, 30))
plt.xlim(-96, 96)
plt.ylim(-96, 96)
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.show()
