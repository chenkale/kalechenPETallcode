
# Reconstruction viewer

import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
import SimpleITK as sitk

name = 'hcenterpsi1'
#name = 'run3_pcut_tcut'

img_path = f'/home/kale-chen/Documents/CASToR/images/'
iteration = 1
save = False
transverseHistogram = True
fit_derenzogaussian = False
saggitalHistogram = False
bilateral_smoothing = False

# Thresholding adjusts maximum intensity to boost low intensity voxels and cap high intensity voxels
auto_threshold = True # If true maxthresh is not used
maxthresh = 0.0001 # Scales intensity below maxthresh to 0-1, above maxthresh becomes 1


#dim = (64, 64, 128)
#dim = (512, 512, 256)
#dim = (256, 256, 128)
dim = (20, 20, 100)
vox = (1, 1, 1)
#vox = (0.5, 0.5, 0.5)
#vox = (1.5, 1.5, 1.5)

trim = (8, 8, 8) # In voxels
trim = (1, 1, 1)
#trim = (24, 24, 16) # In voxels

img = np.fromfile(os.path.join(img_path, name, f'{name}_it{iteration}.img'), dtype = np.float32).reshape(dim, order = 'F')
#img = np.fromfile(os.path.join(img_path, name, f'{name}_sensitivity.img'), dtype = np.float32).reshape(dim, order = 'F')
#img = np.fromfile(os.path.join(img_path, name, f'mumap512.img'), dtype = np.float32).reshape(dim, order = 'F')
cmap = 'jet'
#cmap = 'gray' # Gray scale



'''
Shouldn't need to change anything below here
'''

# Try to fit this
def derenzogaussian(x, a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, a5, b5, c5, bg):
    return a1 * np.exp(-(x - b1)**2 / (2 * c1**2)) + \
    a2 * np.exp(-(x - b2)**2 / (2 * c2**2)) + \
    a3 * np.exp(-(x - b3)**2 / (2 * c3**2)) + \
    a4 * np.exp(-(x - b4)**2 / (2 * c4**2)) + \
    a5 * np.exp(-(x - b5)**2 / (2 * c5**2)) + \
    bg

def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (2 * c**2))

# Flip y axis, this is needed for image to be displayed correctly (specifically xy)
img = np.flip(img, axis=1)
img = np.flip(img, axis=2)

# Voxel intensity histogram, for use in auto thresholding
fig = plt.figure(3, figsize = (10, 5))
ax = fig.add_subplot()
ax.hist(img[trim[0]:-trim[0], trim[1]:-trim[1], trim[2]:-trim[2]].flatten(), bins = 100, histtype = 'step', color = 'blue')
if auto_threshold:
    auto_thresh = np.percentile(img.flatten(), 99.99)
    ax.vlines(auto_thresh, ymin = 0, ymax = np.max(plt.gca().get_yticks()), color = 'red', linestyles = 'dashed', label = f'Auto threshold: {auto_thresh:.2e}')
    maxthresh = auto_thresh
else:
    ax.vlines(maxthresh, ymin = 0, ymax = np.max(plt.gca().get_yticks()), color = 'red', linestyles = 'dashed', label = f'Set threshold: {maxthresh:.2e}')
ax.set_xlabel('Intensity')
ax.set_ylabel('Frequency')
ax.set_title('Voxel intensities')
ax.set_yscale('log')
ax.legend()
xticks = ax.get_xticks()
xlabels = [f"{x:.1e}" for x in xticks]
ax.set_xticklabels(xlabels)
if save:
    plt.savefig(os.path.join(img_path, name, f'{name}_intensities_it{iteration}.png'))

if bilateral_smoothing:
    img = sitk.GetImageFromArray(img)
    
    img_filtered = sitk.Bilateral(
        img,
        domainSigma=2.0,   # spatial (in voxels unless spacing set)
        rangeSigma=0.01
    )
    img = sitk.GetArrayFromImage(img_filtered)

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


# Transverse (xy) - Add histogram below heatmap, x-axis aligned
if transverseHistogram:
    fig = plt.figure(0, figsize=(5, 6))
    gs = GridSpec(2, 1, height_ratios=[5, 1], hspace=0)
    axxyhm = fig.add_subplot(gs[0])   # Heatmap
    axxyhm.grid(True, linewidth=0.5)
    axhist = fig.add_subplot(gs[1], sharex=axxyhm) # Histogram (share x axis)
else:
    fig = plt.figure(0)
    axxyhm = fig.add_subplot()

imgxy = applythresh(img.sum(axis=2) / dim[2])  # Image with threshold cast
imgxy = imgxy[trim[0]:-trim[0], trim[1]:-trim[1]]

imgx = imgxy.sum(axis=1)  # Sum over y, gives x-projection

# Heatmap
if transverseHistogram:
    axxyhm.set_xlim(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0])
    axxyhm.set_ylim(ylower + trim[1]*vox[1], yupper - trim[1]*vox[1])
    xyhm = axxyhm.imshow(
    imgxy.T / np.max(imgxy), 
    cmap=cmap, 
    extent=(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0], ylower + trim[1]*vox[1], yupper - trim[1]*vox[1]),
    vmin=0, vmax=1, 
    aspect='auto'
)
else:
    axxyhm.set_xlim(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0])
    axxyhm.set_ylim(ylower + trim[1]*vox[1], yupper - trim[1]*vox[1])
    xyhm = axxyhm.imshow(
    imgxy.T/ np.max(imgxy), 
    cmap=cmap, 
    extent=(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0], ylower + trim[1]*vox[1], yupper - trim[1]*vox[1]),
    vmin=0, vmax=1)
axxyhm.set_ylabel("y (mm)")
axxyhm.set_title(f'Transverse Sum (collapse z) it {iteration}')
#axxyhm.set_title('Sensitivity map selected channels')

# Histogram under heatmap
# Bin edges for x
if transverseHistogram:
    x_bins = np.linspace(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0], dim[0]+1-2*trim[0])
    bin_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    fit_success = False
    try:
        if fit_derenzogaussian:
            p0 = [np.max(imgx)*0.6, 0, 2,
                  np.max(imgx)*0.4, -5, 2,
                  np.max(imgx)*0.4, +5, 2,
                  np.max(imgx)*0.2, -10, 2,
                  np.max(imgx)*0.2, +10, 2,
                  np.max(imgx)*0.05]
            popt, pcov = curve_fit(derenzogaussian, bin_centers, imgx, p0 = p0)
        else:
            popt, pcov = curve_fit(gaussian, bin_centers, imgx)
        axhist.plot(
            np.linspace(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0], 100), 
            gaussian(np.linspace(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0], 100), *popt), 
            color='red',
            linewidth=0.5
        )
        fit_success = True
    except:
        print('Fit failed')
        popt = None

    axhist.step(bin_centers, imgx, color='blue', where='mid')
    axhist.set_xlim(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0])
    axhist.set_xlabel("x (mm)")
    
    plt.setp(axhist.get_yticklabels(), visible=False)
    #axhist.legend(['x projection', f'μ={popt[1]:.2f}, σ={abs(popt[2]):.2f}'])
    if fit_success:
        axhist.legend([f'μ={popt[1]:.3f}, σ={abs(popt[2]):.3f}'])
    else:
        pass
    axhist.tick_params(axis='y', which='both', left=False, right=False)
    plt.setp(axxyhm.get_xticklabels(), visible=False)  # Hide upper plot ticks (for cleanliness)

if save:
    plt.savefig(os.path.join(img_path, name, f'{name}_xy_it{iteration}.png'))

# Coronal (xz)
fig = plt.figure(1)
axxzhm = fig.add_subplot() # Heatmap
imgxz = applythresh(img.sum(axis = 1) / dim[1]) # Image with threshold cast
imgxz = imgxz[trim[0]:-trim[0], trim[2]:-trim[2]]
xzhm = axxzhm.imshow(imgxz.T / np.max(imgxz), cmap=cmap, extent=(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0], zlower + trim[2]*vox[2], zupper - trim[2]*vox[2]), vmin = 0, vmax = 1) # Display heatmap
axxzhm.set_aspect('equal')
axxzhm.set_xlim(xlower + trim[0]*vox[0], xupper - trim[0]*vox[0]) 
axxzhm.set_ylim(zlower + trim[2]*vox[2], zupper - trim[2]*vox[2])
axxzhm.set_xlabel("x (mm)")
axxzhm.set_ylabel("z (mm)")
axxzhm.set_title('Coronal Sum (collapse y)')
if save:
    plt.savefig(os.path.join(img_path, name, f'{name}_xz_it{iteration}.png'))

# Sagittal (yz)
if saggitalHistogram:
    fig = plt.figure(2, figsize = (6, 5))
    gs = GridSpec(1, 2, width_ratios=[5, 1], hspace=0)
    axyzhm = fig.add_subplot(gs[0])
    axhist = fig.add_subplot(gs[1], sharey=axyzhm)
else:
    fig = plt.figure(2)
    axyzhm = fig.add_subplot()

imgyz = applythresh(img.sum(axis=0) / dim[0])
imgyz = imgyz[trim[1]:-trim[1], trim[2]:-trim[2]]
imgz = imgyz.sum(axis=0)


# Display heatmap for sagittal projection
yzhm = axyzhm.imshow(
    imgyz.T / np.max(imgyz),
    cmap=cmap,
    extent=(ylower + trim[1]*vox[1], yupper - trim[1]*vox[1], zlower + trim[2]*vox[2], zupper - trim[2]*vox[2]),
    vmin=0, vmax=1
)
axyzhm.set_xlim(ylower+trim[1]*vox[1], yupper-trim[1]*vox[1])
axyzhm.set_ylim(zlower+trim[2]*vox[2], zupper-trim[2]*vox[2])
axyzhm.set_aspect('equal')
axyzhm.set_xlabel("y (mm)")
axyzhm.set_ylabel("z (mm)")
axyzhm.set_title('Sagittal Sum (collapse x)')

if saggitalHistogram:
    z_bins = np.linspace(zlower + trim[2]*vox[2], zupper - trim[2]*vox[2], dim[2]+1)
    z_bin_centers = 0.5 * (z_bins[:-1] + z_bins[1:])
    # Step 1: Prepare the bar heights and bin centers for the horizontal bar plot
    bar_heights = imgz[::-1]
    bar_positions = z_bin_centers
    bar_height_size = z_bins[1] - z_bins[0]

    # Step 2: Plot the horizontal bars with the specified appearance
    axhist.stairs(bar_heights, z_bins, color='blue', alpha=0.6, orientation='horizontal')


    axhist.set_ylim(zlower + trim[2]*vox[2], zupper - trim[2]*vox[2])
    plt.setp(axhist.get_xticklabels(), visible=False)

if save:
    plt.savefig(os.path.join(img_path, name, f'{name}_yz_it{iteration}.png'))

if not save:
    plt.show()

