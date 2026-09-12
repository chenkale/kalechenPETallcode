
# Reconstruction viewer

from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit

#names = ['0932',
#         '0923',
#         '0916',
#         '1732',
#         '1723',
#         '1716',
#         '2532',
#         '2523',
#         '2516']
names = ['pr0', 'pR0', 'pR1', 'pR2', 'pR3', 'pR4', 'pR5', 'pR6', 'pR7']
#names = ['pra0', 'pRa0', 'pRa1', 'pRa2', 'pRa3', 'pRa4', 'pRa5', 'pRa6', 'pRa7']
#names = ['hcircpgd0','hcircpgd1','hcircpgd2','hcircpgd3','hcircpgd4','hcircpgd5','hcircpgd6','hcircpgd7']

tag = '1ns'
#centers = [(-78, 73, 0),
#          (3, 76, 0),
#          (75, 73, 0),
#          (-78, 0, 0),
#          (4, 0, 0),    
#          (73, -2, 0),
#          (-77, -78, 0),
#          (4, -78, 0),
#          (74, -76, 0)]
centers = [(-2, 11, 0), 
        (-3, 96, 0), 
        (-71, 68, 0), 
        (-99, 0, 0), 
        (-71, -69, 0),
        (-3, -97, 0),
        (65,-70, 0),
        (93, -2, 0),
        (66, 67, 0)
]
#centers = [(0, 120, 0), (-85, 85, 0), (-120, 0, 0), (-85, -85, 0), (0, -120, 0), (85, -85, 0), (120, 0, 0), (85, 85, 0)]


extraname = ''

img_path = f'/home/kale-chen/Documents/CASToR/images/'
#iterations = range(1, 31)
iterations = range(1, 2)

save = True
fit = True
out_path = f'/home/kale-chen/Documents/PET/Sensitivity Normalization/writeup/plots/spatialresolution/'
csv_path = f'/home/kale-chen/Documents/PET/Spatial Resolution/Data2/resolutionsiterated.csv'
mode = 'MLEM'

#dim = (64, 64, 64)
#vox = (0.5, 0.5, 0.5)
if '2' in tag:
    dim = (10, 10, 50)
    vox = (2, 2, 2)
elif '1' in tag:
    dim = (32, 32, 100)
    vox = (1, 1, 1)
#wsize = (dim[0] - 2, dim[1] - 2, int(5 / vox[2])) # Window size in voxels, centered on viewcenter
wsize = (dim[0], dim[1], 10)

tick_interval = wsize[0]*vox[0] / 8 # Tick interval in mm

cmap = 'jet'

hmapsums = []


def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (2 * c**2))
def extendEdges(array):
    array = np.append(array, array[-1])
    array = np.insert(array, 0, array[0])
    return array
vmax = 0
sumtrim = 0
for iteration in iterations:
    for i in range(len(names)):
        name = names[i]
        img = np.fromfile(os.path.join(img_path, name+tag, f'{name}{tag}_it{iteration}.img'), dtype = np.float32).reshape(dim, order = 'F')
        pixextentx = (int(dim[0] / 2 - wsize[0] / 2), int(dim[0] / 2 + wsize[0] / 2))
        pixextenty = (int(dim[1] / 2 - wsize[1] / 2), int(dim[1] / 2 + wsize[1] / 2))
        pixextentz = (int(dim[2] / 2 - wsize[2] / 2), int(dim[2] / 2 + wsize[2] / 2))
        print(pixextentx, pixextenty, pixextentz)
        img = img[pixextentx[0]:pixextentx[1], pixextenty[0]:pixextenty[1], pixextentz[0]:pixextentz[1]]
        print(img.shape)
        #img = np.flip(img, axis=0)
        img = np.flip(img, axis=1)
        img = np.flip(img, axis=2)

        imgxy = img.sum(axis=2) 
        imgx = imgxy.sum(axis=1)  # Sum over y, gives x-projection
        imgy = imgxy.sum(axis=0)  # Sum over x, gives y-projection

        # Transverse (xy) - Add histogram above heatmap, x-axis aligned
        fig = plt.figure(i, figsize=(5,5))
        gs = GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1], hspace=0, wspace=0)
        axxyhm = fig.add_subplot(gs[1, 0])   # Heatmap
        axxhist = fig.add_subplot(gs[0, 0], sharex=axxyhm) # x hist
        axyhist = fig.add_subplot(gs[1, 1], sharey=axxyhm) # y hist
        infobox = fig.add_subplot(gs[0, 1]) # infobox

        # Heatmap
        extentx = (centers[i][0] - wsize[0]*vox[0]/2, centers[i][0] + wsize[0]*vox[0]/2)
        extenty = (centers[i][1] - wsize[1]*vox[1]/2, centers[i][1] + wsize[1]*vox[1]/2)
        axxyhm.set_xlim(extentx)
        axxyhm.set_ylim(extenty)
        if vmax == 0:
            vmax = np.max(imgxy) * 0.8
        xyhm = axxyhm.imshow(
            imgxy.T,
            cmap=cmap, 
            extent=(extentx[0], extentx[1], extenty[0], extenty[1]),
            aspect='auto',
            #vmax = vmax
        )
        axxyhm.set_xticks(np.linspace(extentx[0], extentx[1] + 0.0001, 5))
        axxyhm.set_yticks(np.linspace(extenty[0], extenty[1] + 0.0001, 5))
        axxyhm.grid(True, linewidth=0.5, color = 'white')
        axxyhm.set_xlabel("x (mm)")
        axxyhm.set_ylabel("y (mm)")

        # x histogram
        x_bins = np.linspace(extentx[0] - vox[0], extentx[1] + vox[0], wsize[0]+3)
        bin_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
        axxhist.step(bin_centers, extendEdges(imgx), color='red', where='mid')
        bin_centers = bin_centers[1:-1] # Adjustment for resolution calculations
        if fit:
            # FWHM and FWTM calculations
            # Start from the peak and look left for L (low x), right for R (high x)
            peak_idx = np.argmax(imgx)
            # Find L (left) crossing for FWTM (10%)
            i = peak_idx
            while i > 0 and imgx[i] >= 0.1 * np.max(imgx):
                i -= 1
            # Interpolate crossing for FWTMxL
            x1, y1 = bin_centers[i], imgx[i]
            x2, y2 = bin_centers[i+1], imgx[i+1]
            FWTMxL = x2 + (0.1 * np.max(imgx) - y2) * (x1 - x2) / (y1 - y2)
            axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
            # Find L (left) crossing for FWHM (50%)
            i = peak_idx
            while i > 0 and imgx[i] >= 0.5 * np.max(imgx):
                i -= 1
            x1, y1 = bin_centers[i], imgx[i]
            x2, y2 = bin_centers[i+1], imgx[i+1]
            FWHMxL = x2 + (0.5 * np.max(imgx) - y2) * (x1 - x2) / (y1 - y2)
            axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)

            # Find R (right) crossing for FWTM (10%)
            i = peak_idx
            while i < len(imgx) - 1 and imgx[i] >= 0.1 * np.max(imgx):
                i += 1
            x1, y1 = bin_centers[i], imgx[i]
            x2, y2 = bin_centers[i-1], imgx[i-1]
            FWTMxR = x2 + (0.1 * np.max(imgx) - y2) * (x1 - x2) / (y1 - y2)
            axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
            # Find R (right) crossing for FWHM (50%)
            i = peak_idx
            while i < len(imgx) - 1 and imgx[i] >= 0.5 * np.max(imgx):
                i += 1
            x1, y1 = bin_centers[i], imgx[i]
            x2, y2 = bin_centers[i-1], imgx[i-1]
            FWHMxR = x2 + (0.5 * np.max(imgx) - y2) * (x1 - x2) / (y1 - y2)
            axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
            axxhist.plot([FWHMxL, FWHMxR], [0.5*np.max(imgx), 0.5*np.max(imgx)], color='black', linewidth=0.8, linestyle='--')
            axxhist.plot([FWTMxL, FWTMxR], [0.1*np.max(imgx), 0.1*np.max(imgx)], color='black', linewidth=0.8, linestyle='--')
            # Parabolic fit using 3 bins surrounding max
            center_bin = np.argmax(imgx)
            x1, y1 = bin_centers[center_bin-1], imgx[center_bin-1]
            x2, y2 = bin_centers[center_bin], imgx[center_bin]
            x3, y3 = bin_centers[center_bin+1], imgx[center_bin+1]
            a, b, c = np.polyfit([x1, x2, x3], [y1, y2, y3], 2)
            x_fit = np.linspace(x1, x3, 100)
            y_fit = a * x_fit**2 + b * x_fit + c
            pos = -b / (2 * a)
            axxhist.plot(x_fit, y_fit, color='black', linewidth=0.8)

            axxhist.set_ylim(0, np.max(imgx)*1.1)
            axxhist.legend(handles=[Patch(alpha=0, label=f'FWHM = {(FWHMxR - FWHMxL):.2f}mm\nFWTM = {(FWTMxR - FWTMxL):.2f}mm\npos = {pos:.1f}mm')], 
            fontsize='small', frameon=False, loc='upper right', bbox_to_anchor=(1, 1))
        # y histogram
        y_bins = np.linspace(extenty[0] - vox[1], extenty[1] + vox[1], wsize[1]+3)
        bin_centers = 0.5 * (y_bins[:-1] + y_bins[1:]) + vox[1] / 2
        bin_centers = bin_centers[::-1]
        axyhist.step(extendEdges(imgy), bin_centers, color='red', where='pre')
        bin_centers = bin_centers[1:-1] - vox[1] / 2 # Adjustment for resolution calculations
        if fit:
            # FWHM and FWTM calculations
            # Find L (left) crossing for FWTM (10%)
            peak_idx = np.argmax(imgy)
            i = peak_idx
            while i > 0 and imgy[i] >= 0.1 * np.max(imgy):
                i -= 1
            y1, x1 = bin_centers[i], imgy[i]
            y2, x2 = bin_centers[i+1], imgy[i+1]
            FWTMyL = y2 + (0.1 * np.max(imgy) - x2) * (y1 - y2) / (x1 - x2)
            axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)

            # Find L (left) crossing for FWHM (50%)
            i = peak_idx
            while i > 0 and imgy[i] >= 0.5 * np.max(imgy):
                i -= 1
            y1, x1 = bin_centers[i], imgy[i]
            y2, x2 = bin_centers[i+1], imgy[i+1]
            FWHMyL = y2 + (0.5 * np.max(imgy) - x2) * (y1 - y2) / (x1 - x2)
            axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)

            # Find R (right) crossing for FWTM (10%)
            i = peak_idx
            while i < len(imgy) - 1 and imgy[i] >= 0.1 * np.max(imgy):
                i += 1
            y1, x1 = bin_centers[i], imgy[i]
            y2, x2 = bin_centers[i-1], imgy[i-1]
            FWTMyR = y2 + (0.1 * np.max(imgy) - x2) * (y1 - y2) / (x1 - x2)
            axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)

            # Find R (right) crossing for FWHM (50%)
            i = peak_idx
            while i < len(imgy) - 1 and imgy[i] >= 0.5 * np.max(imgy):
                i += 1
            y1, x1 = bin_centers[i], imgy[i]
            y2, x2 = bin_centers[i-1], imgy[i-1]
            FWHMyR = y2 + (0.5 * np.max(imgy) - x2) * (y1 - y2) / (x1 - x2)
            axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
            axyhist.plot([0.5*np.max(imgy), 0.5*np.max(imgy)], [FWHMyL, FWHMyR], color='black', linewidth=0.8, linestyle='--')
            axyhist.plot([0.1*np.max(imgy), 0.1*np.max(imgy)], [FWTMyL, FWTMyR], color='black', linewidth=0.8, linestyle='--')
            # Parabolic fit using 3 bins surrounding max
            center_bin = np.argmax(imgy)
            y1, x1 = bin_centers[center_bin-1], imgy[center_bin-1]
            y2, x2 = bin_centers[center_bin], imgy[center_bin]
            y3, x3 = bin_centers[center_bin+1], imgy[center_bin+1]
            a, b, c = np.polyfit([y1, y2, y3], [x1, x2, x3], 2)
            y_fit = np.linspace(y1, y3, 100)
            x_fit = a * y_fit**2 + b * y_fit + c
            pos = -b / (2 * a)
            axyhist.plot(x_fit, y_fit, color='black', linewidth=0.8)
            axyhist.set_xlim(0, np.max(imgy)*1.1)
            axyhist.legend(handles=[Patch(alpha=0, label=f'FWHM = \n{(FWHMyL - FWHMyR):.2f}mm\nFWTM = \n{(FWTMyL - FWTMyR):.2f}mm\npos = \n{pos:.1f}mm')], 
            fontsize='small', frameon=False, loc='upper right', bbox_to_anchor=(1, 1))

            print(FWHMxL, FWHMxR, FWTMxL, FWTMxR)
            print(FWHMyL, FWHMyR, FWTMyL, FWTMyR)
        # Corner infobox
        infobox.set_axis_off()
        infobox.text(0.5, 0.7, f'hmap sum:\n{np.sum(imgxy[sumtrim:imgxy.shape[0]-sumtrim, sumtrim:imgxy.shape[1]-sumtrim]):.4f}', fontsize='small', ha='center', va='center')
        hmapsums.append(np.sum(imgxy[sumtrim:imgxy.shape[0]-sumtrim, sumtrim:imgxy.shape[1]-sumtrim]))
        # Tick handling
        plt.setp(axxhist.get_xticklabels(), visible=False)
        plt.setp(axxhist.get_xticklines(), visible=False)
        plt.setp(axxhist.get_yticklabels(), visible=False)
        plt.setp(axxhist.get_yticklines(), visible=False)
        plt.setp(axyhist.get_xticklabels(), visible=False)
        plt.setp(axyhist.get_xticklines(), visible=False)
        plt.setp(axyhist.get_yticklabels(), visible=False)
        plt.setp(axyhist.get_yticklines(), visible=False)

        fig.suptitle(f'{name}, z-region ({(pixextentz[0] - dim[2]/2)*vox[2]}, {(pixextentz[1] - dim[2]/2)*vox[2]}) mm, it {iteration}')


        if save:
            plt.savefig(f'/home/kale-chen/Documents/PET/Sensitivity Normalization/images/{name}{tag}.png')

            if not os.path.exists(out_path):
                os.mkdir(out_path)
            plt.savefig(os.path.join(out_path, f'{name}{tag}_xy_it{iteration}_spatialresolution{extraname}.pdf'), format='pdf', dpi=300)
       
            if not os.path.exists(csv_path):
                with open(csv_path, 'w') as f:
                    f.write('name,tag,iteration,mode,FWTMx,FWHMx,FWTMy,FWHMy\n')
            with open(csv_path, 'a') as f:
                f.write(f'{name},{tag},{iteration},{mode},{np.abs(FWTMxR - FWTMxL):.3f},{np.abs(FWHMxR - FWHMxL):.3f},{np.abs(FWTMyR - FWTMyL):.3f},{np.abs(FWHMyR - FWHMyL):.3f}\n')
        else:
            plt.show()
print(imgxy.shape, imgxy[sumtrim:-sumtrim, sumtrim:-sumtrim].shape)
print(hmapsums)
print(np.mean(hmapsums), np.std(hmapsums), np.mean(hmapsums) / hmapsums[0], np.std(hmapsums) / np.mean(hmapsums))
print(hmapsums / hmapsums[0])