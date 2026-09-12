
# Reconstruction viewer

from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec

names = ['0932',
         '0923',
         '0916',
         '1732',
         '1723',
         '1716',
         '2532',
         '2523',
         '2516']
tag = ''
centers = [(-78, 73, 0),
          (3, 76, 0),
          (75, 73, 0),
          (-78, 0, 0),
          (4, 0, 0),    
          (73, -2, 0),
          (-77, -78, 0),
          (4, -78, 0),
          (74, -76, 0)]


img_path = f'/home/kale-chen/Documents/PET/Spatial Resolution/Data2/'
save = False
out_path = f'/home/kale-chen/Documents/PET/Spatial Resolution/images/'
csv_path = f'/home/kale-chen/Documents/PET/Spatial Resolution/Data2/resolutions.csv'
mode = 'FBP'

auto_threshold = True # If true maxthresh is not used

dim = (201, 201, 1)
vox = (160 * 2 / 200, 160 * 2 / 200, 1)
bandrange = np.ceil(3 / vox[0]) / 2
print(bandrange)

wsize = (20, 20, 1) #mm
tick_interval = 5

cmap = 'hot'

def extendEdges(array):
    array = np.append(array, array[-1])
    array = np.insert(array, 0, array[0])
    return array

for i in range(len(names)):
    name = names[i]
    img = np.fromfile(os.path.join(img_path, f'{name}{tag}_fbp.bin'), dtype = np.float64).reshape(dim, order = 'F')
    realcenter = (centers[i][0], centers[i][1], centers[i][2])
    print(img.shape)
    #img = np.flip(img, axis=0)
    #img = np.flip(img, axis=1)
    img = np.flip(img, axis=2)

    imgxy = img.sum(axis=2)  # Image with threshold cast
    heatmapmaxindex = np.where(imgxy == np.max(imgxy))
    centerpix = (heatmapmaxindex[0][0], heatmapmaxindex[1][0])
    # Sum over a band of y pixels with width bandrange and centered on the center y

    imgx = imgxy[:, int(centerpix[1] - bandrange):int(centerpix[1] + bandrange)].sum(axis=1)  # Sum over selected band in y, gives x-projection
    imgy = imgxy[int(centerpix[0] - bandrange):int(centerpix[0] + bandrange), :].sum(axis=0)  # Sum over selected band in x, gives y-projection
    # Transverse (xy) - Add histogram above heatmap, x-axis aligned
    fig = plt.figure(i, figsize=(5,5))
    gs = GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1], hspace=0, wspace=0)
    axxyhm = fig.add_subplot(gs[1, 0])   # Heatmap
    axxhist = fig.add_subplot(gs[0, 0], sharex=axxyhm) # x hist
    axyhist = fig.add_subplot(gs[1, 1], sharey=axxyhm) # y hist
    infobox = fig.add_subplot(gs[0, 1]) # infobox

    # Heatmap
    extentx = (-dim[0] / 2 * vox[0], dim[0] / 2 * vox[0])
    extenty = (-dim[1] / 2 * vox[1], dim[1] / 2 * vox[1])
    axxyhm.set_xlim(extentx)
    axxyhm.set_ylim(extenty)
    xyhm = axxyhm.imshow(
        imgxy.T,
        cmap=cmap, 
        extent=(extentx[0], extentx[1], extenty[0], extenty[1]),
        aspect='auto'
    )
    viewextentx = (realcenter[0] - wsize[0]/2, realcenter[0] + wsize[0]/2)
    viewextenty = (realcenter[1] - wsize[1]/2, realcenter[1] + wsize[1]/2)
    axxyhm.set_xlim(viewextentx)
    axxyhm.set_ylim(viewextenty)
    axxyhm.set_xticks(np.linspace(viewextentx[0], viewextentx[1] + 0.0001, 5))
    axxyhm.set_yticks(np.linspace(viewextenty[0], viewextenty[1] + 0.0001, 5))
    axxyhm.grid(True, linewidth=0.5, color = 'white')
    axxyhm.set_xlabel("x (mm)")
    axxyhm.set_ylabel("y (mm)")

    # x histogram
    x_bins = np.linspace(extentx[0] - vox[0], extentx[1] + vox[0], dim[0]+3)
    bin_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
    axxhist.step(bin_centers, extendEdges(imgx), color='red', where='mid')
    bin_centers = bin_centers[1:-1] # Adjustment for resolution calculations
    # FWHM and FWTM calculations
    i = 0
    while imgx[i] < 0.1 * np.max(imgx):
        i += 1
    x1, y1, x2, y2 = bin_centers[i-1], imgx[i-1], bin_centers[i], imgx[i]
    axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
    FWTMxL = x1 + (0.1*np.max(imgx) - y1) * (x2 - x1) / (y2 - y1)
    while imgx[i] < 0.5 * np.max(imgx):
        i += 1
    x1, y1, x2, y2 = bin_centers[i-1], imgx[i-1], bin_centers[i], imgx[i]
    FWHMxL = x1 + (0.5*np.max(imgx) - y1) * (x2 - x1) / (y2 - y1)
    axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
    j = len(bin_centers) - 4
    while imgx[j] < 0.1 * np.max(imgx):
        j -= 1
    x1, y1, x2, y2 = bin_centers[j], imgx[j], bin_centers[j+1], imgx[j+1]
    FWTMxR = x1 + (0.1*np.max(imgx) - y1) * (x2 - x1) / (y2 - y1)
    axxhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
    while imgx[j] < 0.5 * np.max(imgx):
        j -= 1
    x1, y1, x2, y2 = bin_centers[j], imgx[j], bin_centers[j+1], imgx[j+1]
    FWHMxR = x1 + (0.5*np.max(imgx) - y1) * (x2 - x1) / (y2 - y1)
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
    y_bins = np.linspace(extenty[0] - vox[1], extenty[1] + vox[1], dim[1]+3)
    bin_centers = 0.5 * (y_bins[:-1] + y_bins[1:]) + vox[1] / 2
    bin_centers = bin_centers[::-1]
    axyhist.step(extendEdges(imgy), bin_centers, color='red', where='pre')
    bin_centers = bin_centers[1:-1] - vox[1] / 2 # Adjustment for resolution calculations
    # FWHM and FWTM calculations
    i = 0
    while imgy[i] < 0.1 * np.max(imgy):
        i += 1
    y1, y2, x1, x2 = bin_centers[i-1], bin_centers[i], imgy[i-1], imgy[i]
    FWTMyL = y1 + (0.1 * np.max(imgy) - x1) * (y2 - y1) / (x2 - x1)
    axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
    while imgy[i] < 0.5 * np.max(imgy):
        i += 1
    y1, y2, x1, x2 = bin_centers[i-1], bin_centers[i], imgy[i-1], imgy[i]
    FWHMyL = y1 + (0.5 * np.max(imgy) - x1) * (y2 - y1) / (x2 - x1)
    axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
    j = len(bin_centers) - 4
    while imgy[j] < 0.1 * np.max(imgy):
        j -= 1
    y1, y2, x1, x2 = bin_centers[j], bin_centers[j+1], imgy[j], imgy[j+1]
    FWTMyR = y1 + (0.1 * np.max(imgy) - x1) * (y2 - y1) / (x2 - x1)
    axyhist.plot([x1, x2], [y1, y2], color='black', linewidth=0.8)
    while imgy[j] < 0.5 * np.max(imgy):
        j -= 1
    y1, y2, x1, x2 = bin_centers[j], bin_centers[j+1], imgy[j], imgy[j+1]
    FWHMyR = y1 + (0.5 * np.max(imgy) - x1) * (y2 - y1) / (x2 - x1)
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
    # Corner infobox
    infobox.set_axis_off()
    infobox.text(0.5, 0.7, f'hmap sum:\n{np.sum(imgxy):.3f}', fontsize='small', ha='center', va='center')
    # Tick handling
    plt.setp(axxhist.get_xticklabels(), visible=False)
    plt.setp(axxhist.get_xticklines(), visible=False)
    plt.setp(axxhist.get_yticklabels(), visible=False)
    plt.setp(axxhist.get_yticklines(), visible=False)
    plt.setp(axyhist.get_xticklabels(), visible=False)
    plt.setp(axyhist.get_xticklines(), visible=False)
    plt.setp(axyhist.get_yticklabels(), visible=False)
    plt.setp(axyhist.get_yticklines(), visible=False)

    fig.suptitle(f'{name} FBP, binsize = {dim[0] - 1}')


    if save:
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        plt.savefig(os.path.join(out_path, f'{name}{tag}fbp_xy_spatialresolution.png'))
        if not os.path.exists(csv_path):
            with open(csv_path, 'w') as f:
                f.write('name,tag,iteration,mode,FWTMx,FWHMx,FWTMy,FWHMy\n')
        with open(csv_path, 'a') as f:
            f.write(f'{name},{tag},0,{mode},{np.abs(FWTMxR - FWTMxL):.3f},{np.abs(FWHMxR - FWHMxL):.3f},{np.abs(FWTMyR - FWTMyL):.3f},{np.abs(FWHMyR - FWHMyL):.3f}\n')
    else:
        plt.show()