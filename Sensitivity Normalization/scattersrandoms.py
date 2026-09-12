import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import pandas as pd

name = 'pointsrp'
dirname = 'TPPT2026'
path = f'/home/kale-chen/Documents/PET/{dirname}/data/{name}.dat'
geomaskedoutpath = f'/home/kale-chen/Documents/PET/{dirname}/data/{name}g.bin'
intertimesmaskpath = f'/home/kale-chen/Documents/PET/{dirname}/data/{name}_intertimesmask.bin'
#npos, duration, acceptance_radius = 180, 600, 20
npos, duration, acceptance_radius = 8, 480, 20
plotsdir = f'/home/kale-chen/Documents/PET/{dirname}/writeup/plots'
geommode = 'circle1'
starttimeoffset = 0

makehistogram = True

mapdf = pd.read_csv('/home/kale-chen/Documents/PET/TPPT_Scanner_map_adjusted.csv', usecols=[0,1,2,3,4,5], header=None)
map = mapdf.to_numpy()[:, 0:3]

# Function for creating geometry acceptance mask for source position.
def get_geomask(sourcex, sourcey, within = 20):
    x1y1 = map[3072:6144, :2] - np.array([sourcex, sourcey])
    x2y2 = map[0:3072, :2]   - np.array([sourcex, sourcey])
    x1 = x1y1[:, 0][:, None]  # shape (3072, 1)
    y1 = x1y1[:, 1][:, None]
    x2 = x2y2[:, 0][None, :]  # shape (1, 3072)
    y2 = x2y2[:, 1][None, :]
    dx = x2 - x1
    dy = y2 - y1
    nx = -dy
    ny = dx
    norm = np.hypot(nx, ny)
    nx = nx / norm
    ny = ny / norm
    s = x1 * nx + y1 * ny
    table = np.abs(s) < within
    return np.asarray(table, dtype = np.bool)

boundangles = np.array([-1.12011132679, -2.02148132679, -4.26170398, -5.16307398])
boundangles = -boundangles
def sourceposlist(npos, r = 10):
    sourceposlist = []
    if geommode == 'center':
        print('center sourceposlist')
        sourceposlist.append((0, 0))
        return sourceposlist
    elif geommode == 'circle':
        print('circle sourceposlist')
        for n in range(npos):
            theta = boundangles[0] / (npos / 4) * n
            if n >= npos / 4:
                theta += boundangles[1] - boundangles[0]
            if n >= 3*npos / 4:
                theta += boundangles[1] - boundangles[0]
            sourceposlist.append((np.cos(theta) * r, -np.sin(theta) * r))
        return sourceposlist
    elif geommode == 'circle1':
        print('circle1 sourceposlist')
        for n in range(npos):
            theta = 2 * np.pi / npos * n + np.pi / 2
            sourceposlist.append((np.cos(theta) * r, np.sin(theta) * r))
        return sourceposlist
    elif geommode == 'line':
        print('line sourceposlist')
        for n in range(npos):
            sourceposlist.append((0, 120 - n * 240 / npos))
        return sourceposlist
    elif geommode == 'custom':
        print('custom sourceposlist')
        return [(0, 0), (5, 0), (10, 0), (20, 0), (40, 0), (80, 0)]

# Prepare histograms
plotrange = (-50000, 200000)
bins = np.linspace(plotrange[0], plotrange[1], 201)
#plotrange = (-100000, 5000000)
#bins = np.linspace(plotrange[0], plotrange[1], 501)
bincenters = bins[:-1] + (bins[1] - bins[0]) / 2
prehistL = np.zeros(len(bincenters), dtype = np.int32)
posthistL = np.zeros(len(bincenters), dtype = np.int32)
itimehistL = np.zeros(len(bincenters), dtype = np.int32)
prehistR = np.zeros(len(bincenters), dtype = np.int32)
posthistR = np.zeros(len(bincenters), dtype = np.int32)
itimehistR = np.zeros(len(bincenters), dtype = np.int32)


# Geomasked data for a circle run, circle position is hardcoded.
sourcepositions = sourceposlist(npos)
print(sourcepositions[0:10])
time0 = 0
if os.path.exists(geomaskedoutpath):
    os.remove(geomaskedoutpath)
if os.path.exists(intertimesmaskpath):
    os.remove(intertimesmaskpath)
totalpreevents, totalpostevents, totalitimeevents = 0, 0, 0
prevtimeL, prevtimeR = 0, 0
with open(geomaskedoutpath, 'ab') as f:
    for n in range(npos):
        timewindow = ((n * duration / npos + starttimeoffset), ((n+1) * duration / npos + starttimeoffset))
        getgeomask = get_geomask(sourcepositions[n][0], sourcepositions[n][1], within = acceptance_radius)
        #with open(f'/home/kale-chen/Documents/PET/Sensitivity Normalization/geomask{n}.bin', 'ab') as g:
            #getgeomask.tofile(g)
        neventsthispos, neventsthisposcut, neventsthisposcut2 = 0, 0, 0
        for chunk in tqdm(pd.read_csv(path, sep = '\t', chunksize = 1000000, usecols = [2, 4, 7, 9], header = None)):
            if time0 == 0:
                time0 = chunk.iloc[0, 0]
                print(time0)
            chunk.columns = ['TimeL', 'IDL', 'TimeR', 'IDR']
            chunk['TimeL'] = chunk['TimeL'] - time0
            chunk['TimeR'] = chunk['TimeR'] - time0
            timemask = ((chunk['TimeL']) / 1e12 >= timewindow[0]) & ((chunk['TimeL']) / 1e12 < timewindow[1])
            chunk = chunk[timemask]
            neventsthispos += len(chunk)
            if len(chunk) == 0:
                continue
            chunk = np.asarray(chunk)
            chunk[:, 1] = (chunk[:, 1] - 131072) % 3072
            chunk[:, 3] = chunk[:, 3] % 3072
            geommask = getgeomask[chunk[:, 1], chunk[:, 3]]
            chunk = chunk[geommask]
            towrite = np.zeros((len(chunk), 4), dtype = np.int16)
            towrite[:, 0] = np.array(chunk[:, 1], dtype = np.int16)
            towrite[:, 1] = np.array(chunk[:, 3], dtype = np.int16)
            towrite[:, 2] = np.array(chunk[:, 0] - chunk[:, 2], dtype = np.int16)
            towrite[:, 3] = np.array(chunk[:, 0] // 1e12, dtype = np.int16)
            towrite.tofile(f)
            neventsthisposcut += len(chunk)
            totalpreevents += len(chunk)
            intertimes = np.zeros((len(chunk), 2), dtype = np.int64)
            intertimes[0, 0] = chunk[0, 0] - prevtimeL
            intertimes[0, 1] = chunk[0, 2] - prevtimeR
            intertimes[1:, 0] = chunk[1:, 0] - chunk[:-1, 0]
            intertimes[1:, 1] = chunk[1:, 2] - chunk[:-1, 2]
            intertimesmask = (intertimes[:, 0] == 0) | (intertimes[:, 1] == 0)
            intertimesmask = np.asarray(intertimesmask, dtype = np.bool)
            #intertimesmask[:-1] = intertimesmask[:-1] | intertimesmask[1:]
            intertimesmask = ~intertimesmask
            neventsthisposcut2 += np.sum(intertimesmask)
            with open(intertimesmaskpath, 'ab') as g:
                intertimesmask.tofile(g)
            prevtimeL = chunk[-1, 0]
            prevtimeR = chunk[-1, 2]
            if makehistogram:
                posthistL += np.histogram(intertimes[:, 0], bins = bins)[0]
                posthistR += np.histogram(intertimes[:, 1], bins = bins)[0]
                totalpostevents += len(chunk)
                itimehistL += np.histogram(intertimes[intertimesmask, 0], bins = bins)[0]
                itimehistR += np.histogram(intertimes[intertimesmask, 1], bins = bins)[0]
                totalitimeevents += np.sum(intertimesmask)
        print(timewindow, neventsthispos, neventsthisposcut, neventsthisposcut2)
print(totalpreevents)

if makehistogram:
    totalpreevents = 0
    prevtimeL, prevtimeR = 0, 0
    for chunk in tqdm(pd.read_csv(path, sep = '\t', chunksize = 1000000, usecols = [2, 7], header = None)):
        chunk.columns = ['TimeR', 'TimeL']
        chunk = np.asarray(chunk)
        intertimes = np.zeros((len(chunk), 2), dtype = np.int64)
        intertimes[0, 0] = chunk[0, 1] - prevtimeL
        intertimes[0, 1] = chunk[0, 0] - prevtimeR
        intertimes[1:, 0] = chunk[1:, 1] - chunk[:-1, 1]
        intertimes[1:, 1] = chunk[1:, 0] - chunk[:-1, 0]
        prehistL += np.histogram(intertimes[:, 0], bins = bins)[0]
        prehistR += np.histogram(intertimes[:, 1], bins = bins)[0]
        prevtimeL = chunk[-1, 1]
        prevtimeR = chunk[-1, 0]
        totalpreevents += len(chunk)


    ylim = (1e1, 1e6)
    plt.figure(0, figsize = (6, 3))
    plt.bar(bincenters, prehistL, width = bincenters[1] - bincenters[0] + 0.001, color = 'blue', alpha = 0.4, label = 'Left arc')
    plt.bar(bincenters, prehistR, width = bincenters[1] - bincenters[0] + 0.001, color = 'red', alpha = 0.4, label = 'Right arc')
    plt.yscale('log')
    plt.ylim(1e1, 1e6)
    plt.xlim(plotrange[0], plotrange[1])
    plt.vlines([-10000, 10000], *plt.ylim(), color = 'black', linestyle = 'dotted', linewidth = 0.8, label = 'coinc window')
    plt.vlines([100000], *plt.ylim(), color = 'black', linestyle = 'dotted', linewidth = 0.8, label = '100 ns')
    plt.xlabel('Inter-event time (ps)')
    plt.legend(loc = 'lower right')
    text = f'Total events: {totalpreevents}'
    plt.text(0.97, 0.94, text, color='black', fontsize=10, ha='right', va='top',
            transform=plt.gca().transAxes)
    plt.title(f'Inter-event times, high intensity center source')
    plt.ticklabel_format(
        axis='x',
        style='sci',
        scilimits=(0, 0),
        useMathText=False
    )
    plt.savefig(os.path.join(plotsdir, f'{name}itimes.pdf'), format='pdf', dpi=300, bbox_inches='tight')
    plt.figure(1, figsize = (6, 3))
    plt.bar(bincenters, posthistL, width = bincenters[1] - bincenters[0] + 0.001, color = 'blue', alpha = 0.4, label = 'Left arc')
    plt.bar(bincenters, posthistR, width = bincenters[1] - bincenters[0] + 0.001, color = 'red', alpha = 0.4, label = 'Right arc')
    plt.yscale('log')
    plt.xlim(plotrange[0], plotrange[1])
    plt.ylim(1e1, 1e6)
    plt.vlines([-10000, 10000], *plt.ylim(), color = 'black', linestyle = 'dotted', linewidth = 0.8, label = 'coinc window')
    plt.vlines([100000], *plt.ylim(), color = 'black', linestyle = 'dotted', linewidth = 0.8, label = '100 ns')
    plt.xlabel('Inter-event time (ps)')
    plt.legend(loc = 'lower right')
    text = f'Total events: {totalpostevents}'
    plt.text(0.97, 0.94, text, color='black', fontsize=10, ha='right', va='top',
            transform=plt.gca().transAxes)
    plt.title(f'Inter-event times, high intensity center source gcut')
    plt.ticklabel_format(
        axis='x',
        style='sci',
        scilimits=(0, 0),
        useMathText=False
    )
    plt.savefig(os.path.join(plotsdir, f'{name}itimesg.pdf'), format='pdf', dpi=300, bbox_inches='tight')
    
    plt.figure(2, figsize = (6, 3))
    meantrueL = np.mean(itimehistL[bincenters > 105000])
    meandeadL = np.mean(itimehistL[(bincenters < 95000) & (bincenters > 15000)])
    meantrueR = np.mean(itimehistR[bincenters > 105000])
    meandeadR = np.mean(itimehistR[(bincenters < 95000) & (bincenters > 15000)])
    print(meantrueL, meandeadL, meantrueR, meandeadR)
    masktrue = bincenters > 100000
    maskdead = (bincenters < 100000) & (bincenters > 10000)
    print(np.sum(masktrue), np.sum(maskdead))
    correctionL = (meantrueL - meandeadL) * np.sum(maskdead)
    correctionR = (meantrueR - meandeadR) * np.sum(maskdead)
    print(f'{correctionL:.2f}, {correctionR:.2f}', f'{(correctionL + correctionR) / 2:.2f}', f'{(correctionL + correctionR) / (2*totalitimeevents):.6f}')
    plt.bar(bincenters, itimehistL, width = bincenters[1] - bincenters[0] + 0.001, color = 'blue', alpha = 0.4, label = 'Left arc')
    plt.bar(bincenters, itimehistR, width = bincenters[1] - bincenters[0] + 0.001, color = 'red', alpha = 0.4, label = 'Right arc')
    plt.yscale('log')
    plt.ylim(1e1, 1e6)
    plt.xlim(plotrange[0], plotrange[1])
    plt.vlines([-10000, 10000], *plt.ylim(), color = 'black', linestyle = 'dotted', linewidth = 0.8, label = 'coinc window')
    plt.vlines([100000], *plt.ylim(), color = 'black', linestyle = 'dotted', linewidth = 0.8, label = '100 ns')
    plt.hlines((meantrueL + meantrueR) / 2, 10000, 200000, color = 'black', linestyle = 'dotted', linewidth = 0.8)
    plt.hlines((meandeadL + meandeadR) / 2, 10000, 100000, color = 'black', linestyle = 'dotted', linewidth = 0.8)
    plt.xlabel('Inter-event time (ps)')
    plt.legend(loc = 'lower right')
    deadtimecorrection = (totalitimeevents + (correctionL + correctionR) / 2) / totalitimeevents
    text = f'Total events: {totalitimeevents}\nMissing events: {(correctionL + correctionR) / 2:.0f}\nDeadtime correction: {deadtimecorrection:.5f}'
    plt.text(0.97, 0.94, text, color='black', fontsize=10, ha='right', va='top',
            transform=plt.gca().transAxes)
    plt.title(f'Inter-event times, high intensity center source, dcut')
    plt.ticklabel_format(
        axis='x',
        style='sci',
        scilimits=(0, 0),
        useMathText=False
    )
    plt.savefig(os.path.join(plotsdir, f'{name}itimesgd.pdf'), format='pdf', dpi=300, bbox_inches='tight')
    plt.close('all')
    np.column_stack((prehistL, bincenters)).astype(np.float32).tofile(
        os.path.join(plotsdir, f'{name}prehistL.bin')
    )
    np.column_stack((prehistR, bincenters)).astype(np.float32).tofile(
        os.path.join(plotsdir, f'{name}prehistR.bin')
    )
    np.column_stack((posthistL, bincenters)).astype(np.float32).tofile(
        os.path.join(plotsdir, f'{name}posthistL.bin')
    )
    np.column_stack((posthistR, bincenters)).astype(np.float32).tofile(
        os.path.join(plotsdir, f'{name}posthistR.bin')
    )
    np.column_stack((itimehistL, bincenters)).astype(np.float32).tofile(
        os.path.join(plotsdir, f'{name}itimehistL.bin')
    )
    np.column_stack((itimehistR, bincenters)).astype(np.float32).tofile(
        os.path.join(plotsdir, f'{name}itimehistR.bin')
    )