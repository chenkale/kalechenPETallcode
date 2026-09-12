import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from tqdm import tqdm

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
save = True

nbins = 200

n_theta = nbins
n_s = nbins
s_max = 160 
pixel_size = 160 * 2 / nbins  # mm
data_dir = 'Spatial Resolution/Data2/'

theta_edges = np.linspace(0, np.pi, n_theta + 1)
s_edges = np.linspace(-s_max, s_max, n_s + 1)
mapdf = pd.read_csv('/home/kale-chen/Documents/PET/TPPT_Scanner_map_adjusted.csv', usecols=[0,1,2], header=None)
map = mapdf.to_numpy()
def lor_to_sinogram_coords(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1

    # Normal vector to the LOR
    nx = -dy
    ny = dx

    norm = np.hypot(nx, ny)
    nx /= norm
    ny /= norm

    theta = np.arctan2(ny, nx)

    # Signed distance from origin
    s = x1 * nx + y1 * ny

    # Put theta in [0, pi)
    if theta < 0:
        theta += np.pi
        s = -s
    elif theta >= np.pi:
        theta -= np.pi
        s = -s

    return theta, s

def ramp_filter_sinogram(sinogram, ds):
    n_theta, n_s = sinogram.shape

    freqs = np.fft.fftfreq(n_s, d=ds)

    # ramp filter
    ramp = np.abs(freqs)

    filtered = np.zeros_like(sinogram)

    for i in range(n_theta):
        P = np.fft.fft(sinogram[i])
        filtered[i] = np.real(np.fft.ifft(P * ramp))

    return filtered

def backproject(filtered_sino, theta_centers, s_centers, X, Y):
    image = np.zeros_like(X, dtype=np.float64)

    for i, theta in enumerate(theta_centers):
        s = X * np.cos(theta) + Y * np.sin(theta)

        projection = np.interp(
            s.ravel(),
            s_centers,
            filtered_sino[i],
            left=0,
            right=0
        ).reshape(X.shape)

        image += projection

    image *= np.pi / len(theta_centers)
    radius = np.sqrt(X**2 + Y**2)
    r_max = np.min([np.max(np.abs(X)), np.max(np.abs(Y))])
    mask = radius <= r_max
    image *= mask

    return image

def sinogram_heatmap_nice(data, fig_num, cmap='hot', name = 'Sinogram'):
    fig = plt.figure(fig_num, figsize=(5, 4.5))
    gs = GridSpec(1, 2, width_ratios=[9, 1], hspace=0, wspace=0)
    axhm = fig.add_subplot(gs[0, 0])
    axcbar = fig.add_subplot(gs[0, 1])
    hm = axhm.imshow(data, cmap=cmap, aspect='auto', extent=[-s_max, s_max, 0, 180])
    axhm.set_xticks(np.linspace(-s_max, s_max, 5))
    axhm.set_yticks(np.linspace(0, 180, 5))
    axhm.set_xlabel('Radius (mm)')
    axhm.set_ylabel('Angle (deg)')
    axhm.set_title(name)
    cbar = fig.colorbar(hm, cax=axcbar)
    cbar.set_ticks(np.linspace(np.nanmin(data), np.nanmax(data), 5))
    cbar.set_ticklabels([f"{tick:.2e}" for tick in np.linspace(np.nanmin(data), np.nanmax(data), 5)]) 
    return hm

def recon_heatmap_nice(data, fig_num, cmap='hot', name = 'FBP Recon'):
    fig = plt.figure(fig_num, figsize=(5, 4.5))
    gs = GridSpec(1, 2, width_ratios=[9, 1], hspace=0, wspace=0)
    axhm = fig.add_subplot(gs[0, 0])
    axcbar = fig.add_subplot(gs[0, 1])
    hm = axhm.imshow(data, cmap=cmap, aspect='auto', extent=[-s_max, s_max, -s_max, s_max])
    axhm.set_xticks(np.linspace(-s_max, s_max, 5))
    axhm.set_yticks(np.linspace(-s_max, s_max, 5))
    axhm.set_xlabel('x (mm)')
    axhm.set_ylabel('y (mm)')
    axhm.set_title(name)
    cbar = fig.colorbar(hm, cax=axcbar)
    cbar.set_ticks(np.linspace(np.nanmin(data), np.nanmax(data), 5))
    cbar.set_ticklabels([f"{tick:.2e}" for tick in np.linspace(np.nanmin(data), np.nanmax(data), 5)]) 
    return hm

for i, name in enumerate(names, start = 0):
    num_rows = os.path.getsize(os.path.join(data_dir, name + '.bin')) // 12
    data = np.fromfile(os.path.join(data_dir, name + '.bin'), dtype=np.int16).reshape(num_rows, 6)

    # Fill sinogram
    sinogram = np.zeros((n_theta, n_s), dtype=np.float32)
    for row in tqdm(data):
        x1, y1 = map[row[1] + 3072, 0], map[row[1] + 3072, 1]
        x2, y2 = map[row[3], 0], map[row[3], 1]

        theta, s = lor_to_sinogram_coords(x1, y1, x2, y2)

        itheta = np.searchsorted(theta_edges, theta) - 1
        is_ = np.searchsorted(s_edges, s) - 1

        if 0 <= itheta < n_theta and 0 <= is_ < n_s:
            sinogram[itheta, is_] += 1

    # Ramp filter sinogram
    filtered = ramp_filter_sinogram(sinogram, 1)

    # Backproject
    coords = np.arange(-s_max, s_max + 0.00001, pixel_size)
    X, Y = np.meshgrid(coords, coords, indexing='xy')
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    s_centers = 0.5 * (s_edges[:-1] + s_edges[1:])
    ds = s_centers[1] - s_centers[0]
    recon = backproject(filtered,theta_centers,s_centers,X,Y)
    reconraw = backproject(sinogram,theta_centers,s_centers,X,Y)


    if save:
        recon.tofile(os.path.join(data_dir, name + tag + '_fbp.bin'))
        print(recon.shape)
    else:
        sinogram_heatmap_nice(sinogram, 0, name = f'Sinogram {name}')
        recon_heatmap_nice(recon, 1, name = f'FBP Recon {name}')
        plt.show()