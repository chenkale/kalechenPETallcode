import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, a, mu, c):
    return a * np.exp(-(x - mu)**2 / (2 * c**2))

def fitphotopeak(energies):
    energies = np.array(energies)
    hist, bins = np.histogram(energies, bins=100, range=(0, 40))
    bincenters = (bins[:-1] + bins[1:]) / 2
    fr = (np.argmax(hist[10:])-10 + 10, np.argmax(hist[10:])+10 + 10)
    popt, _ = curve_fit(gaussian, bincenters[fr[0]:fr[1]], hist[fr[0]:fr[1]], p0=[max(hist), bincenters[fr[0] + np.argmax(hist[fr[0]:fr[1]])], 2])
    return popt, fr

nominalcoeffs = (8.00000, 1.04676, 1.02734, 0.31909)
def energycal(q, a, b, c, d):
    return a * b ** (q ** c) + d * q - a

def calculatescalar(loc511):
    return 511 / energycal(loc511, *nominalcoeffs)

def obtainscalars(photopeaks): # Takes in list of photopeaks
    escalars = np.zeros(len(photopeaks))
    for i in range(len(photopeaks)):
        escalars[i] = calculatescalar(photopeaks[i])
    return escalars

def convertenergy(q, escalar, coeffs = nominalcoeffs):
    return energycal(q, *coeffs) * escalar