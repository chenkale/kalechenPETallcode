#Kyle Klein, 6/13/23
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ChannelPairStyles.mplstyle')


#attentuation of tungsten and LYSO and BGO at 511 keV
#interaction probsbility is 1-the exponential decay with attenuation
def Iprob(x, mu):
    e = 2.718281828
    y = 1-e**(-1*mu*x)
    return y

#x is in mm in this case
x_f = np.linspace(0,110, 1000)
#mass attenuation coefficient of tungsten * density of 19.28 g/cm^3 gives an absorbtion length of 0.26568 mm^-1

#approximate lengths for LYSO calculated using the mass attenuation length for Lu and LuO3 density of 9.02 g/cm^3
#BGO has an effectiv atomic number of 73 so the mass attenuation coeff is taken from tantalum. Density of BGO is 7.13 g/cm^3
plt.plot(x_f, Iprob(x_f, 3.048), color='b', label='Tungsten at 140 keV')
plt.plot(x_f, Iprob(x_f, 1.35), color='b', label='LYSO at 140 keV (approx.)', linestyle="--")
plt.plot(x_f, Iprob(x_f, 1.091), color = 'b', label='BGO at 140 keV (approx.)', linestyle="dotted")

#source on mass attneuation coefficient: https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z74.html
plt.plot(x_f, Iprob(x_f, 0.26568), color='r', label='Tungsten at 511 keV')
#attenaution length from LYSO taken directly from Geant4 Simulation
plt.plot(x_f, Iprob(x_f, 0.12), color='r', label='LYSO at 511 keV', linestyle="--")
plt.plot(x_f, Iprob(x_f, 0.1), color = 'r', label='BGO at 511 keV', linestyle="dotted")


plt.plot(x_f, Iprob(x_f, 0.1936), color='g', label='Tungsten at 662 keV')
plt.plot(x_f, Iprob(x_f, 0.09088), color='g', label='LYSO at 662 keV (approx.)', linestyle="--")
plt.plot(x_f, Iprob(x_f, 0.07057), color = 'g', label='BGO at 662 keV (approx.)', linestyle="dotted")
# #higher energy gammas (1-3MeV)
# plt.plot(x_f, Iprob(x_f, 0.12769), color='b', label='Tungsten at 1 MeV')
# plt.plot(x_f, Iprob(x_f, 0.08348), color='r', label='Tungsten at 2 MeV')
# plt.plot(x_f, Iprob(x_f, 0.078566), color='g', label='Tungsten at 3 MeV')
#
# plt.plot(x_f, Iprob(x_f, 0.05843), color='b', label='LYSO at 1 MeV (approx.)', linestyle="--")
# plt.plot(x_f, Iprob(x_f, 0.03955), color='r', label='LYSO at 2 MeV (approx.)', linestyle="--")
# plt.plot(x_f, Iprob(x_f, 0.03633), color='g', label='LYSO at 3 MeV (approx.)', linestyle="--")


plt.legend(fontsize=20)
plt.ylabel("Interaction Probability")
plt.xlabel("Depth [mm]")
plt.ylim(0,1)
plt.xscale('log')
plt.xlim(2*10**(-1), 100)
plt.xticks(fontsize=25)
plt.tick_params(axis='both', which='minor',width=2, length=8)
plt.tick_params(axis='both', which='major',width=3, length=15)
plt.show()
