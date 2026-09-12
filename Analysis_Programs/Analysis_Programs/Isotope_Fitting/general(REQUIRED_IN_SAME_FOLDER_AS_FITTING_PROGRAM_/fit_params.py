#This is supposed to be an improvement on the previously-implimented way of specifying
#the name of PES-isotopes and their respective mean lifetimes. This new implimentation will boil
#down to tossing everything in a big dictionary where the key will be the name fo the isotope,
#and the value will be the mean lifetime in untis of seconds

#ADDITIONAL ISOTOPES SHOULD BE ADDDED FOLLOWNG THIS FORMAT
from numpy import log
ln2 = log(2)

### PMMA
tau_C10 = 19.3011/ln2
tau_C11 = 20.3402*60/ln2
tau_N13 = 9.965*60/ln2
tau_O15 = 122.266/ln2

### Methylene Blue Testing
tau_Si27 = 4.117/ln2
tau_P30 = 2.5*60/ln2
tau_P29 = 4.102/ln2
tau_S30 = 1.1798/ln2
tau_S31 = 2.5534/ln2
tau_Cl34m = 31.99*60/ln2
tau_K38 = 7.651*60/ln2

### Metals
tau_Co55 = 17.53*60 * 60/ln2
tau_Co57 = 271.79*60 * 24 * 60/ln2
tau_Cu58 = 3.2 * 60 * 60/ln2
tau_Cu59 = 3.33 * 60 * 60/ln2
tau_Co58 = 70.86 * 24 * 60 * 60/ln2
tau_Co58m1 = 8.853 * 60 * 60/ln2
tau_Co58m2 = 1.25 / 1000000/ln2
tau_Ni56 = 6.075 * 24 * 60 * 60 /ln2
tau_Ni57 = 35.6 * 60 * 60 /ln2
tau_Cu60 = 23.7 *60/ln2
tau_Cu62 = 9.672*60/ln2
tau_Zn60 = 2.38*60/ln2
tau_Zn61 = 89.1/ln2
tau_Zn63 = 38.47*60/ln2
tau_Co54m = 88.8/ln2
tau_Co60 = 5.27*60 * 365 * 24 * 60/ln2
tau_Cu61 = 3.33 * 60 * 60/ln2
tau_Cu64 = 12.7 * 60 * 60/ln2
tau_Mn52m = 21.1 * 60 /ln2
tau_Mn54 = 312.2 * 24 * 60 * 60/ln2
tau_Ga68 = 67.71 * 60 /ln2
tau_Cu63 = 30.1 * 60 /ln2
tau_Cu59 = 81.5/ln2
tau_Cu58 = 3.2402/ln2
tau_Ni55 = 209 / 1000 / ln2


Isotopes_Lifetimes_Dict = {
    "C10": tau_C10,
    "C11": tau_C11,
    "N13": tau_N13,
    "O15": tau_O15,
    "Si27": tau_Si27,
    "P29": tau_P29,
    "P30": tau_P30,
    "S30": tau_S30,
    "S31": tau_S31,
    "Cl34m": tau_Cl34m,
    "K38": tau_K38,
    "Co55": tau_Co55,
    "Co57": tau_Co57,
    "Co58": tau_Co58,
    "Co58m1": tau_Co58m1,
    "Co58m2": tau_Co58m2,
    "Ni56": tau_Ni56,
    "Ni57": tau_Ni57,
    "Cu58": tau_Cu58,
    "Cu59": tau_Cu59,
    "Cu60": tau_Cu60,
    "Cu61": tau_Cu61,
    "Cu62": tau_Cu62,
    "Cu64": tau_Cu64,
    "Zn60": tau_Zn60,
    "Zn61": tau_Zn61,
    "Zn63": tau_Zn63,
    "Co54m": tau_Co54m,
    "Co60": tau_Co60,
    "Mn52m": tau_Mn52m,
    "Mn54": tau_Mn54,
    "Ga68": tau_Ga68,
    "Cu63": tau_Cu63,
    "Cu59": tau_Cu59,
    "Cu58": tau_Cu58,
    "Ni55": tau_Ni55
}