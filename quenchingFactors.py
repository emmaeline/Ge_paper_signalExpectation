import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import os

# this is a Ge quenching factor function that inputs an energy (in keV) and
# calculates the percent quenching, and returns the electron
# equivalent energy
def QF_StandLind(x):
    k = 0.157
    M = 72.63  # average of isotopes
    Z = 32
    eConst = np.sqrt(1439.9764)  # units of e in sqrt(keV * fm)
    a = 0.6261 * 52900 * Z**(-1/3)  # Bohr radius in fm = 52900
    ep = x * (a * M)/(Z * Z * eConst**2 * (M + M))
    g_ep = 3 * ep**(0.15) + 0.7 * ep**(0.6) + ep
    QF_val = k * g_ep / (1 + k * g_ep)
    elecEquiv = (QF_val * x)
    return elecEquiv

def QF_TUNL(x):
    k = 0.142
    M = 72.63  # average of isotopes
    Z = 32
    eConst = np.sqrt(1439.9764)  # units of e in sqrt(keV * fm)
    a = 0.6261 * 52900 * Z**(-1/3)  # Bohr radius in fm = 52900
    ep = x * (a * M)/(Z * Z * eConst**2 * (M + M))
    g_ep = 3 * ep**(0.15) + 0.7 * ep**(0.6) + ep
    QF_val = k * g_ep / (1 + k * g_ep)
    elecEquiv = (QF_val * x)
    return elecEquiv

def QF_Scholz(x):
    k = 0.1789
    M = 72.63  # average of isotopes
    Z = 32
    eConst = np.sqrt(1439.9764)  # units of e in sqrt(keV * fm)
    a = 0.6261 * 52900 * Z**(-1/3)  # Bohr radius in fm = 52900
    ep = x * (a * M)/(Z * Z * eConst**2 * (M + M))
    g_ep = 3 * ep**(0.15) + 0.7 * ep**(0.6) + ep
    QF_val = k * g_ep / (1 + k * g_ep)
    elecEquiv = (QF_val * x)
    return elecEquiv

# this outputs the quenching factor percentage of the input energy
def percentQF(x,kVal):
    k = 0.157
    M = 72.63  # average of isotopes
    Z = 32
    eConst = np.sqrt(1439.9764)  # units of e in sqrt(keV * fm)
    a = 0.6261 * 52900 * Z**(-1/3)  # Bohr radius in fm = 52900
    ep = x * (a * M)/(Z * Z * eConst**2 * (M + M))
    g_ep = 3 * ep**(0.15) + 0.7 * ep**(0.6) + ep
    QF_val = kVal * g_ep / (1 + kVal * g_ep)
    return QF_val * 100