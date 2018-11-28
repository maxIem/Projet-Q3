import numpy as np
from matplotlib import pyplot as plt

#######################################

x = np.linspace(273.15,1400,1000)

X = np.array([298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400])

# Calcule disponible dans les notebooks jupyter
# Valeurs proviennent de table JANAF
CpmolCH4 = np.array([35.639, 35.708, 40.5, 46.342, 52.227, 57.794, 62.932, 67.601, 71.795, 75.529, 78.833, 81.744, 84.305])
CpmolH2O = np.array([33.59, 33.596, 34.262, 35.226, 36.325, 37.495, 38.721, 39.987, 41.268, 42.536, 43.768, 44.945, 46.054])
CpmolCO2 = np.array([37.129, 37.221, 41.325, 44.627, 47.321, 49.564, 51.434, 52.999, 54.308, 55.409, 56.342, 57.137, 57.802])
CpmolCO = np.array([29.142, 29.142, 29.342, 29.974, 30.443, 31.171, 31.899, 32.577, 33.183, 33.710, 34.175, 34.572, 34.920])
CpmolH2 = np.array([28.836, 28.849, 29.181, 29.26, 29.327, 29.441, 29.624, 29.881, 30.205, 30.581, 30.992, 31.423, 31.861])
CpmolO2 = np.array([29.376, 29.385, 30.106, 31.091, 32.090, 32.981, 33.733, 34.355, 34.870, 35.3, 35.667, 35.988, 36.277])
CpmolN2 = np.array([29.124, 29.125, 29.249, 29.58, 30.11, 30.754, 31.433, 32.090, 32.697, 33.241, 33.723, 34.147, 34.518])

yh_CH4 = np.polyfit(X,CpmolCH4,3)
yh_H2O = np.polyfit(X,CpmolH2O,3)
yh_CO2 = np.polyfit(X,CpmolCO2,3)
yh_CO = np.polyfit(X,CpmolCO,3)
yh_H2 = np.polyfit(X,CpmolH2,2)
yh_O2 = np.polyfit(X,CpmolO2,3)
yh_N2 = np.polyfit(X,CpmolN2,3)

yh_CH4 = np.polyint(yh_CH4)
yh_H2O = np.polyint(yh_H2O)
yh_CO2 = np.polyint(yh_CO2)
yh_CO = np.polyint(yh_CO)
yh_H2 = np.polyint(yh_H2)
yh_O2 = np.polyint(yh_O2)
yh_N2 = np.polyint(yh_N2)

#######################################


def CpmolCH4(Tf, Ti):
    return np.polyval(yh_CH4,Tf)-np.polyval(yh_CH4,Ti)

def CpmolH2O(Tf, Ti):
    return np.polyval(yh_H2O,Tf)-np.polyval(yh_H2O,Ti)

def CpmolCO2(Tf, Ti):
    return np.polyval(yh_CO2,Tf)-np.polyval(yh_CO2,Ti)

def CpmolCO(Tf, Ti):
    return np.polyval(yh_CO,Tf)-np.polyval(yh_CO,Ti)

def CpmolH2(Tf, Ti):
    return np.polyval(yh_H2,Tf)-np.polyval(yh_H2,Ti)

def CpmolO2(Tf, Ti):
    return np.polyval(yh_O2,Tf)-np.polyval(yh_O2,Ti)

def CpmolN2(Tf, Ti):
    return np.polyval(yh_N2,Tf)-np.polyval(yh_N2,Ti)
