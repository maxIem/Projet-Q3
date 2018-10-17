import numpy as np
from Variables import getVariable


def waterGasShift(flux):
    CO_flux, H2O_flux, CO2_flux, H2_flux = flux
    delta_flux = min(CO_flux, H2O_flux)
    return np.array([CO_flux-delta_flux, H2O_flux-delta_flux, CO2_flux+delta_flux, H2_flux+delta_flux])
