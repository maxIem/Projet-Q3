import numpy as np


#######################################
# Fichier contenant la fonction neccessaire
# afin d'appliquer la fonction WGS (complete)
# au systeme
#######################################


# Applique la reaction WGS (CO + H2O -> CO2 + H2)
# recoit [CO + H2O -> CO2 + H2]
# et renvoit avec x = min(CO, H2O)
# [CO-x + H2O-x -> CO2+x + H2+x]
def waterGasShift(flux):
    CO_flux, H2O_flux, CO2_flux, H2_flux = flux
    delta_flux = min(CO_flux, H2O_flux)
    return np.array([CO_flux-delta_flux, H2O_flux-delta_flux, CO2_flux+delta_flux, H2_flux+delta_flux])
