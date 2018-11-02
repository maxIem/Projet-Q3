import numpy as np

from Variables import getVariable

#######################################
# Fichier contenant la fonction neccessaire
# afin d'appliquer la reaction de combustion
# (complete) au systeme
#######################################


# Applique la reaction combustion (CH4 + 2O2 -> CO2 + 2H2O)
# recoit [CH4 + 2O2 -> CO2 + 2H2O]
# et renvoit avec x = min(CH4, O2/2)
# [CH4-x + O2-2x -> CO2+x + H20+2x]
def combustion(flux):
    CH4_flux, O2_flux, CO2_flux, H2O_flux = flux
    delta_flux = min(CH4_flux, O2_flux/2)
    return np.array([CH4_flux-delta_flux, O2_flux-2*delta_flux, CO2_flux+delta_flux, H2O_flux+2*delta_flux])
