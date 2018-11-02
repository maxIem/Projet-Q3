import numpy as np

from Reaction_ATR import ATR
from Reaction_WGS import waterGasShift
from Variables import getVariable

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# les flux en sortie du systeme avec des parametres fixes
#######################################

temperature, pression, ratio, ratioO2, flux = getVariable()      # Importe les variables depuis Variables.py
product_Flux_Tab = []                                   # Initialise le tableau des flux en sorite de forme [CH4, H2O, CO, H2, CO2]
#######################################

# Resous le systeme pour des parametres donnes en appelant
# Reaction_Vaporeformage.py/Vaporeformage(Temperature,Pression,Ratio,Flux)
# et reourne les flux en sorties du reateur sous forme d'array [ CH4, H2O, CO ,H2, CO2]
def vaporeformageFluxSortant(Temperature,Pression,Ratio,RatioO2,flux):
    sol = ATR(Temperature,Pression,Ratio,RatioO2,flux)
    # flux_Sortant = [ CH4, H2O, CO ,H2, CO2]
    flux_Sortant = np.array([ flux*(1-ratioO2/2)-sol[0] , flux*(1+ratio+ratioO2) - sol[0] - sol[1] , sol[0] - sol[1] , 3*sol[0] + sol[1], ratioO2/2 + sol[1]])
    #print(flux_Sortant)
    #print('################################################### SMR')
    # waterGasShift(CO, H2O, CO2, H2)
    wgs_sortant = waterGasShift([flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3]])
    flux_Sortant[2] = wgs_sortant[0]                     # CO
    flux_Sortant[1] = wgs_sortant[1]                     # H2O
    flux_Sortant[4] = wgs_sortant[2]                     # CO2
    flux_Sortant[3] = wgs_sortant[3]                     # H2
    #print(flux_Sortant)
    #print('################################################### WGS')
    flux_Sortant[1] = 0                                  # Condensation H2O
    #print(flux_Sortant)
    #print('################################################### Condensation')
    flux_Sortant[4] = 0                                  # Absorption CO2
    #print(flux_Sortant)
    #print('################################################### Absorption')
    #print('Purete H2 %.2f%%' % (100*flux_Sortant[3]/np.sum(flux_Sortant)))
    return flux_Sortant

# Recoit en arg les solutions du Systeme Vaporeformage
# et reourne les flux en sorties du reacteur sous forme d'array [ CH4, H2O, CO ,H2, CO2]
def vaporeformageDegreAvancement(arg):
    flux_Sortant = np.array([ flux*(1-ratioO2/2)-arg[0] , flux*(1+ratio+ratioO2) - arg[0] - arg[1] , arg[0] - arg[1] , 3*arg[0] + arg[1], ratioO2/2 + arg[1]])
    wgs_sortant = waterGasShift([flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3]])
    flux_Sortant[2] = wgs_sortant[0]                     # CO
    flux_Sortant[1] = wgs_sortant[1]                     # H2O
    flux_Sortant[4] = wgs_sortant[2]                     # CO2
    flux_Sortant[3] = wgs_sortant[3]                     # H2
    #print(flux_Sortant)
    #print('################################################### WGS')
    flux_Sortant[1] = 0                                  # Condensation H2O
    #print(flux_Sortant)
    #print('################################################### Condensation')
    flux_Sortant[4] = 0                                  # Absorption CO2
    #print(flux_Sortant)
    #print('################################################### Absorption')
    return flux_Sortant

#vaporeformageFluxSortant(temperature, pression, ratio, flux)
