import numpy as np

from Reactions import Classique, ATR
from Reaction_WGS import waterGasShift
from Variables import getVariableSMR, getVariableATR

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# les flux en sortie du systeme avec des parametres fixes
#######################################


# Resous le systeme pour des parametres donnes en appelant
# Reaction_Vaporeformage.py/Vaporeformage(Temperature,Pression,Ratio,Flux)
# et reourne les flux en sorties du reateur sous forme d'array [ CH4, H2O, CO ,H2, CO2]
def SMRFluxSortant(Temperature,Pression,Ratio,Flux):
    sol = Classique(Temperature,Pression,Ratio,Flux)
    # flux_Sortant = [ CH4, H2O, CO ,H2, CO2]
    flux_Sortant = np.array([ Flux-sol[0] , Ratio*Flux - sol[0] - sol[1] , sol[0] - sol[1] , 3*sol[0] + sol[1], sol[1]])
    #print(flux_Sortant)
    #print('################################################### SMR')
    # waterGasShift(CO, H2O, CO2, H2)
    wgs_sortant = waterGasShift([flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3]])
    flux_Sortant[2] = wgs_sortant[0]                     # CO
    flux_Sortant[1] = wgs_sortant[1]                     # H2O
    flux_Sortant[4] = wgs_sortant[2]                     # CO2
    flux_Sortant[3] = wgs_sortant[3]                     # H2
    #print(flux_Sortant, '\n ################################################### WGS')
    flux_Sortant[1] = 0                                  # Condensation H2O
    #print(flux_Sortant, '\n ################################################### Condensation')
    flux_Sortant[4] = 0                                  # Absorption CO2
    #print(flux_Sortant, '\n ################################################### Absorption')
    #print('Purete H2 %.2f%%' % (100*flux_Sortant[3]/np.sum(flux_Sortant)))
    return flux_Sortant                                  # flux_Sortant = [ CH4, H2O, CO ,H2, CO2]
#######################################

# Recoit en arg les solutions du Systeme Vaporeformage
# et reourne les flux en sorties du reacteur sous forme d'array [ CH4, H2O, CO ,H2, CO2]
def SMRDegreAvancement(arg, ratio):
    flux = getVariableSMR()[3]
    flux_Sortant = np.array([ flux-arg[0] , ratio*flux - arg[0] - arg[1] , arg[0] - arg[1] , 3*arg[0] + arg[1], arg[1]])
    wgs_sortant = waterGasShift([flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3]])
    #print(wgs_sortant)
    flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3] = wgs_sortant                     # CO
    #print('################################################### WGS')
    flux_Sortant[1] = 0                                  # Condensation H2O
    #print(flux_Sortant)
    #print('################################################### Condensation')
    flux_Sortant[4] = 0                                  # Absorption CO2
    #print(flux_Sortant)
    #print('################################################### Absorption')
    return flux_Sortant
#######################################

def SMRDebitTonne(debit=1000):
    temperature, pression, ratio, flux = getVariableSMR()      # Importe les variables depuis Variables.py
    molH = SMRFluxSortant(temperature, pression, ratio, flux)[3]
    f = (debit/(24*60*60))/(molH*2*10**(-3))
    print(f)
    print(SMRFluxSortant(temperature, pression, ratio, f))
#######################################

# Resous le systeme pour des parametres donnes en appelant
# Reaction_Vaporeformage.py/Vaporeformage(Temperature,Pression,Ratio,Flux)
# et reourne les flux en sorties du reateur sous forme d'array [ CH4, H2O, CO ,H2, CO2]
def ATRFluxSortant(Temperature,Pression,Ratio,RatioO2,flux):
    sol = ATR(Temperature,Pression,Ratio,RatioO2,flux)
    # flux_Sortant = [ CH4, H2O, CO ,H2, CO2]
    flux_Sortant = np.array([ flux*(1-RatioO2/2)-sol[0] , flux*(1+Ratio+RatioO2) - sol[0] - sol[1] , sol[0] - sol[1] , 3*sol[0] + sol[1], RatioO2/2 + sol[1]])
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
#######################################

# Recoit en arg les solutions du Systeme Vaporeformage
# et reourne les flux en sorties du reacteur sous forme d'array [ CH4, H2O, CO ,H2, CO2]
def ATRDegreAvancement(arg, ratio):
    ratioO2, flux = getVariableATR()[3:5]      # Importe les variables depuis Variables.py
    flux_Sortant = np.array([ flux*(1-ratioO2/2)-arg[0] , flux*(1+ratio+ratioO2) - arg[0] - arg[1] , arg[0] - arg[1] , 3*arg[0] + arg[1], ratioO2/2 + arg[1]])
    wgs_sortant = waterGasShift([flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3]])
    print(wgs_sortant)
    flux_Sortant[2], flux_Sortant[1], flux_Sortant[4], flux_Sortant[3] = wgs_sortant                     # CO
    #print(flux_Sortant)
    #print('################################################### WGS')
    flux_Sortant[1] = 0                                  # Condensation H2O
    #print(flux_Sortant)
    #print('################################################### Condensation')
    flux_Sortant[4] = 0                                  # Absorption CO2
    #print(flux_Sortant)
    #print('################################################### Absorption')
    return flux_Sortant
#######################################

def ATRDebitTonne(debit=1000):
    temperature, pression, ratio, ratioO2, flux = getVariableATR()
    molH = ATRFluxSortant(temperature, pression, ratio, ratioO2, flux)[3]
    f = (debit/(24*60*60))/(molH*2*10**(-3))
    print(f)
    print(ATRFluxSortant(temperature, pression, ratio, ratioO2, f))


#print(SMRFluxSortant(getVariableSMR()[0], getVariableSMR()[1], getVariableSMR()[2], getVariableSMR()[3]))
#print(SMRDegreAvancement([0.65,0.27]))
#SMRDebitTonne()

#print(ATRFluxSortant(getVariableATR()[0], getVariableATR()[1], getVariableATR()[2], getVariableATR()[3], getVariableATR()[4]))
#print(ATRDegreAvancement([0.65,0.27]))
#ATRDebitTonne()
