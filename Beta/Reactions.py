import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from Variables import getVariableSMR,getVariableATR

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# les degres d'avancement systeme avec des parametres fixes.
#######################################
def variablesSMR():
    temperature, pression, ratio, flux = getVariableSMR()                       # Importe les variables depuis Variables.py
    KSMR = 10**(-(11650/temperature) + 13.076)                                  # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
    KWGS = 10**((1910/temperature) - 1.764)                                     # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage
    return [temperature, pression, ratio, flux, KSMR, KWGS]
#######################################

# systemGuess   : valeur utilise pour resoudre le systeme
# KSMR          : valeur de K pour la reaction SMR, depend de la temperature
# KWGS          : valeur de K pour la reaction WGS, depend de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# flux          : flux de CH4 en mol/s
# SystemGuess correspond aux premieres valeurs pour resoudre le systeme
def equationsClassique(systemGuess, KSMR, KWGS, pression, ratio, flux):
    x, y = systemGuess
    return np.array([KSMR*(flux-x)*(ratio*flux-x-y)*((ratio+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(x-y)*(ratio*flux-x-y)-y*(3*x+y)])                                     # Retourne array [equation avancement SMR, equation avancement WGS]
#######################################

# Resous le systeme pour des parametres donnes
# temperature   : valeur de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# flux          : flux de CH4 en mol/s
def Classique(temperature, p, k, flux):
    KSMR = 10**(-(11650/temperature) + 13.076)
    KWGS = 10**((1910/temperature) - 1.764)
    result_System = fsolve(equationsClassique, np.array([flux / 2, flux / 5]), args=(KSMR, KWGS, p, k, flux))
    return(result_System)
#######################################

# Retourne la valeur du systeme en focntion des degre d'avancement
# Utiliser afin de verifier si les degres d'avancement sont correctes
def verifClassique(z):
    x,y = z
    temperature, pression, ratio, flux, KSMR, KWGS = getVariableSMR()
    KSMR = 10**(-(11650/temperature) + 13.076)
    KWGS = 10**((1910/temperature) - 1.764)
    return np.array([KSMR*(flux-x)*(ratio*flux-x-y)*((ratio+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(x-y)*(ratio*flux-x-y)-y*(3*x+y)])
#######################################

#######################################
################# ATR #################
#######################################

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# les degres d'avancement systeme avec des parametres fixes.
#######################################
def variablesATR():
    temperature, pression, ratio, ratioO2, flux = getVariableATR()              # Importe les variables depuis Variables.py
    #CH4_flux, O2_flux, CO2_flux, H2O_flux = combustion([flux, ratioO2*flux,0,0])
    KSMR = 10**(-(11650/temperature) + 13.076)                                  # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
    KWGS = 10**((1910/temperature) - 1.764)                                     # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage
    return [temperature, pression, ratio, ratioO2, flux, KSMR, KWGS]
#######################################

# systemGuess   : valeur utilise pour resoudre le systeme
# KSMR          : valeur de K pour la reaction SMR, depend de la temperature
# KWGS          : valeur de K pour la reaction WGS, depend de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# flux          : flux de CH4 en mol/s
# SystemGuess correspond aux premieres valeurs pour resoudre le systeme
def equationsATR(systemGuess,KSMR,KWGS,pression,ratio,ratioO2,flux):
    x, y = systemGuess
    return np.array([KSMR*(flux*(1-ratioO2/2)-x)*(flux*(ratio+ratioO2)-x-y)*(flux*(1+ratio+ratioO2)+2*x)**2 - ((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(flux*(ratio+ratioO2)-x-y)*(x-y) - ((ratioO2/2)+y)*(3*x+y)])
#######################################

# Resous le systeme pour des parametres donnes
# temperature   : valeur de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# kO2           : ratio O2/CH4
# flux          : flux de CH4 en mol/s
def ATR(temperature,p,k,kO2,flux):
    KSMR = 10**(-(11650/temperature) + 13.076)
    KWGS = 10**((1910/temperature) - 1.764)
    return(fsolve(equationsATR,np.array([0.5,0]), args=(KSMR,KWGS,p,k,kO2,flux)))
#######################################

# Retourne la valeur du systeme en focntion des degre d'avancement
# Utiliser afin de verifier si les degres d'avancement sont correctes
def verifATR(z):
    x,y = z
    temperature, pression, ratio, ratioO2, flux, KSMR, KWGS = getVariableATR()
    KSMR = 10**(-(11650/temperature) + 13.076)
    KWGS = 10**((1910/temperature) - 1.764)
    return np.array([KSMR*(flux*(1-ratioO2/2)-x)*(flux*(ratio+ratioO2)-x-y)*(flux*(1+ratio+ratioO2)+2*x)**2 - ((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(flux*(ratio+ratioO2)-x-y)*(x-y) - ((ratioO2/2)+y)*(3*x+y)])
#######################################
