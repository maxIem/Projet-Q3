import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from Reaction_Combustion import combustion
from Variables import getVariable

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# les degres d'avancement systeme avec des parametres fixes.
#######################################

temperature, pression, ratio, ratioO2, flux = getVariable()      # Importe les variables depuis Variables.py
#CH4_flux, O2_flux, CO2_flux, H2O_flux = combustion([flux, ratioO2*flux,0,0])
KSMR = 10**(-(11650/temperature) + 13.076)              # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/temperature) - 1.764)                 # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage

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
    return np.array([KSMR*(flux*(1-ratioO2/2)-x)*(flux*(ratio+ratioO2)-x-y)*(flux*(1+ratio+ratioO2)+2*x)**2 - ((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(flux*(ratio+ratioO2)-x-y)*(x-y) - ((ratioO2/2)+y)*(3*x+y)])
#######################################
