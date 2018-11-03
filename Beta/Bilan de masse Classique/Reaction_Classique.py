import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from Variables import getVariable

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# les degres d'avancement systeme avec des parametres fixes.
#######################################

temperature, pression, ratio, flux = getVariable()      # Importe les variables depuis Variables.py
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
def equationsClassique(systemGuess, KSMR, KWGS, pression, ratio, flux):
    x, y = systemGuess
    return np.array([KSMR*(flux-x)*(ratio*flux-x-y)*((ratio+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(x-y)*(ratio*flux-x-y)-y*(3*x+y)])                      # Retourne array [equation avancement SMR, equation avancement WGS]
#######################################

# Resous le systeme pour des parametres donnes
# temperature   : valeur de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# flux          : flux de CH4 en mol/s
def Classique(temperature, p, k, flux):
    KSMR = 10**(-(11650/temperature) + 13.076)
    KWGS = 10**((1910/temperature) - 1.764)
    result_System = fsolve(equationsClassique, np.array([flux / 2, flux / 2]), args=(KSMR, KWGS, p, k, flux))
    return(result_System)
#######################################

# Retourne la valeur du systeme en focntion des degre d'avancement
# Utiliser afin de verifier si les degres d'avancement sont correctes
def verifClassique(z):
    x,y = z
    return np.array([KSMR*(flux-x)*(ratio*flux-x-y)*((ratio+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(x-y)*(ratio*flux-x-y)-y*(3*x+y)])
#######################################

#
# ATTENTION
# Laisser les '#' devant les appele de fonctions
# qui suivent afin d'executer Bilan de masse.py
# car lorsque python importe ce fichier il l'execute
#

# Verification des solutions
#######################################
#sol = Classique(1100,30.0,2.5,1.0)
#print(sol)
#print(verifClassique(sol))
#######################################
