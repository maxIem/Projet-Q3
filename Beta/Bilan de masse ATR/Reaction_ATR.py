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

temperature_Tab = np.arange(700,1401)                   # Tableau contenant les temperatures de 700 a 1400 K
SMR_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
WGS_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab
SMR_T_Tab[-1] = flux/2                                  # Initialise les deux premieres valeurs de fsolve pour optimiser le temps de calcul
WGS_T_Tab[-1] = flux/2                                  # Initialise les deux premieres valeurs de fsolve pour optimiser le temps de calcul

pression_Tab = np.arange(10,40,0.2)                     # Tableau contenant les pressions de 10 a 40 bar
SMR_P_Tab = np.zeros(len(pression_Tab))                 # Tableau contenant les degre d'avancement SMR pour chaque pression de pression_Tab
WGS_P_Tab = np.zeros(len(pression_Tab))                 # Tableau contenant les degre d'avancement WGS pour chaque pression de pression_Tab

ratio_Tab = np.arange(1,4,0.1)                          # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
SMR_K_Tab = np.zeros(len(ratio_Tab))                    # Tableau contenant les degre d'avancement SMR pour chaque ratio de ratio_Tab
WGS_K_Tab = np.zeros(len(ratio_Tab))                    # Tableau contenant les degre d'avancement WGS pour chaque ratio de ratio_Tab

#######################################

# systemGuess   : valeur utilise pour resoudre le systeme
# KSMR          : valeur de K pour la reaction SMR, depend de la temperature
# KWGS          : valeur de K pour la reaction WGS, depend de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# flux          : flux de CH4 en mol/s
# SystemGuess correspond aux premieres valeurs pour resoudre le systeme
def equationsVaporeformage(systemGuess,KSMR,KWGS,pression,ratio,flux):
    x, y = systemGuess
    return np.array([,
    ])                      # Retourne array [equation avancement SMR, equation avancement WGS]
#######################################