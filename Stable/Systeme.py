import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


#######################################

flux = 1.0                                              # Debit d’entree de CH4 : normalisé a 1 mol/s
p = 30.0                                                # Pression en bar
k = 2.5                                                 # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
TSMR = 1100.0                                           # Temperature en kelvin dans le SMR : variable entre 700K et 1400K
TWGS = 480.0                                            # Temperature en kelvin dans le WGS : variable entre 700K et 1400K
KSMR = 10**(-(11650/TSMR) + 13.076)                     # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/TWGS) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage

temperature_Tab = np.arange(700,1401)                   # Tableau contenant les temperatures de 700 a 1400
SMR_Tab = np.zeros(len(temperature_Tab))                # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
WGS_Tab = np.zeros(len(temperature_Tab))                # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab

#######################################

# systemGuess : valeur utilise pour resoudre le systeme
def equationsVaporeformage(systemGuess):                # SystemGuess correspond aux premieres valeurs pour resoudre le systeme
    x, y = systemGuess
    return np.array([KSMR*(flux-x)*(k*flux-x-y)*((k+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*p**2,
    KWGS*(x-y)*(k*flux-x-y)-y*(3*x+y)])                 # Retourne array [equation avancement SMR, equation avancement WGS]

#######################################

def plotVaporeformage():                                # Plot le graphe des degres d'avancement, axe x = temperatures, axe y = degre avancement SMR et WGS
    plt.plot(temperature_Tab,SMR_Tab/flux*100,label='SMR')
    plt.plot(temperature_Tab,WGS_Tab/flux*100,label='WGS')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Degre d\'avancement [%]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
    plt.grid(axis='both')
    plt.legend()
    plt.show()


#######################################
#######################################

for i in range(700,1401):                               # Resous le systeme pour toutes les temperatures
    KSMR = 10**(-(11650/i) + 13.076)
    result_System = fsolve(equationsVaporeformage, np.array([SMR_Tab[i-700-1],WGS_Tab[i-700-1]]))
    SMR_Tab[i-700] = result_System[0]
    WGS_Tab[i-700] = result_System[1]

plotVaporeformage()
