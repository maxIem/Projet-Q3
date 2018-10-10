import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


#######################################

flux = 1.0                                              # Debit d’entree de CH4 : normalisé a 1 mol/s
p = 30.0                                                # Pression en bar
k = 2.5                                                 # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
TSMR = 1100.0                                           # Temperature en kelvin dans le SMR : variable entre 700K et 1400K
KSMR = 10**(-(11650/TSMR) + 13.076)                     # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/TSMR) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage

temperature_Tab = np.arange(700,1401)                   # Tableau contenant les temperatures de 700 a 1400 K
SMR_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
WGS_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab

pression_Tab = np.arange(10,40,0.2)                     # Tableau contenant les pressions de 10 a 40 bar
SMR_P_Tab = np.zeros(len(pression_Tab))                 # Tableau contenant les degre d'avancement SMR pour chaque pression de pression_Tab
WGS_P_Tab = np.zeros(len(pression_Tab))                 # Tableau contenant les degre d'avancement WGS pour chaque pression de pression_Tab

#######################################


# systemGuess   : valeur utilise pour resoudre le systeme
# KSMR          : valeur de K pour la reaction SMR, depend de la temperature
# KWGS          : valeur de K pour la reaction WGS, depend de la temperature
# p             : valeur de la pression
# k             : ratio H20/CH4
# flux          : flux de CH4 en mol/s
# SystemGuess correspond aux premieres valeurs pour resoudre le systeme
def equationsVaporeformage(systemGuess,KSMR,KWGS,p,k,flux):
    x, y = systemGuess
    return np.array([KSMR*(flux-x)*(k*flux-x-y)*((k+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*p**2,
    KWGS*(x-y)*(k*flux-x-y)-y*(3*x+y)])                 # Retourne array [equation avancement SMR, equation avancement WGS]

#######################################

# Plot le graphe des degres d'avancement, axe x = temperatures, axe y = degre avancement SMR et WGS
def VaporeformageTvariable():
    i = 0
    for T in temperature_Tab:                           # Resous le systeme pour toutes les temperatures
        KSMR = 10**(-(11650/T) + 13.076)
        KWGS = 10**((1910/T) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage
        result_System = fsolve(equationsVaporeformage,
            np.array([SMR_T_Tab[i-1],WGS_T_Tab[i-1]]), args=(KSMR,KWGS,p,k,flux))
        SMR_T_Tab[i] = result_System[0]
        WGS_T_Tab[i] = result_System[1]
        i+=1
    plt.plot(temperature_Tab,SMR_T_Tab,label='SMR')
    plt.plot(temperature_Tab,WGS_T_Tab,label='WGS')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Degré d\'avancement [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
    plt.grid(axis='both')
    plt.legend()
    plt.show()

#######################################

# Resous le systeme pour une temperature donnee
def Vaporeformage(temperature):
    KSMR = 10**(-(11650/temperature) + 13.076)
    result_System = fsolve(equationsVaporeformage,np.array([flux,flux]), args=(KSMR,KWGS,p,k,flux))
    return(result_System)

#######################################


# Plot le graphe des degres d'avancement, axe x = pression, axe y = degre avancement SMR et WGS
def VaporeformagePvariable():
    i = 0
    for pression in pression_Tab:                       # Resous le systeme pour toutes les pressions
        result_System = fsolve(equationsVaporeformage,
            np.array([SMR_P_Tab[i-1],WGS_P_Tab[i-1]]), args=(KSMR,KWGS,pression,k,flux))
        SMR_P_Tab[i] = result_System[0]
        WGS_P_Tab[i] = result_System[1]
        i+=1
    plt.plot(pression_Tab,SMR_P_Tab,label='SMR')
    plt.plot(pression_Tab,WGS_P_Tab,label='WGS')
    plt.xlabel('Pression [bar] à %d K' % TSMR)
    plt.ylabel('Conversion [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
    plt.grid(axis='both')
    plt.legend()
    plt.show()

#######################################

VaporeformageTvariable()
VaporeformagePvariable()
#print(Vaporeformage(1100))
