import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from Variables import getVariable

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
    return np.array([KSMR*(flux-x)*(ratio*flux-x-y)*((ratio+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*pression**2,
    KWGS*(x-y)*(ratio*flux-x-y)-y*(3*x+y)])                      # Retourne array [equation avancement SMR, equation avancement WGS]
#######################################

# Plot le graphe des degres d'avancement, axe x = temperatures, axe y = degre avancement SMR et WGS
def VaporeformageTvariable():
    i = 0
    for T in temperature_Tab:                                # Resous le systeme pour toutes les temperatures
        KSMR = 10**(-(11650/T) + 13.076)
        KWGS = 10**((1910/T) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage
        result_System = fsolve(equationsVaporeformage,
            np.array([SMR_T_Tab[i-1],WGS_T_Tab[i-1]]), args=(KSMR,KWGS,pression,ratio,flux))
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

# Plot le graphe des degres d'avancement, axe x = Ratio H2O/CH4, axe y = degre avancement SMR et WGS
def VaporeformageKVariable():
    i = 0
    for ratio in ratio_Tab:                                   # Resous le systeme pour toutes les temperatures
        result_System = fsolve(equationsVaporeformage,
            np.array([SMR_T_Tab[i-1],WGS_T_Tab[i-1]]), args=(KSMR,KWGS,pression,ratio,flux))
        SMR_K_Tab[i] = result_System[0]
        WGS_K_Tab[i] = result_System[1]
        i+=1
    plt.plot(ratio_Tab,SMR_K_Tab,label='SMR')
    plt.plot(ratio_Tab,WGS_K_Tab,label='WGS')
    plt.xlabel('Ratio H2O/CH4 à %d K' % temperature)
    plt.ylabel('Degré d\'avancement [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
    plt.grid(axis='both')
    plt.legend()
    plt.show()
#######################################

# Plot le graphe des degres d'avancement, axe x = pression, axe y = degre avancement SMR et WGS
def VaporeformagePvariable():
    i = 0
    for pression in pression_Tab:                       # Resous le systeme pour toutes les pressions
        result_System = fsolve(equationsVaporeformage,
            np.array([SMR_P_Tab[i-1],WGS_P_Tab[i-1]]), args=(KSMR,KWGS,pression,ratio,flux))
        SMR_P_Tab[i] = result_System[0]
        WGS_P_Tab[i] = result_System[1]
        i+=1
    plt.plot(pression_Tab,SMR_P_Tab,label='SMR')
    plt.plot(pression_Tab,WGS_P_Tab,label='WGS')
    plt.xlabel('Pression [bar] à %d K' % temperature)
    plt.ylabel('Conversion [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
    plt.grid(axis='both')
    plt.legend()
    plt.show()
#######################################

# Resous le systeme pour des parametres donnes
def Vaporeformage(temperature,p,k,flux):
    KSMR = 10**(-(11650/temperature) + 13.076)
    KWGS = 10**((1910/temperature) - 1.764)
    result_System = fsolve(equationsVaporeformage,np.array([flux/2,flux/2]), args=(KSMR,KWGS,p,k,flux))
    return(result_System)
#######################################

# Retourne la valeur du systeme en focntion des degre d'avancement
def verifVaporeformage(z):
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
# Plot des graphs des degres d'avancement en fonction de T, p et k
#######################################
#VaporeformageTvariable()
#VaporeformagePvariable()
#VaporeformageKVariable()
#######################################

# Verification des solutions
#######################################
#sol = Vaporeformage(1100,30.0,2.5,1.0)
#print(sol)
#print(verifVaporeformage(sol))
#######################################
