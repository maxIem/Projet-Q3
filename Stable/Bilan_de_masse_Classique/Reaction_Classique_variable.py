import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve

from Reaction_Classique import equationsClassique
from Variables import getVariable

#######################################
# Fichier contenant les differentes fonctions necessaire afin de calculer
# et plotter les degres d'avancement du systeme avec soit une temperature,
# un ratio H2O/CH4, ou les deux variables (ainsi que la pression, mais inutile
# dans le cadre du projet).
#######################################

temperature, pression, ratio, flux = getVariable()      # Importe les variables depuis Variables.py
KSMR = 10**(-(11650/temperature) + 13.076)              # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/temperature) - 1.764)                 # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage

temperature_Tab = np.arange(700,1400,25)                # Tableau contenant les temperatures de 700 a 1400 K
SMR_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
WGS_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab
SMR_T_Tab[-1], WGS_T_Tab[-1] = [flux/2, flux/2]         # Initialise les deux premieres valeurs de fsolve pour optimiser le temps de calcul

pression_Tab = np.arange(10,40,0.2)                     # Tableau contenant les pressions de 10 a 40 bar
SMR_P_Tab = np.zeros(len(pression_Tab))                 # Tableau contenant les degre d'avancement SMR pour chaque pression de pression_Tab
WGS_P_Tab = np.zeros(len(pression_Tab))                 # Tableau contenant les degre d'avancement WGS pour chaque pression de pression_Tab
SMR_P_Tab[-1], WGS_P_Tab[-1] = [flux/2, flux/2]         # Initialise les deux premieres valeurs de fsolve pour optimiser le temps de calcul


ratio_Tab = np.arange(1,4,0.1)                          # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
SMR_K_Tab = np.zeros(len(ratio_Tab))                    # Tableau contenant les degre d'avancement SMR pour chaque ratio de ratio_Tab
WGS_K_Tab = np.zeros(len(ratio_Tab))                    # Tableau contenant les degre d'avancement WGS pour chaque ratio de ratio_Tab
SMR_K_Tab[-1], WGS_K_Tab[-1] = [flux/2, flux/2]         # Initialise les deux premieres valeurs de fsolve pour optimiser le temps de calcul

#######################################

# Plot le graphe des degres d'avancement, axe x = temperatures,
# axe y = degre avancement SMR et WGS
def ClassiqueTVariable(plot):
    i = 0
    for T in temperature_Tab:                                # Resous le systeme pour toutes les temperatures
        KSMR = 10**(-(11650/T) + 13.076)
        KWGS = 10**((1910/T) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage
        result_System = fsolve(equationsClassique,
            np.array([SMR_T_Tab[i-1],WGS_T_Tab[i-1]]), args=(KSMR,KWGS,pression,ratio,flux))
        SMR_T_Tab[i] = result_System[0]
        WGS_T_Tab[i] = result_System[1]
        i+=1
    if plot:
        plt.plot(temperature_Tab,SMR_T_Tab,label='SMR')
        plt.plot(temperature_Tab,WGS_T_Tab,label='WGS')
        plt.xlabel('Temperature [K]')
        plt.ylabel('Degré d\'avancement [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
        plt.grid(axis='both')
        plt.legend()
        plt.show()
    else:
        return [SMR_T_Tab, WGS_T_Tab]
#######################################

# Plot le graphe des degres d'avancement, axe x = Ratio H2O/CH4,
# axe y = degre avancement SMR et WGS
def ClassiqueKVariable(plot):
    i = 0
    for ratio in ratio_Tab:                                   # Resous le systeme pour toutes les temperatures
        result_System = fsolve(equationsClassique,
            np.array([SMR_K_Tab[i-1],WGS_K_Tab[i-1]]), args=(KSMR,KWGS,pression,ratio,flux))
        SMR_K_Tab[i] = result_System[0]
        WGS_K_Tab[i] = result_System[1]
        i+=1
    if plot:
        plt.plot(ratio_Tab,SMR_K_Tab,label='SMR')
        plt.plot(ratio_Tab,WGS_K_Tab,label='WGS')
        plt.xlabel('Ratio H2O/CH4 à %d K' % temperature)
        plt.ylabel('Degré d\'avancement [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
        plt.grid(axis='both')
        plt.legend()
        plt.show()
    else:
        return [SMR_K_Tab, WGS_K_Tab]
#######################################

# Plot le graphe des degres d'avancement, axe x = pression,
# axe y = degre avancement SMR et WGS
def ClassiquePVariable(plot):
    i = 0
    for pression in pression_Tab:                       # Resous le systeme pour toutes les pressions
        result_System = fsolve(equationsClassique,
            np.array([SMR_P_Tab[i-1],WGS_P_Tab[i-1]]), args=(KSMR,KWGS,pression,ratio,flux))
        SMR_P_Tab[i] = result_System[0]
        WGS_P_Tab[i] = result_System[1]
        i+=1
    if plot:
        plt.plot(pression_Tab,SMR_P_Tab,label='SMR')
        plt.plot(pression_Tab,WGS_P_Tab,label='WGS')
        plt.xlabel('Pression [bar] à %d K' % temperature)
        plt.ylabel('Degré d\'avancement [mol/s]')               # Le degre d'avancement est exprime en pourcentage du flux d'entree
        plt.grid(axis='both')
        plt.legend()
        plt.show()
    else:
        return [SMR_P_Tab, WGS_P_Tab]
#######################################

# Si plot==True
# Plot le graphe des degres d'avancement, axe x = temperatures,
# axe y = Ratio H2O/CH4, axe z = degre avancement SMR et WGS
# Sinon renvoit les axes z, SMR et WGS
def VaporeformageTKVariable(plot):
    SMR_Tab = np.ones((len(temperature_Tab),len(ratio_Tab))) * flux/2
    WGS_Tab = np.ones((len(temperature_Tab),len(ratio_Tab))) * flux/2
    i = 0
    for T in temperature_Tab:
        KSMR = 10**(-(11650/T) + 13.076)
        KWGS = 10**((1910/T) - 1.764)
        j = 0
        for ratio in ratio_Tab:                                   # Resous le systeme pour toutes les temperatures
            result_System = fsolve(equationsClassique,
                np.array([SMR_Tab[i-1][j-1],WGS_Tab[i-1][j-1]]), args=(KSMR,KWGS,pression,ratio,flux))
            if result_System[0]<0 or result_System[0]>1  or result_System[1]<0 or result_System[1]>0.5 :
                result_System = fsolve(equationsClassique, np.array([1,1]), args=(KSMR,KWGS,pression,ratio,flux))
                if result_System[0]<0 or result_System[0]>1  or result_System[1]<0 or result_System[1]>0.5 :
                    result_System = fsolve(equationsClassique, np.array([0.8,0.2]), args=(KSMR,KWGS,pression,ratio,flux))
            SMR_Tab[i][j] = result_System[0]
            WGS_Tab[i][j] = result_System[1]
            j+=1
        i+=1
    if plot==True:
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.suptitle('Degré d\'avancement du systeme Vaporeformage en fonction de la température et de la pression')
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        X, Y = np.meshgrid(temperature_Tab, ratio_Tab)
        surf = ax.plot_surface(X, Y, SMR_Tab.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Ratio H20/CH4')
        ax.set_zlabel('Degré d\'avancement Vaporeformage');

        ax = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax.plot_surface(X, Y, WGS_Tab.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Ratio H20/CH4')
        ax.set_zlabel('Degré d\'avancement WaterGasShift');
        plt.show()
    else:
        return np.array([SMR_Tab,WGS_Tab])
#######################################

# Plot des graphs des degres d'avancement en fonction de T, p et k
#######################################
#ClassiqueTVariable(True)
#ClassiqueKVariable(True)
#ClassiquePVariable(True)
#VaporeformageTKVariable(True)
#######################################
