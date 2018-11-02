import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve

from Bilan_de_masse import (vaporeformageDegreAvancement,
                            vaporeformageFluxSortant)
from Reaction_ATR import equationsATR
from Reaction_ATR_variable import VaporeformageTKVariable
from Variables import getVariable

#######################################
# Fichier contenant les differentes fonctions necessaire afin de plotter
# les flux en sortie du systeme avec soit une temperature, un ratio H2O/CH4,
# ou les deux variables.
#######################################

temperature, pression, ratio, ratioO2, flux = getVariable()      # Importe les variables depuis Variables.py
KSMR = 10**(-(11650/temperature) + 13.076)              # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/temperature) - 1.764)                 # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage

temperature_Tab = np.arange(700,1400,25)                # Tableau contenant les temperatures de 700 a 1400 K
CH4_T_Tab = np.zeros(len(temperature_Tab))              # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
H2_T_Tab = np.zeros(len(temperature_Tab))               # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab

ratio_Tab = np.arange(1,4,0.1)                          # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
CH4_K_Tab = np.zeros(len(ratio_Tab))                    # Tableau contenant les degre d'avancement SMR pour chaque ratio de ratio_Tab
H2_K_Tab = np.zeros(len(ratio_Tab))                     # Tableau contenant les degre d'avancement WGS pour chaque ratio de ratio_Tab

#######################################

# Plot le graphe des flux en sorties, axe x = temperatures, axe y = nombres de moles
def Bilan_de_masse_T_Variable():
    i=0
    H2_Purete_Tab = np.zeros(len(temperature_Tab))
    for T in temperature_Tab:
        x = vaporeformageFluxSortant(T,pression, ratio, ratioO2, flux)
        CH4_T_Tab[i] = x[0]
        H2_T_Tab[i] = x[3]
        H2_Purete_Tab[i] = x[3]/np.sum(x)
        i+=1
    plt.plot(temperature_Tab,CH4_T_Tab, label='CH4')
    plt.plot(temperature_Tab,H2_T_Tab, label='H2')
    plt.plot(temperature_Tab,H2_Purete_Tab, label='purete')
    plt.xlabel('Temperature [K]')
    plt.grid(axis='both')
    plt.legend()
    plt.show()

# Plot le graphe des flux en sorties, axe x = ratio H2O/CH4, axe y = nombres de moles
def Bilan_de_masse_K_Variable():
    i=0
    H2_Purete_Tab = np.zeros(len(ratio_Tab))
    for k in ratio_Tab:
        x = vaporeformageFluxSortant(temperature,pression, k, ratioO2, flux)
        CH4_K_Tab[i] = x[0]
        H2_K_Tab[i] = x[3]
        H2_Purete_Tab[i] = x[3]/np.sum(x)
        i+=1
    plt.plot(ratio_Tab,CH4_K_Tab, label='CH4')
    plt.plot(ratio_Tab,H2_K_Tab, label='H2')
    plt.plot(ratio_Tab,H2_Purete_Tab, label='Purete')
    plt.xlabel('Ratio H2O/CH4 à %d K' % temperature)
    plt.grid(axis='both')
    plt.legend()
    plt.show()

# Plot le graphe des flux en sorties, axe x = ratio H2O/CH4, axe y = ratio H2O/CH4, axe z = nombres de moles
def Bilan_de_masse_TK_Variable():
    CH_Tab = np.ones((len(temperature_Tab),len(ratio_Tab)))
    H_Tab = np.ones((len(temperature_Tab),len(ratio_Tab)))
    SMR, WGS = VaporeformageTKVariable(False)
    for i in range (0,len(SMR)):
        for j in range (0,len(SMR[0])):
            x = vaporeformageDegreAvancement(np.array([SMR[i][j],WGS[i][j]]))
            CH_Tab[i][j] = x[0]
            H_Tab[i][j] = x[3]

    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.suptitle('Flux de sortie du réacteur en fonction de la temperature et du ratio H2O/CH4')
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    X, Y = np.meshgrid(temperature_Tab, ratio_Tab)
    surf = ax.plot_surface(X, Y, CH_Tab.T, cmap='viridis', edgecolor='none')
    fig.colorbar(surf, shrink=0.5, aspect=10)
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Ratio H20/CH4')
    ax.set_zlabel('mol/s de CH4');

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax.plot_surface(X, Y, H_Tab.T, cmap='viridis', edgecolor='none')
    fig.colorbar(surf, shrink=0.5, aspect=10)
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Ratio H20/CH4')
    ax.set_zlabel('mol/s de H2');
    plt.show()

    fig = plt.figure()
    fig.suptitle('Pureté du flux de H2 gazeux en sortie du réacteur ATR en fonction de la température et du ratio H2O/CH4')
    ax = fig.gca(projection='3d')
    H_Purete_Tab = 100*H_Tab/(CH_Tab+H_Tab)
    surf = ax.plot_surface(X, Y, H_Purete_Tab.T, cmap='viridis', edgecolor='none')
    fig.colorbar(surf, shrink=0.5, aspect=10)
    ax.set_xlabel('Température [K]')
    ax.set_ylabel('Ratio H20/CH4')
    ax.set_zlabel('Pureté du flux [%]');
    plt.show()


#Bilan_de_masse_T_Variable()
#Bilan_de_masse_K_Variable()
Bilan_de_masse_TK_Variable()
