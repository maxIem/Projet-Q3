import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve

from Bilan_de_masse import *
from Reactions import equationsATR, equationsClassique
from Reactions_Variable import *
from Variables import *

#######################################
# Fichier contenant les differentes fonctions necessaire afin de plotter
# les flux en sortie du systeme avec soit une temperature, un ratio H2O/CH4,
# ou les deux variables.
#######################################
def variablesSMR(arg):
    temperature, pression, ratio, flux = getVariableSMR()                       # Importe les variables depuis Variables.py
    SMR_tTabMin, SMR_tTabMax, SMR_tTablength, SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength = getBorneSMR()

    if arg.lower()=='t':
        temperature_Tab = np.linspace(SMR_tTabMin,
                SMR_tTabMax,SMR_tTablength)                 # Tableau contenant les temperatures de 700 a 1400 K
        CH4_T_Tab = np.zeros(len(temperature_Tab))                              # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
        H2_T_Tab = np.zeros(len(temperature_Tab))                               # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab
        return [temperature, pression, ratio, flux,
                temperature_Tab, CH4_T_Tab, H2_T_Tab]

    elif arg.lower()=='k':
        ratio_Tab = np.linspace(SMR_ratioTabMin,SMR_ratioTabMax,
                SMR_ratioTablength)                                             # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
        CH4_K_Tab = np.zeros(len(ratio_Tab))                                    # Tableau contenant les degre d'avancement SMR pour chaque ratio de ratio_Tab
        H2_K_Tab = np.zeros(len(ratio_Tab))                                     # Tableau contenant les degre d'avancement WGS pour chaque ratio de ratio_Tab
        return [temperature, pression, ratio, flux,
                ratio_Tab, CH4_K_Tab, H2_K_Tab]

    elif arg.lower()=='tk':
        temperature_Tab = np.linspace(SMR_tTabMin,
                SMR_tTabMax,SMR_tTablength)                 # Tableau contenant les temperatures de 700 a 1400 K
        ratio_Tab = np.linspace(SMR_ratioTabMin,SMR_ratioTabMax,
                SMR_ratioTablength)                                             # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
        return [temperature, pression, ratio, flux,
                temperature_Tab, ratio_Tab]
#######################################

def variablesATR(arg):
    temperature, pression, ratio, ratioO2, flux = getVariableATR()              # Importe les variables depuis Variables.py
    SMR_tTabMin, SMR_tTabMax, SMR_tTablength, SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength = getBorneSMR()

    if arg.lower()=='t':
        temperature_Tab = np.linspace(SMR_tTabMin,
                SMR_tTabMax,SMR_tTablength)                 # Tableau contenant les temperatures de 700 a 1400 K
        CH4_T_Tab = np.zeros(len(temperature_Tab))                              # Tableau contenant les degre d'avancement SMR pour chaque temperatures de temperature_Tab
        H2_T_Tab = np.zeros(len(temperature_Tab))                               # Tableau contenant les degre d'avancement WGS pour chaque temperatures de temperature_Tab
        return [temperature, pression, ratio, ratioO2, flux,
                temperature_Tab, CH4_T_Tab, H2_T_Tab]
    elif arg.lower()=='k':
        ratio_Tab = np.linspace(SMR_ratioTabMin,SMR_ratioTabMax,
                SMR_ratioTablength)                                             # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
        CH4_K_Tab = np.zeros(len(ratio_Tab))                                    # Tableau contenant les degre d'avancement SMR pour chaque ratio de ratio_Tab
        H2_K_Tab = np.zeros(len(ratio_Tab))                                     # Tableau contenant les degre d'avancement WGS pour chaque ratio de ratio_Tab
        return [temperature, pression, ratio, ratioO2, flux,
                ratio_Tab, CH4_K_Tab, H2_K_Tab]
    elif arg.lower()=='tk':
        temperature_Tab = np.linspace(SMR_tTabMin,
                SMR_tTabMax,SMR_tTablength)                 # Tableau contenant les temperatures de 700 a 1400 K
        ratio_Tab = np.linspace(SMR_ratioTabMin,SMR_ratioTabMax,
                SMR_ratioTablength)                                             # Tableau contenant les ratio de H2O/CH4 1 a 4 bar
        return [temperature, pression, ratio, ratioO2, flux,
                temperature_Tab, ratio_Tab]
#######################################

# Plot le graphe des flux en sorties, axe x = temperatures, axe y = nombres de moles
def Bilan_de_masse_T_Variable(reaction, plot):
    i=0
    if reaction.lower()==('atr'):
        temperature, pression, ratio, ratioO2, flux, temperature_Tab, CH4_T_Tab, H2_T_Tab = variablesATR('T')
        reactionFsolve = ATRDegreAvancement
        H2_Purete_Tab = np.zeros(len(temperature_Tab))
        SMR, WGS = ATRTVariable(False)
    else :
        temperature, pression, ratio, flux, temperature_Tab, CH4_T_Tab, H2_T_Tab = variablesSMR('T')
        reactionFsolve = SMRDegreAvancement
        H2_Purete_Tab = np.zeros(len(temperature_Tab))
        SMR, WGS = ClassiqueTVariable(False)
    for i in range (0,len(SMR)):
        x = reactionFsolve(np.array([SMR[i],WGS[i]]))
        CH4_T_Tab[i] = x[0]
        H2_T_Tab[i] = x[3]
        H2_Purete_Tab[i] = x[3]/np.sum(x)
    if plot:
        plt.plot(temperature_Tab,CH4_T_Tab, label='CH4')
        plt.plot(temperature_Tab,H2_T_Tab, label='H2')
        plt.plot(temperature_Tab,H2_Purete_Tab, label='purete')
        plt.xlabel('Temperature [K]')
        plt.grid(axis='both')
        plt.legend()
        plt.show()
    else:
        return [CH4_T_Tab, H2_T_Tab, H2_Purete_Tab]

# Plot le graphe des flux en sorties, axe x = ratio H2O/CH4, axe y = nombres de moles
def Bilan_de_masse_K_Variable(reaction, plot):
    i=0
    if reaction.lower()==('atr'):
        temperature, pression, ratio, ratioO2, flux, ratio_Tab, CH4_K_Tab, H2_K_Tab = variablesATR('K')
        reactionFsolve = ATRDegreAvancement
        H2_Purete_Tab = np.zeros(len(ratio_Tab))
        SMR, WGS = ATRKVariable(False)
    else :
        temperature, pression, ratio, flux, ratio_Tab, CH4_K_Tab, H2_K_Tab = variablesSMR('K')
        reactionFsolve = SMRDegreAvancement
        H2_Purete_Tab = np.zeros(len(ratio_Tab))
        SMR, WGS = ClassiqueKVariable(False)
    for i in range (0,len(SMR)):
        x = reactionFsolve(np.array([SMR[i],WGS[i]]))
        CH4_K_Tab[i] = x[0]
        H2_K_Tab[i] = x[3]
        H2_Purete_Tab[i] = x[3]/np.sum(x)
    if plot:
        plt.plot(ratio_Tab,CH4_K_Tab, label='CH4')
        plt.plot(ratio_Tab,H2_K_Tab, label='H2')
        plt.plot(ratio_Tab,H2_Purete_Tab, label='Purete')
        plt.xlabel('Ratio H2O/CH4 à %d K' % temperature)
        plt.grid(axis='both')
        plt.legend()
        plt.show()
    else:
        return [CH4_K_Tab, H2_K_Tab, H2_Purete_Tab]

# Si plot==True
# Plot le graphedes flux en sortiest, axe x = temperatures,
# axe y = Ratio H2O/CH4, axe z = nombres de moles CH4 et H2
# Sinon renvoit les axes z
def Bilan_de_masse_TK_Variable(reaction, plot):
    if reaction.lower()==('atr'):
        temperature, pression, ratio, ratioO2, flux, temperature_Tab, ratio_Tab = variablesATR('TK')
        reactionFsolve = ATRDegreAvancement
        CH_Tab = np.ones((len(temperature_Tab),len(ratio_Tab)))
        H_Tab = np.ones((len(temperature_Tab),len(ratio_Tab)))
        SMR, WGS = ATRTKVariable(False)
    else :
        temperature, pression, ratio, flux, temperature_Tab, ratio_Tab = variablesSMR('TK')
        reactionFsolve = SMRDegreAvancement
        CH_Tab = np.ones((len(temperature_Tab),len(ratio_Tab)))
        H_Tab = np.ones((len(temperature_Tab),len(ratio_Tab)))
        SMR, WGS = ClassiqueTKVariable(False)
    for i in range (0,len(SMR)):
        for j in range (0,len(SMR[0])):
            x = reactionFsolve(np.array([SMR[i][j],WGS[i][j]]))
            CH_Tab[i][j] = x[0]
            H_Tab[i][j] = x[3]
    H_Purete_Tab = 100*H_Tab/(CH_Tab+H_Tab)
    if plot:
        plt.rcParams.update({'figure.autolayout': True})
        fig = plt.figure()
        fig.suptitle('Flux de CH4 en sortie du réacteur')# en fonction de la temperature et du ratio H2O/CH4
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(temperature_Tab, ratio_Tab)
        surf = ax.plot_surface(X, Y, CH_Tab.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='mol/s de CH4')
        plt.show()

        fig = plt.figure()
        fig.suptitle('Flux de H2 en sortie du réacteur')
        #ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, H_Tab.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='mol/s de H2')
        plt.show()

        fig = plt.figure()
        fig.suptitle('Pureté du flux de H2 gazeux en sortie du réacteur')
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, H_Purete_Tab.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='Pureté du flux [%]')
        plt.show()
    else:
        return [CH_Tab,H_Tab, H_Purete_Tab]

def debitTonne(reaction, plot, debit=1000):
    if reaction.lower()=='atr':
        temperature, pression, ratio, ratioO2, flux, temperature_Tab, ratio_Tab = variablesATR('TK')
    else :
        temperature, pression, ratio, flux, temperature_Tab, ratio_Tab = variablesSMR('TK')      # Importe les variables depuis Variables.py
    molH = Bilan_de_masse_TK_Variable(reaction, False)[1]
    f = (debit/(24*60*60))/(molH*2*10**(-3))
    if plot:
        fig = plt.figure()
        fig.suptitle('Flux de CH4 necessaire afin de produire %.2f tonne(s) de H2 par jour' % float(debit/(10**3)))
        #ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(temperature_Tab, ratio_Tab)
        surf = ax.plot_surface(X, Y, f.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='mol/s de H2')
        plt.show()
        #print(function(temperature, pression, ratio, ratioO2, f))
    else :
        return f
#Bilan_de_masse_T_Variable('SMR', True)
#Bilan_de_masse_K_Variable('SMR', True)
#Bilan_de_masse_TK_Variable('SMR', True)
#debitTonne('SMR',True)

#Bilan_de_masse_T_Variable('ATR', True)
#Bilan_de_masse_K_Variable('ATR', True)
#Bilan_de_masse_TK_Variable('ATR', True)
#debitTonne('ATR',True)
