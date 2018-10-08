from numpy import *
from sympy import *
from sympy.solvers import solve
from sympy.solvers.solveset import nonlinsolve
from matplotlib import pyplot as plt


#######################################


CpGazF = 1200                   # Capacité calorifique des gaz du four en J/(kg.K)
CpGazP = 2900                   # Capacité calorifique des gaz du procede en J/(kg.K)
dHSMR = 224000                  # Enthalpie de la réaction SMR en J/mol
dHWGS = -34000                  # Enthalpie de la réaction WGS en J/mol
dHCH4 = -803000                 # Enthalpie de la combustion du CH4 en J/mol
MmCH4 = 0.01604                 # Masse molaire du CH4 en kg/mol
MmH2O = 0.018                   # Masse molaire du H2O en kg/mol
MmO2 = 0.032                    # Masse molaire du O2 en kg/mol
MmN2 = 0.028                    # Masse molaire du N2 en kg/mol

tExt = 293.15                   # Temperature des gaz injectés directement de l'exterieur: air en Kelvin

#######################################


# tIn  : Temperature d'entree des gaz
# tOut : Temperature de sortie des gaz
# n    : Nombre de moles de CH4 en entree
# k    : Rapport molaire H2O/CH4 dans l’alimentation du reacteur
# Energie necessaire pour le SMR classique en joule
def besoinEnergieClassique(tIn, tOut, n, k):
    E1 = dHSMR * n                                          # Energie neccesaire a la reaction
    E2 = CpGazP * (tOut-tIn) * (n*MmCH4 + n*k*MmH2O)        # Energie necessaire pour elever la temperature des gaz
    return(E1 + E2)

#######################################


# e : Energie en joule
# Methane necessaire pour une certaine quantite d'energie dans un four classique en moles
def besoinMethaneEnergie(e):
    E1 = CpGazP * (1223 - tExt) * (2*(20/100)*MmO2 + 2*(78/100)*MmN2 + 1*MmCH4) #*n   # Energie necessaire pour augmenter la temperature des gaz a bruler jusque 1223K
    E2 = dHCH4
    n = -e / (E1 + E2)
    return(n)

#######################################


# tIn  : temperature d'entree des gaz
# tOut : temperature de sortie des gaz
# n    : nombre de moles de CH4 en entree par seconde
# Surplu de methane necessaire pour le SMR autotherme en moles
def besoinMethaneAutotherme(tIn, tOut, n):
    k = 1.15                                                # Rapport molaire H2O/CH4 : fixe pour Autotherme
    l = 0.6                                                 # Rapport molaire O2/CH4 : fixe pour Autotherme

    E1 = dHSMR * n                                          # Energie neccesaire a la reaction
    E2 = CpGazP * (1223 - tIn) * (2*l*MmO2 + 1*MmCH4) #*x   # Energie necessaire pour augmenter la temperature des gaz a bruler jusque 1223K
    E3 = CpGazP * (tOut - tIn) * (n*MmCH4 + n*k*MmH2O)      # Energie necessaire pour elever la temperature des gaz
    E4 = dHCH4 #*x                                          # Energie degagee par la combustion du methane

    x = (-E1-E3) / (E2+E4)
    return (x)

#######################################


temp = arange(700, 1401)
mol = arange(1, 30)
plt.figure()
plt.plot(temp, besoinMethaneAutotherme(300, temp, 10), '-r',label='Autotherme')
plt.plot(temp, besoinMethaneEnergie(besoinEnergieClassique(300, temp, 10, 2.5)), '-b',label='Classique')
plt.title('Consommation en methane pour un flux de 10 mol/s')
plt.xlabel('Temperature de sortie des gaz (K)')
plt.ylabel('Quantite de methane (mol/s)')
plt.legend(loc='upper right')
plt.grid(axis='both')
plt.show()
