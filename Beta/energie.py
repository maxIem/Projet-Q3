import numpy
from sympy import *
from sympy.solvers import solve
from sympy.solvers.solveset import nonlinsolve

###############################


CpGazF = 1200                   # Capacité calorifique des gaz du four en J/(kg.K)
CpGazP = 2900                   # Capacité calorifique des gaz du procede en J/(kg.K)
dHSMR = 224000                  # Enthalpie de la réaction SMR en J/mol
dHWGS = -34000                  # Enthalpie de la réaction WGS en J/mol
dHCH4 = -803000                 # Enthalpie de la combustion du CH4
MmCH4 = 0.01604                 # Masse molaire du CH4 en kg/mol
MmH2O = 0.018                   # Masse molaire du H2O en kg/mol


###############################


# tIn  : Temperature d'entree des gaz
# tOut : Temperature de sortie des gaz
# n    : Nombre de moles de CH4 en entree
# k    : Rapport molaire H2O/CH4 dans l’alimentation du reacteur

def besoinEnergieClassique(tIn, tOut, n, k):                    # Energie necessaire pour le SMR besoinEnergieClassique en joule
    E1 = dHSMR * n                                          # Energie neccesaire a la reaction
    E2 = CpGazP * (tOut-tIn) * (n*MmCH4 + n*k*MmH2O)        # Energie necessaire pour elever la temperature des gaz
    return(E1 + E2)


# tIn  : temperature d'entree des gaz
# tOut : temperature de sortie des gaz
# n    : nombre de moles de CH4 en entree
# k    : Rapport molaire H2O/CH4 dans l’alimentation du reacteur

def besoinMethaneAutotherme(tIn, tOut, n, k):                    # Methane necessaire pour le SMR autotherme en kg
    E1 = dHSMR * n                                          # Energie neccesaire a la reaction
    E2 = CpGazP * (1223 - tIn) #* x                         # Energie necessaire pour augmenter la temperature du methane a bruler jusque 1223K
    E3 = CpGazP * (tOut - tIn) * (n*MmCH4 + n*k*MmH2O)      # Energie necessaire pour elever la temperature des gaz
    E4 = dHCH4 #* x                                         # Energie degagee par la combustion du methane

    x = (-E1-E3) / (E2+E4)
    return (E1, E2, E3, E4, x)


#print(besoinEnergieClassique(693, 1100, 1, 2.5))
print(besoinMethaneAutotherme(300, 1300, 1, 2.5))
