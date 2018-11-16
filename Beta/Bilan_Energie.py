from matplotlib import pyplot as plt
from numpy import *
from Variables import *
from Reactions import *
#######################################
pression, ratio, fluxCH4 = getVariableSMR()[1:4]

CpGazF = 1200                                                                   # Capacite calorifique des gaz du four en J/(kg.K)
CpGazP = 2900                                                                   # Capacite calorifique des gaz du procede en J/(kg.K)
CpmolCH4 = 35.69                                                                # Capacite calorifique du CH4 en J/(mol.K)
CpmolH2O = 37.47                                                                # Capacite calorifique de H2O(g) en J/(mol.K)
CpmolCO = 29.1                                                                  # Capacite calorifique du CO en J/(mol.K)
CpmolCO2 = 37.135                                                               # Capacite calorifique du CO2 en J/(mol.K)
CpmolH2 = 28.82                                                                 # Capacite calorifique du H2 en J/(mol.K)

dHSMR = 224000                                                                  # Enthalpie de la réaction SMR en J/mol @ 948K
dHWGS = -34000                                                                  # Enthalpie de la réaction WGS en J/mol @ 948K
dHCH4 = -803000                                                                 # Enthalpie de la combustion du CH4 en J/mol @ 1223K

MmCH4 = 0.01604                                                                 # Masse molaire du CH4 en kg/mol
MmH2O = 0.018                                                                   # Masse molaire du H2O en kg/mol
MmCO = 0.028                                                                    # Masse molaire du CO en kg/mol
MmH2 = 0.002                                                                    # Masse molaire du H2 en kg/mol
MmCO2 = 0.044                                                                   # Masse molaire du CO2 en kg/mol
MmO2 = 0.032                                                                    # Masse molaire du O2 en kg/mol
MmN2 = 0.028                                                                    # Masse molaire du N2 en kg/mol

tExt = 293.15                                                                   # Temperature des gaz injectes directement de l'exterieur: air en Kelvin
tExtCH4 = 273.15 + 30.5                                                         # Temperature du flux de CH4
tExtH2O = 273.15 + 100                                                          # Temperature du flux de H2O
tInVap = 693                                                                    # Temperature d'entree des gaz dans le vaporeformage en K
tExtVap = getVariableSMR()[0]                                                   # Temperature de sortie des gaz dans le vaporeformage en K
tInWGS = 570                                                                    # Temperature d'entree des gaz dans le reacteur WGS en K
tExtWGS = 480                                                                   # Temperature de dortie des gaz dans le reacteur WGS en K

#######################################


def Enthalpie(x, y, tExtVap):
    global dHSMR; global dHWGS
    #print(dHSMR, dHWGS)
    nTot = fluxCH4*(ratio+1)+2*x
    dHSMR1 = dHSMR + (tExtVap - 948) * CpGazP * MmCH4
    dHWGS1 = dHWGS + CpGazP * (tExtVap - 948) * MmCO
    # Ou alors delta Cp Cp produit - Cp reactif (avec ou sans coef stochiometrique)
    dHWGS2 = dHWGS + CpGazP * (tExtWGS - 948) * MmCO
    #print(dHSMR, dHWGS1, dHWGS2)
    return [dHSMR1, dHWGS1, dHWGS2]
#######################################

def Bilan_Energie(tExtVap=tExtVap, ratio=ratio, retourne=False):
    deltaE = CpmolCH4*(tInVap-tExtCH4) * fluxCH4  + CpmolH2O*(tInVap-tExtH2O)*ratio*fluxCH4
                                                                                # Rechauffement des gazs provenant de l'exterieur

    x, y = Classique(tExtVap, pression, ratio, fluxCH4)                         # Calcul des degres d'avancement
    dHSMR, dHWGS1, dHWGS2 = Enthalpie(x, y, tExtVap)                            # Calcul des enthalpie en fct de la temperature
    deltaE += CpmolCH4*(tExtVap-tInVap) * fluxCH4  + CpmolH2O*(tExtVap-tInVap)*ratio*fluxCH4 + dHSMR * x + dHWGS1 * y
                                                                                # Rechauffement des gazs dans le reacteur et calcul des energie des reactions

    deltaE += (tInWGS - tExtVap) * ( (fluxCH4-x)*CpmolCH4 + (ratio*fluxCH4 - x - y)*CpmolH2O  + (x-y)*CpmolCO + (3*x+y)*CpmolH2 + y*CpmolCO2 )
                                                                                # Refroidissement des gazs provenant du reacteur

    deltaE += dHWGS2*(x-y) + (tExtWGS-tInWGS) * ( (ratio*fluxCH4 - 2*y)*CpmolH2O + y*CpmolCO2 + 4*x*CpmolH2  )
                                                                                # Refroidissement des gazs dans le reacteur et calcul des energie des reactions
    if retourne:
        return deltaE
    else:
        print(deltaE)

def Bilan_Energie_TVariable():
    bilan = zeros(700)
    for T in range(700,1400):
        bilan[T-700] = Bilan_Energie(tExtVap=T,retourne=True)
    print(bilan)

def Bilan_Energie_KVariable():
    bilan = zeros(30)
    i=0
    for K in arange(1,4,0.1):
        bilan[i] = Bilan_Energie(ratio=K, retourne=True)
        i+=1
    print(bilan)

def Bilan_Energie_TKVariable():
    bilan = zeros((700,30))
    for T in range(700,1400):
        i=0
        for K in arange(1,4,0.1):
            bilan[T-700][i] = Bilan_Energie(tExtVap=T, ratio=K, retourne=True)
            i+=1
    print(bilan)

#Bilan_Energie()
#Bilan_Energie(ratio=2.5)
#Bilan_Energie_TVariable()
#Bilan_Energie_KVariable()
Bilan_Energie_TKVariable()
