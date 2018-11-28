import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from Variables import *
from Reactions import *
from Reactions_Variable import *
from Bilan_de_masse import *
from CP import *
#######################################
temperature, pression, ratio, fluxCH4 = getVariableSMR()
SMR_tTabMin, SMR_tTabMax, SMR_tTablength, SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength = getBorneSMR()

CpmolH2OL = 75.351 	                                                            # Capacite calorifique de H2O(l) en J/(mol.K)

CpmolH2Vap = 40660                                                              # Chaleur latente de vaporisation H2O

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

tExt = 283.15                                                                   # Temperature des gaz injectes directement de l'exterieur: air en Kelvin
tExtCH4 = 273.15 + 30.5                                                         # Temperature du flux de CH4
tExtH2O = 283.15                                                                # Temperature du flux de H2O
tExtH2OG = 373.15                                                               # Temperature du flux de H2O(g)
tInVap = 693                                                                    # Temperature d'entree des gaz dans le vaporeformage en K
tExtVap = temperature                                                           # Temperature de sortie des gaz dans le vaporeformage en K
tInWGS = 570                                                                    # Temperature d'entree des gaz dans le reacteur WGS en K
tExtWGS = 480                                                                   # Temperature de dortie des gaz dans le reacteur WGS en K
tFinal = 350

#######################################


# e : Energie en joule
# Methane necessaire pour une certaine quantite d'energie dans un four classique en moles
def besoinMethaneEnergie(e, plot=False):
    SMR_tTab = array([linspace(SMR_tTabMin,SMR_tTabMax,SMR_tTablength),]*SMR_ratioTablength).T
    CpCombustion = CpmolCO2(SMR_tTab, 1223) + 2*CpmolH2O(SMR_tTab, 1223) + CpmolCH4(SMR_tTab, 1223) + CpmolO2(SMR_tTab, tExtCH4)
    E1 = CpmolCH4(SMR_tTab, tExtCH4) + 2*((21/100)*CpmolO2(SMR_tTab, tExt) + (79/100)*CpmolN2(SMR_tTab, tExt))
                                                                                # Energie necessaire pour augmenter la temperature des gaz a bruler jusque 1223K
    E2 = dHCH4 + CpCombustion
    n = -e / (E1 + E2)

    if plot:
        fig = plt.figure()
        fig.suptitle('Nombre de moles de CH4 nécessaire \n au réacteur SMR')
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(linspace(SMR_tTabMin,SMR_tTabMax,SMR_tTablength), linspace(SMR_ratioTabMin,SMR_ratioTabMax,SMR_ratioTablength))
        surf = ax.plot_surface(X, Y, n.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='Nbre de mol de CH4 [mol]')
        plt.show()
    else:
        return(n)

#######################################


def Enthalpie(x, y, tExtVap, ratio):
    global dHSMR; global dHWGS

    Tf, Ti = [tExtVap, 948]
    dHSMR1 = (dHSMR + 3*CpmolH2(Tf, Ti) + CpmolCO(Tf, Ti) - (CpmolCH4(Tf, Ti) + CpmolH2O(Tf, Ti)) )*x
    dHWGS1 = (dHWGS +  CpmolCO2(Tf, Ti) + CpmolH2(Tf, Ti) - (CpmolCO(Tf, Ti)  + CpmolH2O(Tf, Ti)) )*y
    dHWGS2 = (dHWGS +  CpmolCO2(Tf, Ti) + CpmolH2(Tf, Ti) - (CpmolCO(Tf, Ti)  + CpmolH2O(Tf, Ti)) )*minimum(x-y, ratio*fluxCH4 -x -y)

    return [dHSMR1, dHWGS1, dHWGS2]
#######################################

def Bilan_Energie(tExtVap=tExtVap, ratio=ratio, fluxCH4=fluxCH4, retourne=True):
    deltaE = ((tExtH2OG - tExtH2O)*CpmolH2OL + CpmolH2Vap)*ratio*fluxCH4
    deltaE += CpmolCH4(tInVap, tExtCH4) * fluxCH4  + CpmolH2O(tInVap, tExtH2OG)*ratio*fluxCH4
                                                                                # Rechauffement des gazs provenant de l'exterieur

    x, y = Classique(tExtVap, pression, ratio, fluxCH4)                         # Calcul des degres d'avancement
    dHSMR, dHWGS1, dHWGS2 = Enthalpie(x, y, tExtVap, ratio)                     # Calcul des enthalpie en fct de la temperature
    deltaE += CpmolCH4(tExtVap, tInVap) * fluxCH4  + CpmolH2O(tExtVap, tInVap)*ratio*fluxCH4 + dHSMR + dHWGS1
                                                                                # Rechauffement des gazs dans le reacteur et calcul des energie des reactions

    deltaE += (fluxCH4-x)*CpmolCH4(tInWGS, tExtVap) + (ratio*fluxCH4 - x - y)*CpmolH2O(tInWGS, tExtVap)  + (x-y)*CpmolCO(tInWGS, tExtVap)
    deltaE += (3*x+y)*CpmolH2(tInWGS, tExtVap) + y*CpmolCO2(tInWGS, tExtVap)
                                                                                # Refroidissement des gazs provenant du reacteur
    X = minimum(x-y, ratio*fluxCH4 -x -y)
    deltaE += dHWGS2 + (fluxCH4-x)*CpmolCH4(tFinal, tInWGS) + (ratio*fluxCH4 -x -y -X)*CpmolH2O(tFinal, tInWGS) + (x-y-X)*CpmolCO(tFinal, tInWGS)
    deltaE += (3*x+y+X)*CpmolH2(tFinal, tInWGS) + (y+X)*CpmolCO2(tFinal, tInWGS)
                                                                                # Refroidissement des gazs dans le reacteur et calcul des energie des reactions
    deltaE += -(ratio*fluxCH4 -x -y -X)*CpmolH2Vap                              # Condensation H2O
    if retourne:
        return deltaE
    else:
        print(deltaE)
#######################################

def Bilan_Energie_Variable(reaction, plot=False):
    if reaction.lower()=='atr':
        pass
    else:
        SMR_tTab = array([linspace(SMR_tTabMin,SMR_tTabMax,SMR_tTablength),]*SMR_ratioTablength).T
        SMR_kTab = array([linspace(SMR_ratioTabMin,SMR_ratioTabMax,SMR_ratioTablength),]*SMR_tTablength)

        deltaE = ((tExtH2OG - tExtH2O)*CpmolH2OL + CpmolH2Vap)*SMR_kTab*fluxCH4
        deltaE += CpmolCH4(tInVap, tExtCH4) * fluxCH4  + CpmolH2O(tInVap, tExtH2OG)*SMR_kTab*fluxCH4
                                                                                # Rechauffement des gazs provenant de l'exterieur

        x, y = TKVariable('SMR', False)                                          # Calcul des degres d'avancement
        dHSMR, dHWGS1, dHWGS2 = Enthalpie(x, y, SMR_tTab, SMR_kTab)              # Calcul des enthalpie en fct de la temperature
        deltaE += CpmolCH4(SMR_tTab, tInVap) * fluxCH4  + CpmolH2O(SMR_tTab, tInVap)*SMR_kTab*fluxCH4 + dHSMR + dHWGS1
                                                                                # Rechauffement des gazs dans le reacteur et calcul des energie des reactions

        deltaE += (fluxCH4-x)*CpmolCH4(tInWGS, SMR_tTab) + (SMR_kTab*fluxCH4 - x - y)*CpmolH2O(tInWGS, SMR_tTab)  + (x-y)*CpmolCO(tInWGS, SMR_tTab) + (3*x+y)*CpmolH2(tInWGS, SMR_tTab) + y*CpmolCO2(tInWGS, SMR_tTab)
                                                                                # Refroidissement des gazs provenant du reacteur
        X = minimum(x-y, ratio*fluxCH4 -x -y)
        deltaE += dHWGS2 +(fluxCH4-x)*CpmolCH4(tFinal, tInWGS) + (x-y-X)*CpmolCO(tFinal, tInWGS) + (SMR_kTab*fluxCH4 -x -y -X)*CpmolH2O(tFinal, tInWGS) + (y+X)*CpmolCO2(tFinal, tInWGS) + (3*x+y+X)*CpmolH2(tFinal, tInWGS)
                                                                                # Refroidissement des gazs dans le reacteur et calcul des energie des reactions
        deltaE += -(SMR_kTab*fluxCH4 -x -y -X)*CpmolH2Vap                       # Condensation H2O
    if plot:
        fig = plt.figure()
        fig.suptitle('Énergie nécessaire au réacteur SMR')
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(linspace(SMR_tTabMin,SMR_tTabMax,SMR_tTablength), linspace(SMR_ratioTabMin,SMR_ratioTabMax,SMR_ratioTablength))
        surf = ax.plot_surface(X, Y, deltaE.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='Énergie [J]')
        plt.show()
    else:
        return deltaE
#######################################

def MBTU(E):
    return E/(1055*10**6)

def Bilan_Energie_TVariable():
    bilan = linspace(SMR_tTabMin, SMR_tTabMax, SMR_tTablength)
    i=0
    for T in bilan:
        bilan[i] = Bilan_Energie(tExtVap=T,retourne=True)
        i+=1
    print(bilan)

def Bilan_Energie_KVariable():
    bilan = linspace(SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength)
    i=0
    for K in bilan:
        bilan[i] = Bilan_Energie(ratio=K, retourne=True)
        i+=1
    print(bilan)

#Bilan_Energie(retourne=False)
#Bilan_Energie(ratio=2.5)
#Bilan_Energie_TVariable()
#Bilan_Energie_KVariable()
#Bilan_Energie_Variable('SMR',True)
#besoinMethaneEnergie(Bilan_Energie_Variable('SMR', plot=False), plot=False)


#f=SMRDebitTonne(debit=100000)
#print(MBTU(Bilan_Energie(fluxCH4=f)+f*803000)*86400)
