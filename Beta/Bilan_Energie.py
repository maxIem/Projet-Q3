import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from Variables import *
from Reactions import *
from Reactions_Variable import *
#######################################
temperature, pression, ratio, fluxCH4 = getVariableSMR()
SMR_tTabMin, SMR_tTabMax, SMR_tTablength, SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength = getBorneSMR()

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
tExtVap = temperature                                                           # Temperature de sortie des gaz dans le vaporeformage en K
tInWGS = 570                                                                    # Temperature d'entree des gaz dans le reacteur WGS en K
tExtWGS = 480                                                                   # Temperature de dortie des gaz dans le reacteur WGS en K
tFinal = 350

#######################################


def Enthalpie(x, y, tExtVap, ratio):
    global dHSMR; global dHWGS

    CpSMR = 3*CpmolH2*3*x + CpmolCO*x - (CpmolCH4*x + CpmolH2O*x)
    dHSMR1 = dHSMR*x + (tExtVap - 948) * CpSMR

    CpWGS = CpmolCO2 + CpmolH2 + (CpmolCO + CpmolH2O)       # + car sort y, devient -y -> -- -> +
    dHWGS1 = (dHWGS + CpWGS * (tExtVap - 948))*y

    CpWGS = CpmolCO2 + CpmolH2 + (CpmolCO + CpmolH2O)       # + car sort coef, devient -coef -> -- -> +   coef = min()
    dHWGS2 = (dHWGS + CpWGS * (tExtWGS - 948)) * minimum(x-y, ratio*fluxCH4 -x -y)

    return [dHSMR1, dHWGS1, dHWGS2]
#######################################

def Bilan_Energie(tExtVap=tExtVap, ratio=ratio, retourne=False):
    deltaE = CpmolCH4*(tInVap-tExtCH4) * fluxCH4  + CpmolH2O*(tInVap-tExtH2O)*ratio*fluxCH4
                                                                                # Rechauffement des gazs provenant de l'exterieur

    x, y = Classique(tExtVap, pression, ratio, fluxCH4)                         # Calcul des degres d'avancement
    dHSMR, dHWGS1, dHWGS2 = Enthalpie(x, y, tExtVap, ratio)                                  # Calcul des enthalpie en fct de la temperature
    deltaE += CpmolCH4*(tExtVap-tInVap) * fluxCH4  + CpmolH2O*(tExtVap-tInVap)*ratio*fluxCH4 + dHSMR + dHWGS1
                                                                                # Rechauffement des gazs dans le reacteur et calcul des energie des reactions

    deltaE += (tInWGS - tExtVap) * ( (fluxCH4-x)*CpmolCH4 + (ratio*fluxCH4 - x - y)*CpmolH2O  + (x-y)*CpmolCO + (3*x+y)*CpmolH2 + y*CpmolCO2 )
                                                                                # Refroidissement des gazs provenant du reacteur
    X = minimum(x-y, ratio*fluxCH4 -x -y)
    deltaE += dHWGS2 + (tFinal-tInWGS) * ((ratio*fluxCH4 -x -y -X)*CpmolH2O + (y+X)*CpmolCO2 + (3*x+y+X)*CpmolH2)
                                                                                # Refroidissement des gazs dans le reacteur et calcul des energie des reactions
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
        deltaE = CpmolCH4*(tInVap-tExtCH4) * fluxCH4  + CpmolH2O*(tInVap-tExtH2O)*SMR_kTab*fluxCH4
                                                                                    # Rechauffement des gazs provenant de l'exterieur

        x, y = TKVariable('SMR', False)                                             # Calcul des degres d'avancement
        dHSMR, dHWGS1, dHWGS2 = Enthalpie(x, y, SMR_tTab, SMR_kTab)                 # Calcul des enthalpie en fct de la temperature
        deltaE += CpmolCH4*(SMR_tTab-tInVap) * fluxCH4  + CpmolH2O*(SMR_tTab-tInVap)*SMR_kTab*fluxCH4 + dHSMR + dHWGS1
                                                                                    # Rechauffement des gazs dans le reacteur et calcul des energie des reactions

        deltaE += (tInWGS - SMR_tTab) * ( (fluxCH4-x)*CpmolCH4 + (SMR_kTab*fluxCH4 - x - y)*CpmolH2O  + (x-y)*CpmolCO + (3*x+y)*CpmolH2 + y*CpmolCO2 )
                                                                                    # Refroidissement des gazs provenant du reacteur
        X = minimum(x-y, ratio*fluxCH4 -x -y)
        deltaE += dHWGS2 + (tFinal-tInWGS) * ((SMR_kTab*fluxCH4 -x -y -X)*CpmolH2O + (y+X)*CpmolCO2 + (3*x+y+X)*CpmolH2)
                                                                                    # Refroidissement des gazs dans le reacteur et calcul des energie des reactions
    if plot:
        fig = plt.figure()
        fig.suptitle('Bilan d\'energie du reacteur SMR')
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(linspace(SMR_tTabMin,SMR_tTabMax,SMR_tTablength), linspace(SMR_ratioTabMin,SMR_ratioTabMax,SMR_ratioTablength))
        surf = ax.plot_surface(X, Y, deltaE.T, cmap='viridis', edgecolor='none')
        fig.colorbar(surf, shrink=0.5, aspect=10)
        ax.set(xlabel='Temperature [K]', ylabel='Ratio H20/CH4', zlabel='Bilan d\'energie [J]')
        plt.show()
    else:
        return deltaE
#######################################


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

#Bilan_Energie()
#Bilan_Energie(ratio=2.5)
#Bilan_Energie_TVariable()
#Bilan_Energie_KVariable()
Bilan_Energie_Variable('SMR',True)
