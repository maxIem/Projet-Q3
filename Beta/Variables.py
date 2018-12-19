import numpy as np

#######################################
# Fichier servant a unifier les variables
# utiliser pour le systeme
#######################################

SMR_flux = 1.0                                                      # Debit d’entree de CH4 : normalisé a 1 mol/s
SMR_pression = 30.0                                                 # Pression en bar
SMR_ratio = 2.5                                                     # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
SMR_temperature = 1100.0                                            # Temperature en kelvin dans le reacteur SMR : variable entre 700K et 1400K
SMR_temperatureTabMin = 700
SMR_temperatureTabMax = 1400
SMR_temperatureTablength = 50
SMR_ratioTabMin = 1
SMR_ratioTabMax = 4
SMR_ratioTablength = 30

ATR_flux = 1.0                                                      # Debit d’entree de CH4 : normalisé a 1 mol/s
ATR_pression = 50.0                                                 # Pression en bar
ATR_ratio = 1.15                                                    # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
ATR_ratioO2 = 0.6                                                   # Rapport molaire O2/CH4 dans l’alimentation du reacteur
ATR_temperature = 1300.0                                            # Temperature en kelvin dans le reacteur ATR

#######################################


# Retourne les variables neccesaire au bilan de masse du reacteur SMR
def getVariableSMR():
    return [SMR_temperature,SMR_pression,SMR_ratio,SMR_flux]

def setVariableSMR(arg):
    global SMR_temperature
    global SMR_pression
    global SMR_ratio
    global SMR_flux
    [SMR_temperature,SMR_pression,SMR_ratio,SMR_flux] = arg

def getBorneSMR():
    return [SMR_temperatureTabMin, SMR_temperatureTabMax, SMR_temperatureTablength, SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength]

def setTBorneSMR(T1, T2):
    global SMR_temperatureTabMin
    global SMR_temperatureTabMax
    SMR_temperatureTabMin, SMR_temperatureTabMax = [T1, T2]

def setKBorneSMR(K1, K2):
    global SMR_ratioTabMin
    global SMR_ratioTabMax
    SMR_ratioTabMin, SMR_ratioTabMax = [K1, K2]

def resetBorneSMR():
    global SMR_temperatureTabMin
    global SMR_temperatureTabMax
    global SMR_ratioTabMin
    global SMR_ratioTabMax
    SMR_ratioTabMin, SMR_ratioTabMax, SMR_ratioTablength = [1, 4, 30]
    SMR_temperatureTabMin, SMR_temperatureTabMax, SMR_temperatureTablength = [700, 1400, 100]

#######################################
################# ATR #################
#######################################

# Retourne les variables neccesaire au bilan de masse du reacteur ATR
def getVariableATR():
    return np.array([ATR_temperature,ATR_pression,ATR_ratio,ATR_ratioO2,ATR_flux])

def setVariableATR(arg):
    global temperature
    global pression
    global ratio
    global ratioO2
    global flux
    [ATR_temperature,ATR_pression,ATR_ratio,ATR_ratioO2,ATR_flux] = arg

def resetVariableATR():
    global ATR_temperature
    global ATR_pression
    global ATR_ratio
    global ATR_ratioO2
    global ATR_flux
    ATR_flux = 1.0                                                      # Debit d’entree de CH4 : normalisé a 1 mol/s
    ATR_pression = 50.0                                                 # Pression en bar
    ATR_ratio = 1.15                                                    # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
    ATR_ratioO2 = 0.6                                                   # Rapport molaire O2/CH4 dans l’alimentation du reacteur
    ATR_temperature = 1300.0                                            # Temperature en kelvin dans le reacteur ATR
