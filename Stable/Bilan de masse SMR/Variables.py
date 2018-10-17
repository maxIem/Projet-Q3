import numpy as np




# Fichier servant a unifier les variables
# dans tous les autres
#######################################
flux = 1.0                                                     # Debit d’entree de CH4 : normalisé a 1 mol/s
pression = 30.0                                                # Pression en bar
ratio = 2.5                                                    # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
temperature = 1100.0                                           # Temperature en kelvin dans le SMR : variable entre 700K et 1400K
KSMR = 10**(-(11650/temperature) + 13.076)                     # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/temperature) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage
#######################################


# Retourne les variables neccesaire au bilan de masse du reacteur
def getVariable():
    return np.array([temperature,pression,ratio,flux])
