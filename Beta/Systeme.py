import numpy as np
from sympy import *
from sympy.solvers import solve
from sympy.solvers.solveset import nonlinsolve
from scipy.optimize import fsolve
#from gekko import GEKKO

###############################

flux = 1.0                                              # Debit d’entree de CH4 : normalisé a 1 mol/s
p = 30.0                                                # Pression en bar
k = 2.5                                                 # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
TSMR = 1100.0                                           # Temperature en kelvin dans le SMR : variable entre 700K et 1400K
TWGS = 480.0                                            # Temperature en kelvin dans le WGS : variable entre 700K et 1400K
KSMR = 10**(-(11650/TSMR) + 13.076)                     # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/TWGS) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS) lors du vaporeformage

###############################
def equationsVaporeformage(systemGuess):                #systemGuess correspond aux premieres valeurs pour resoudre le systeme
    x, y = systemGuess
    return np.array([KSMR*(flux-x)*(k*flux-x-y)*((k+1)*flux+2*x)**2-((x-y)*(3*x+y)**3)*p**2,
    KWGS*(x-y)*(k*flux-x-y)-y*(3*x+y)])                 #retourne array [equation avancement SMR, equation avancement WGS]

sol =  fsolve(equationsVaporeformage, np.array([flux,flux]))    #resous le systeme
print(sol)
