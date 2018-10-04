import numpy
from sympy import *
from gekko import GEKKO

###############################

flux = 1.0                              # Debit d’entree de CH4 normalise a 1 mol/s
p = 30.0                                # Pression en bar
k = 2.5                                 # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
TSMR = 1100.0                           # Temperature en Kelvin variable entre 700 et 1400
TWGS = 480.0                            # Temperature en Kelvin variable entre 1 et 4
KSMR = 10**(-(11650/TSMR) + 13.076)
KWGS = 10**((1910/TWGS) - 1.764)

###############################

def equationVapo():
    m = GEKKO()
    x = m.Var(value=1)
    y = m.Var(value=1)
    m.Equation(KSMR*((k+1)*flux + 2*x - y)**2 * (flux-x) * (k*flux - x - y) - ((x-y) * p**2 * (3*x + y)**3)==0)
    m.Equation(KWGS*(x - y)*(k*flux - x - y) - y * (3*x+y)==0)
    m.solve(disp=False)
    print(x.value,y.value)

equationVapo()
