import numpy
from sympy import *
from sympy.solvers import solve
from sympy.solvers.solveset import nonlinsolve
from gekko import GEKKO

###############################

flux = 1.0                                              # Debit d’entree de CH4 : normalisé a 1 mol/s
p = 30.0                                                # Pression en bar
k = 2.5                                                 # Rapport molaire H2O/CH4 dans l’alimentation du reacteur : variable entre 1 et 4
TSMR = 1100.0                                           # Temperature en kelvin dans le SMR : variable entre 700K et 1400K
TWGS = 480.0                                            # Temperature en kelvin dans le WGS : variable entre 700K et 1400K
KSMR = 10**(-(11650/TSMR) + 13.076)                     # Constante d’equilibre de la reaction Steam Methane Reforming (SMR)
KWGS = 10**((1910/TWGS) - 1.764)                        # Constante d’equilibre de la reaction Water–Gas Shift (WGS)

###############################


def equationVapo():
    """
    Résout le système non lineaire
    """
    m = GEKKO()
    x = m.Var(value=1)
    y = m.Var(value=1)
    m.Equation(KSMR*((k+1)*flux + 2*x - y)**2 * (flux-x) * (k*flux - x - y) - ((x-y) * p**2 * (3*x + y)**3)==0)
    m.Equation(KWGS*(x - y)*(k*flux - x - y) - y * (3*x+y)==0)
    m.solve(disp=False)
    print(x.value,y.value)

equationVapo()


"""
Ne trouve aucune solution
def equation():
    eqSMR = KSMR*((k+1)*flux + 2*x - y)**2 * (flux-x) * (k*flux - x - y) - ((x-y) * p**2 * (3*x + y)**3)
    eqWGS = KWGS*(x - y)*(k*flux - x - y) - y * (3*x+y)
    eqSMR = simplify(eqSMR)
    eqWGS = simplify(eqWGS)
    return nonlinsolve([eqSMR,eqWGS], [x, y])
print(equation())
"""
