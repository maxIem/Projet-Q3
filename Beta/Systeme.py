import numpy
from sympy import *
from sympy.solvers import solve
from sympy.solvers.solveset import nonlinsolve
from gekko import GEKKO
###############################
#x = symbols('x')                #epsilon
#y = symbols('y')                #eta
flux = 1.0
p = 30.0
k = 2.5
TSMR = 1100.0
TWGS = 480.0
KSMR = 10**(-(11650/TSMR) + 13.076)
KWGS = 10**((1910/TWGS) - 1.764)
###############################

def equationVapo():
    """
    Resout le systeme non lineaire
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
