import numpy
from sympy import *
from sympy.solvers import solve
from sympy.solvers.solveset import nonlinsolve
from sympy.solvers import solve
from scipy.optimize import fsolve
###############################
x = symbols('x')                #epsilon
y = symbols('y')                #eta
flux = 1.0
p = 30.0
k = 2.5
TSMR = 1100.0
TWGS = 480.0
KSMR = 10**((-11660/TSMR) + 13.076)
KWGS = 10**((1910/TWGS) - 1.764)
###############################
def equation():
    eqSMR = KSMR*((k+1)*flux + 2*x - y)**2 * (flux-x) * (k*flux - x) - (x-y)*p**2*27*x**3
    eqWGS = KWGS*(k*flux*(x - y) - x**2 + y**2) - y*(3*x + y)
    return nonlinsolve([eqSMR, eqWGS], [x, y])

print(equation())
