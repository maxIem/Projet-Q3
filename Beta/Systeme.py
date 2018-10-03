import numpy
from sympy import *
from sympy.solvers import solve

###############################
epsilon = Symbol('epsilon',real=True)
eta = Symbol('eta',real=True)
x = 1.0
p = 30.0
k = 2.5
TK1 = 1100.0
TK2 = 480.0
###############################

def equationSMR():
    K1 = 10**((-11660/TK1) + 13.076)
    K2 = 10**((1910/TK2) + 1.764)
    eq1 = K1*((k+1)*x + 2*epsilon)**2 * (x-epsilon) * (k*x - epsilon) - p**2 *27*epsilon**4
    eq2 = K2*(epsilon - eta)*(k*x - epsilon - eta) - eta*(3*epsilon + eta)

    return solve([eq1, eq2], [epsilon,eta])

###############################################################################
#x = Symbol('x')
#k = Symbol('k')
#eta = Symbol('eta')
#epsilon = Symbol('epsilon')
#print(simplify(((k+1)*x + 2*epsilon)**2 * (x-epsilon) * (k*x - epsilon) - p**2 * 27 * epsilon**4))
###############################################################################

print(equationSMR())
