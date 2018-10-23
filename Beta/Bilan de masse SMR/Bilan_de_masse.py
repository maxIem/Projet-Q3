import numpy as np
from Reaction_Vaporeformage import Vaporeformage
from Variables import getVariable
from Reaction_WaterGasShift import waterGasShift


#######################################
temperature, pression, ratio, flux = getVariable()      # Importe les variables depuis Variables.py
product_Flux_Tab = []                                   # Initialise le tableau des flux en sorite de forme [CH4, H2O, CO, H2, CO2]
#######################################


# Resous le systeme pour des parametres donnes en appelant
# Systeme.py/Vaporeformage(Temperature,Pression,Ratio,Flux)
# et reourne les flux en sorties sous forme d'array [[reaction 1], [reaction 2]]
def vaporeformageFluxSortant(Temperature,Pression,Ratio,Flux):
    sol = Vaporeformage(Temperature,Pression,Ratio,Flux)
    return np.array([ [ Flux-sol[0] , Ratio*Flux - sol[0] - sol[1] , sol[0] - sol[1] , 3*sol[0] + sol[1] ],
    [ sol[0] - sol[1] , Ratio*Flux - sol[0] - sol[1] , sol[1] , 3*sol[0] + sol[1] ] ])
#######################################


product_Flux_Tab = vaporeformageFluxSortant(temperature,pression,ratio,flux)
print('Flux apres reaction SMR')
print('Flux de CH4 en sortie %.3f [mol/s]' %(product_Flux_Tab[0][0]))
print('Flux de H2O en sortie %.3f [mol/s]' %(product_Flux_Tab[0][1]))
print('Flux de CO en sortie %.3f [mol/s]'  %(product_Flux_Tab[0][2]))
print('Flux de H2 en sortie %.3f [mol/s]'  %(product_Flux_Tab[0][3]))
print('Flux de CO2 en sortie %.3f [mol/s]' %(product_Flux_Tab[1][2]))
print('-----------------------------------')

# Applique la reaction WGS (complete) aux flux sortant
WGS_Tab = waterGasShift([product_Flux_Tab[0][2], product_Flux_Tab[0][1], product_Flux_Tab[1][2], product_Flux_Tab[0][3]])
product_Flux_Tab[0][2] = WGS_Tab[0]
product_Flux_Tab[0][1] = WGS_Tab[1]
product_Flux_Tab[1][2] = WGS_Tab[2]
product_Flux_Tab[0][3] = WGS_Tab[3]

print('Flux apres reaction WGS complete')
print('Flux de CH4 en sortie %.3f [mol/s]' %(product_Flux_Tab[0][0]))
print('Flux de H2O en sortie %.3f [mol/s]' %(product_Flux_Tab[0][1]))
print('Flux de CO en sortie %.3f [mol/s]'  %(product_Flux_Tab[0][2]))
print('Flux de H2 en sortie %.3f [mol/s]'  %(product_Flux_Tab[0][3]))
print('Flux de CO2 en sortie %.3f [mol/s]' %(product_Flux_Tab[1][2]))
print('-----------------------------------')

# Purete = H2/(H2+CH4)
print('Purete H2 %.2f%%' % (100*product_Flux_Tab[0][3]/(product_Flux_Tab[0][3] + product_Flux_Tab[0][0])))
