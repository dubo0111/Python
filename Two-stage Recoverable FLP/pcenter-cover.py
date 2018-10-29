# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 2018

@author: DUBO

p-center model (cover formulation)
"""

import data_generator0 as dg0
data = dg0.data_gen(10,1,2)
p,cd,cdk,sk = data.data()
from gurobipy import *
import itertools

# Preprocessing
cd0 = list(itertools.chain.from_iterable(cd)) # combine lists
cd1 = sorted(set(cd0)) # sort without duplicates
ni = len(cd)
ne = len(cd1)
a = [[[0 for e in range(ne)] for j in range(ni)] for i in range(ni)]
for i in range(ni):
    for j in range(ni):
        for e in range(ne):
            if cd[i][j] <= cd1[e]:
                a[i][j][e] = 1

try:

    # Create a new model
    m = Model("p-center")


    # Create variables
    # z:ordered cost, y:location
    z = m.addVars(ne,vtype=GRB.BINARY, name="z")
    y = m.addVars(ni,vtype=GRB.BINARY, name="y")

    # Set objective to minimize
    m.modelSense = GRB.MINIMIZE
    # Minimize :\sum_e \rho_e*z_e
    m.setObjective(LinExpr(cd1,z.select()))
    # (1) \sum_j a_ije*y_j >= z_e \forall i,e
    for i in range(ni):
        for e in range(ne):
            sum_ay=0
            for j in range(ni):
                sum_ay += a[i][j][e]*y[j]
            m.addConstr(
                    (sum_ay >= z[e]),
                    'ay>e'+str(i)+str(e))
    # (2) \sum y_j = p
    m.addConstr(
            (y.sum() == p),
            'sump')
    # (3) \sum z_e = 1
    m.addConstr(
            (z.sum() == 1),
            'sumz')

    m.addConstr(y[0] == 1)
    
    # save model and optimize
    # m.write('.\model\pcenter.lp')
    # m.params.OutputFlag = 1
    # m._vars = m.getVars()
    m.optimize()

    # Output
#    for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    print('Obj: %g' % m.objVal)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
