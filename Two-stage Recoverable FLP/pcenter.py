# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 2018

@author: DUBO

p-center model
"""
import data_generator as dg
p,cd = dg.ins_small() #INPUT Parameters:p, cost matrix
from gurobipy import *

try:

    # Create a new model
    m = Model("p-center")

    # Number of nodes
    ni = len(cd)

    # Create variables
    # x:allocations y:location L:auxiliary variable
    x = m.addVars(ni,ni,vtype=GRB.BINARY, name="x")
    y = m.addVars(ni,vtype=GRB.BINARY, name="y")
    L = m.addVar(vtype=GRB.CONTINUOUS,obj=1,name="L")

    # Set objective to minimize
    m.modelSense = GRB.MINIMIZE

    # Maximum cost constraints (objective): L>sum(cdx) forall i
#    m.addConstrs(
#        (x.sum(i,'*') <= L for i in range(ni)),
#        "Capacity")
    cdx=[]
    for j in range(ni):
        cdx.append(x[])
#    m.addConstrs(
#        (sum(x['*',j]*cd['*',j] for j in range(ni)) <= L for i in range(ni)),
#        "objective")
    
    m.write('pcenter.lp')

    m.optimize()
#
#     for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
#
#     print('Obj: %g' % m.objVal)
#
except GurobiError as e:
     print('Error code ' + str(e.errno) + ": " + str(e))
#
except AttributeError:
     print('Encountered an attribute error')
