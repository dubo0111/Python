# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 2018

@author: DUBO

p-center model
"""
import data_generator1 as dg
#INPUT Parameters:p, cost matrix
#p,cd = dg.ins_small()
p,cd = dg.ins_big(100)
from gurobipy import *
def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        print('-----------------------')
        print(model.cbGetSolution(model._vars))

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

    # (1) Maximum cost constraints (objective): L>sum(cdx) forall i
    cdx = x.copy()
    for i in range(ni):
        for j in range(ni):
           cdx[i,j]=cd[i][j]
    m.addConstrs(
            (x.prod(cdx,i,'*') <= L for i in range(ni)),
            "epigraph")
    # (2) Constraints sum(y)=p
    m.addConstr(
            (y.sum() == p),
            "p")
    # (3) x<=y forall i,j
    m.addConstrs(
            (x[i,j] <= y[j] for i in range(ni) for j in range(ni)),
            "x<y")
    # (4) sum(x)=1 forall i
    m.addConstrs(
            (x.sum(i,'*') == 1 for i in range(ni)),
            "sumx")
    # save model and optimize
    m.write('.\model\pcenter.lp')
    m.params.OutputFlag = 0
    m._vars = m.getVars()
    m.optimize(mycallback)

    # Output
#    for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    print('Obj: %g' % m.objVal)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
