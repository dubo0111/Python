from gurobipy import *
m = Model()
m = read('dual0.lp')
m.optimize()

#Model.read('dual.lp')
