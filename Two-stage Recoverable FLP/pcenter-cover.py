# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 2018

@author: DUBO

p-center model (cover formulation)
"""
import data_generator0 as dg0
data = dg0.data_gen(50,1,3) # ni,nk,randomseed
p,cd,_,_ = data.data()
from gurobipy import *
import numpy as np
import itertools
import sys
import time


def p_center(cd,p=1):
    m = Model("p-center")
    # Number of nodes
    ni = len(cd)
    # p = 1
    # Create variables
    # x:allocations y:location L:auxiliary variable
    x = m.addVars(ni,ni,vtype=GRB.CONTINUOUS, name="x")
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
    # m.write('.\model\pcenter.lp')
    m.params.OutputFlag = 0
    # m._vars = m.getVars()
    m.optimize()
    # Output
    #    for v in m.getVars():
    #         print('%s %g' % (v.varName, v.x))
    # print('Obj: %g' % m.objVal)
    # print('Runtime: %g' % m.Runtime)
    return m.objVal,m.Runtime

def cover_p_center(a,cd1,ni,ne,p):
    # Create a new model
    m = Model("p-center-cover")

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

    # m.addConstr(y[0] == 1)

    # save model and optimize
    # m.write('.\model\pcenter.lp')
    m.params.OutputFlag = 0
    # m._vars = m.getVars()
    m.optimize()

    # Output
    #    for v in m.getVars():
    #         print('%s %g' % (v.varName, v.x))
    print('Obj: %g' % m.objVal)
    print('Runtime: %g' % m.Runtime)
    # return

# Preprocessing
cd0 = list(itertools.chain.from_iterable(cd)) # combine lists
cd1 = sorted(set(cd0)) # sort without duplicates
ni = len(cd)

t0 = time.time()
# Find UB1
y_initial = [0 for i in range(ni)] # y
y_set = set()
n = p # counter
cd_matrix = np.array(cd)
cd_matrix_1 = np.copy(cd_matrix) # deep copy
# Find y heuristically
while len(y_set) < p:
    (a,b) = np.unravel_index(cd_matrix.argmax(), cd_matrix.shape)
    cd_matrix[a,b] = 0
    y_set.update({a,b})
if len(y_set) > p:
    y_set.remove(b)
for x in y_set:
    y_initial[x] = 1
# Construct x (|y| cluster) by finding closest facility
x = [[0 for j in range(ni)] for i in range(ni)]
for j in range(ni):
    if y_initial[j] == 0:
        cd_matrix_1[:,j] = 1e8 # set j column to Big M
facility = np.argmin(cd_matrix_1, axis=1) # index of maximum value in each row
cost = cd_matrix_1[np.arange(cd_matrix_1.shape[0]),facility] # advanced index returning maximum value of each row
UB1 = max(cost) # Upper Bound UB1
LB2 = UB1/2# Lower Bound
#b = np.min(cd_matrix_1, axis=1)
for i in range(ni):
    x[i][facility[i]] = 1
print("Find UB1--- %s seconds ---" % round((time.time() - t0), 2))

t0 = time.time()
# Find UB2: Cluster
cluster = [[] for n in range(p)]
y_idx = list(y_set)
for i in range(p):
    cluster[i] = np.argwhere(facility == y_idx[i]).ravel()#.tolist()
cd_array = np.array(cd)
L_cluster = [0 for i in range(p)]
for i in range(p):
    L_cluster[i],_ = p_center(cd_array[cluster[i][:,None],cluster[i]])
UB2 = max(L_cluster)
print("Find UB2--- %s seconds ---" % round((time.time() - t0), 2))
print(UB1)
print(UB2)
print(LB2)
# using UB,LB to modify cd1
cd1 = [x for x in cd1 if x>=LB2 and x<= UB1]
ne = len(cd1)
A = [[[0 for e in range(ne)] for j in range(ni)] for i in range(ni)]
for i in range(ni):
    for j in range(ni):
        for e in range(ne):
            if cd[i][j] <= cd1[e]:
                A[i][j][e] = 1
cover_p_center(A,cd1,ni,ne,p)
obj1,t1 = p_center(cd,p)
print('Obj: %g' % obj1)
print('Runtime: %g' % t1)
