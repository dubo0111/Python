# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 2019

@author: DUBO

* Experiments
* Conditional p-center model (cover formulation)
* Giving value_y and sk (and nk scenarios):
* selected nodes (variable each iteration) and forbidden nodes (same each iteration)
"""
import data_generator0 as dg0
from gurobipy import *
import numpy as np
import itertools
import sys
import time
import math
import copy

# Parameters
output = 0

def pqcenter_compare(ni,p,rnd=None,survive=0):

    #---------------------- Problem setting -----------------------
    # p = round(ni / 2) # p
    if p == 1:
        sum_k = 1
    else:
        sum_k = int(round(p)*0.5) # sum_k
    data = dg0.data_gen(ni,1,p,sum_k,rnd) # ni,nk,randomseed
    _,cd,_,sk = data.data()
    if survive == 0:
        survive = round(p/5)+1 # survived facilities
    # Input value_y & sk
    value_y = [0 for j in range(ni)]
    # for j in range(math.ceil(p/2)):
    #     value_y[j] = 1
    y_0 = np.random.choice(range(ni), survive, replace=False)
    for x in y_0:
        value_y[x] = 1
    # sk[0] = [0 for j in range(ni)] #

    # a_ije
    # A = [[[0 for e in range(ne)] for j in range(ni)] for i in range(ni)]
    # for i in range(ni):
    #     for j in range(ni):
    #         for e in range(ne):
    #             if sk[j] == 0: # disrupted node so cannot cover??
    #                 if cd[i][j] <= cd1[e]:
    #                     A[i][j][e] = 1

    def p_center(cd,p=1,LP=0,value_y=0):
        m = Model("p-center")
        # Number of nodes
        ni = len(cd)
        # p = 1
        # Create variables
        # x:allocations y:location L:auxiliary variable
        x = m.addVars(ni,ni,vtype=GRB.CONTINUOUS, name="x")
        if LP == 0:
            y = m.addVars(ni,vtype=GRB.BINARY, name="y")
        elif LP == 1:
            y = m.addVars(ni,vtype=GRB.CONTINUOUS, name="y")
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
        if LP == 0 and p!=1:
            m.params.OutputFlag = output
        else:
            m.params.OutputFlag = 0
        # m._vars = m.getVars()
        if value_y != 0:
            m.update()
            for j in range(ni):
                if value_y[j] == 1:
                    if sk[0][j] != 1:
                        y_name = ''.join(['y[',str(j), ']'])
                        m.getVarByName(y_name).lb = 1
                if sk[0][j] == 1:
                    y_name = ''.join(['y[',str(j), ']'])
                    m.getVarByName(y_name).ub = 0
        # m.params.NumericFocus = 2
        m.Params.TimeLimit = 1000
        m.optimize()
        # Output
        #    for v in m.getVars():
        #         # print('%s %g' % (v.varName, v.x))
        # # print('Obj: %g' % m.objVal)
        # # print('Runtime: %g' % m.Runtime)
        if LP == 0:
            Gap = m.MIPGap
        else:
            Gap = 0
        return m.objVal,m.Runtime,Gap

    def cover_p_center(a,cd1,ni,ne,p,value_y=0):
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
        m.params.OutputFlag = output
        m.Params.TimeLimit = 1000
        # m.params.NumericFocus = 2
        # variable fixing
        if value_y != 0:
            m.update()
            for j in range(ni):
                if value_y[j] == 1:
                    if sk[0][j] != 1:
                        y_name = ''.join(['y[',str(j), ']'])
                        m.getVarByName(y_name).lb = 1
                if sk[0][j] == 1:
                    y_name = ''.join(['y[',str(j), ']'])
                    m.getVarByName(y_name).ub = 0
        # m._vars = m.getVars()
        m.optimize()

        # Output
        #    for v in m.getVars():
        #         # print('%s %g' % (v.varName, v.x))
        # print('Obj: %g' % m.objVal)
        # print('Runtime: %g' % m.Runtime)
        return m.objVal,m.Runtime,m.MIPGap

    # Preprocessing
    cd0 = list(itertools.chain.from_iterable(cd)) # combine lists
    cd1 = sorted(set(cd0)) # sort without duplicates
    ni = len(cd)

    t0 = time.time()
    # Find UB1
    y_initial = copy.deepcopy(value_y)
    # Construct sets
    y_set = set()
    for j in range(ni):
        if value_y[j] == 1:
            y_set.update({j})
    sk_set = set()
    for j in range(ni):
        if sk[0][j] == 1:
            sk_set.update({j})
    clear_list = list(range(ni))
    clear_set = set(clear_list)
    clear_set = clear_set.difference(y_set)
    clear_set = clear_set.difference(sk_set)
    # n = p # counter
    cd_matrix = np.array(cd)
    cd_matrix_1 = np.copy(cd_matrix) # deep copy
    # **********************************************************
    # Find y heuristically
    # Plan 1 : (L no more than 2 times of optimality)
    # (a,b) = np.unravel_index(cd_matrix.argmax(), cd_matrix.shape) # index of maximum cost
    # y_set.update({a,b})
    # cd_matrix[a,b] = 0

    # UB1  Algorithm 1: find maximum cost to node in y_set
    # while len(y_set) < p:
    #     y_curr = list(y_set)
    #     cd_temp = cd_matrix[y_curr,:]
    #     (a,b) = np.unravel_index(cd_temp.argmax(), cd_temp.shape)
    #     cd_matrix[y_curr[a],b] = 0
    #     if sk[0][b] != 1:
    #         y_set.update({b})

    # UB1 Algorithm 2
    while len(y_set) < p:
        (a,b) = np.unravel_index(cd_matrix.argmax(), cd_matrix.shape) # index of maximum cost
        if a in clear_set and b in clear_set:
            y_set.update({a,b})
        elif a in clear_set and b not in clear_set:
            y_set.update({a})
        elif a not in clear_set and b in clear_set:
            y_set.update({b})
        cd_matrix[a,b] = 0
    if len(y_set) > p:
        y_set.remove(b)
    #  Plan 2 (more than 2 times of optimality)?

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

    #b = np.min(cd_matrix_1, axis=1)
    for i in range(ni):
        x[i][facility[i]] = 1
    # print("Find UB1--- %s seconds ---" % round((time.time() - t0), 5))

    t1 = time.time()
    # Find UB2: Cluster
    cluster = [[] for n in range(p)]
    y_idx = list(y_set)
    for i in range(p):
        cluster[i] = np.argwhere(facility == y_idx[i]).ravel()#.tolist()
    cd_array = np.array(cd)
    L_cluster = [0 for i in range(p)]
    for i in range(p):
        L_cluster[i],_,_ = p_center(cd_array[cluster[i][:,None],cluster[i]])
    UB2 = max(L_cluster)
    # print("Find UB2--- %s seconds ---" % round((time.time() - t1), 5))

    t2 = time.time()
    # LB2 = UB1/2  # Wrong Lower Bound
    LB2,_,_ = p_center(cd,p,1)
#    LB2=0
    # print("Find LB--- %s seconds ---" % round((time.time() - t2), 5))

    # **********************************************************
    # print(UB1)
    # print(UB2)
    # print(LB2)

    # using UB,LB to modify cd1
    cd1 = [x for x in cd1 if x>=LB2 and x<= UB2]
    ne = len(cd1)
    A = [[[0 for e in range(ne)] for j in range(ni)] for i in range(ni)]
    for i in range(ni):
        for j in range(ni):
            for e in range(ne):
                if cd[i][j] <= cd1[e]:
                    A[i][j][e] = 1
    obj0,T0,Gap0 = cover_p_center(A, cd1,ni,ne,p,value_y) #
    obj1,T1,Gap1 = p_center(cd,p,0,value_y) #
    # print('Obj: %g' % obj1)
    # print('Runtime: %g' % T1)
    # print(y_set)
    return obj0,T0,Gap0,obj1,T1,Gap1
