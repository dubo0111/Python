# -*- coding: utf-8 -*-
"""
TEST TEST
CHECKING dual subproblem
Similiar to C&CG
"""

import data_generator1 as dg
#INPUT Parameters:p, cost matrix, cost matrix of each scenarios, disruption scenarios
#p,cd = dg.ins_small()
#p,cd = dg.ins_big(5)
p,cd,cdk,sk = dg.ins_k(3,1,5) #(ni,nk,randomseed)
# !!!!! Make sure index match: cdk VS. v_ij(k) [k][i][j]
from gurobipy import *

try:
    # Create master model
    m = Model("master model")

    # Number of nodes
    ni = len(cd)
    nk = len(cdk)

    # weights of two stages
    a1 = 0.5
    a2 = 1 - a1

    # --------- Master problem ---------
    # Create variables
    # x:allocations y:location L:auxiliary variable
    x = m.addVars(ni,ni,vtype=GRB.BINARY, name="x")
    y = m.addVars(ni,vtype=GRB.BINARY, name="y")
    L = m.addVar(vtype=GRB.CONTINUOUS,obj = a1, name="L")

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
    # ---------- Sub problem ----------
    def sub_model(value_y):
        # Create sub model
        m1 = Model('Sub model')
        # ---------- Sub problem ----------
        # v:allocations u:location L3,eta: auxiliary variable
        v = m1.addVars(nk,ni,ni,lb=0,ub=float('inf'), vtype=GRB.CONTINUOUS, name="v")
        u = m1.addVars(nk,ni,lb=0,ub=float('inf'),vtype=GRB.CONTINUOUS, name="u")
        L3 = m1.addVars(nk,vtype=GRB.CONTINUOUS, name="L3")
        eta = m1.addVar(vtype=GRB.CONTINUOUS, obj=a2, name="eta")
        m1.modelSense = GRB.MINIMIZE
        #(5) eta == sum(L3(k)) forall k
        m1.addConstr(
                (eta == L3.sum()),
                "eta=sumk")
        #(6) L3(k) >= c'd'v(k) forall i,k
        cdv = v.copy()
        for k in range(nk):
            for i in range(ni):
                for j in range(ni):
                   cdv[k,i,j]=cdk[k][i][j]
        m1.addConstrs(
                (v.prod(cdv,k,i,'*') <= L3[k] for k in range(nk) for i in range(ni)),
                "sumcdv<L3k")
        #(7) v(k) <= y + u(k) forall k,i,j
        m1.addConstrs(
                (v[k,i,j] <= value_y[j] + u[k,j] for k in range(nk) for i in range(ni) for j in range(ni)),
                "v<y+u")
        #(8) v(k) <= 1 - a_j(k) forall k,i,j
        m1.addConstrs(
                (v[k,i,j] <= 1 - sk[k][j] for k in range(nk) for i in range(ni) for j in range(ni)),
                "v<1-ak")
        #(9) sum(v) = 1 forall i,k
        m1.addConstrs(
                (v.sum(k,i,'*') == 1 for k in range(nk) for i in range(ni)),
                "sumv")
        #(10) u(k) + y <= 1 forall k,j
        m1.addConstrs(
                (u[k,j] + value_y[j] <= 1 for k in range(nk) for j in range(ni)),
                "u+y<1")
        #(11) u(k) + a_j(k) <= 1 forall k,j
        m1.addConstrs(
                (u[k,j] + sk[k][j] <= 1 for k in range(nk) for j in range(ni)),
                "u+ak<1")
        #(12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
        ky = [0 for k in range(nk)]
        for k in range(nk):
            ky[k] = sum([sk[k][j]*value_y[j] for j in range(ni)])
        m1.addConstrs(
                (u.sum(k,'*') + sum(value_y) - ky[k] == p for k in range(nk)),
                "2S-p")
        return(m1)
    def update_master(m,subx):
        m = 1
        return m
    def update_sub(m1,value_y):
        # remove old constraints
        for k in range(nk):
            constr_name = ''.join(['Q(k)[',str(k),']'])
            m1.remove(m1.getConstrByName(constr_name))
        # add new constraints
        c_delta,c_lamda,c_mu,c_nu = c_constr1(value_y)
        m1.addConstrs(
                 (Qk[k] == LinExpr(value_y*ni,gamma.select(k,'*','*')) + \
                 LinExpr(c_delta[k],delta.select(k,'*','*')) + \
                 epsilon.sum(k,'*') + LinExpr(c_lamda,lamda.select(k,'*')) + \
                 LinExpr(c_mu[k],mu.select(k,'*')) + c_nu[k]*nu[k] for k in range(nk)),
                 "Q(k)")
        return m1
    # update coeff
    def c_constr1(value_y):
        c_delta = [[0 for i in range(ni*ni)] for k in range(nk)]
        for k in range(nk):
            c_delta[k] = [1-sk[k][j] for j in range(ni)]*ni
        c_lamda = [1-value_y[j] for j in range(ni)]
        c_mu = [[0 for i in range(ni)] for k in range(nk)]
        for k in range(nk):
            c_mu[k] = [1-sk[k][j] for j in range(ni)]
        c_nu=[]
        for k in range(nk):
            c_nu.append(p + sum([sk[k][j]*value_y[j] for j in range(ni)]) - sum(value_y))
        return c_delta,c_lamda,c_mu,c_nu
    # ----------Benders' Decompisition----------
    iteration = 0
    # Set upperbound: UB and lowerbound: LB
    LB = -float('inf')
    UB = float('inf')
    gap = 1
    while gap >= 1e-5: # stop criteria
        # --- Solve master problem ---
        # update master model m. Adding a new constraint in each iteration.
        if iteration != 0:
            m = update_master(m,subx)
        filename= ''.join(['.\model\_testmaster(',str(iteration),').lp'])
        m.write(filename)
        m.optimize()
        print('........................................\
              ........................................')
        # extract value_y
        value_y = []
        for j in range(ni):
            y_name = ''.join(['y[',str(j),']'])
            y_temp = m.getVarByName(y_name)
            value_y.append(y_temp.x)
        # extract L
        var_L = m.getVarByName('L')
        value_L = var_L.x
        # update LB = objective value of master problem
        obj_master = m.getObjective()
        LB = obj_master.getValue()
        # --- Solve sub problem ---
        # Input: value_y; Output value: variable values of the worst scenarios to construct a new cut
        # update sub model m1
        if iteration == 0:
            m1 = sub_model(value_y)
        else:
            m1 = update_sub(m1,value_y)
        filename = ''.join(['.\model\_testsub(',str(iteration),').lp'])
        m1.write(filename)
        m1.optimize()
        #
        value_Q = []
        for k in range(nk):
            name = ''.join(['L3[',str(k),']'])
            temp = m1.getVarByName(name)
            value_Q.append(temp.x)
        print(max(value_Q)*0.5+LB)

        # extract subproblem variables: subx
        subx = m1.getVars()
        # update UB = ;

        iteration += 1 #
        #
        gap = -(UB-LB)/LB
        #Output
        # for v in m.getVars():
        #      print('%s %g' % (v.varName, v.x))
        # print('Obj: %g' % m.objVal)
except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
