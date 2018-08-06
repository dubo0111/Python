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
p,cd,cdk,sk = dg.ins_k(3,2,3) #(ni,nk,randomseed)
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
        # Create relaxed sub model
        m1 = Model('Sub model')
        # ---------- Sub problem ----------
        # v:allocations u:location L3,eta: auxiliary variable
        v = m1.addVars(nk,ni,ni,lb=0,ub=float('inf'), vtype=GRB.CONTINUOUS, name="v")
        u = m1.addVars(nk,ni,lb=0,ub=float('inf'),vtype=GRB.CONTINUOUS, name="u")
        L3 = m1.addVars(nk,vtype=GRB.CONTINUOUS, name="L3")
        eta = m1.addVar(vtype=GRB.CONTINUOUS, obj=1, name="eta")
        m1.modelSense = GRB.MINIMIZE
        #(5) eta == sum(L3(k)) forall k
        m1.addConstr(
                (eta == L3.sum()),
                "eta=sumL")
        #(6) L3(k) >= c'd'v(k) forall i,k
        cdv = v.copy()
        for k in range(nk):
            for i in range(ni):
                for j in range(ni):
                   cdv[k,i,j]=cdk[k][i][j]
        m1.addConstrs(
                (v.prod(cdv,k,i,'*') <= L3[k] for k in range(nk) for i in range(ni)),
                "beta")
        #(7) v(k) <= y + u(k) forall k,i,j
        m1.addConstrs(
                (v[k,i,j] <= value_y[j] + u[k,j] for k in range(nk) for i in range(ni) for j in range(ni)),
                "gamma")
        #(8) v(k) <= 1 - a_j(k) forall k,i,j
        m1.addConstrs(
                (v[k,i,j] <= 1 - sk[k][j] for k in range(nk) for i in range(ni) for j in range(ni)),
                "delta")
        #(9) sum(v) = 1 forall i,k
        m1.addConstrs(
                (v.sum(k,i,'*') == 1 for k in range(nk) for i in range(ni)),
                "epsilon")
        #(10) u(k) + y <= 1 forall k,j
        m1.addConstrs(
                (u[k,j] + value_y[j] <= 1 for k in range(nk) for j in range(ni)),
                "lamda")
        #(11) u(k) + a_j(k) <= 1 forall k,j
        m1.addConstrs(
                (u[k,j] + sk[k][j] <= 1 for k in range(nk) for j in range(ni)),
                "mu")
        #(12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
        ky = [0 for k in range(nk)]
        for k in range(nk):
            ky[k] = sum([sk[k][j]*value_y[j] for j in range(ni)])
        m1.addConstrs(
                (u.sum(k,'*') + sum(value_y) - ky[k] == p for k in range(nk)),
                "nu")
        return(m1)
    def update_master(m,m1):
        # extract maximum L3
        value_L3 = []
        for k in range(nk):
            L3_name = ''.join(['L3[',str(k),']'])
            L3_temp = m.getVarByName(L3_name)
            value_L3.append(y_temp.x)
        # maximum L3 and its index (worst k)
        max_k = max([[v,i] for i,v in enumerate(value_L3)])[1]
        # get dual variable value for gerating Benders cut
        cm1 = m1.getConstrs()
        num_c = len(cm1)
        dual_value = []
        constrname = []
        for i in range(num_c):
            dual_value.append(cm1[i].getAttr('Pi'))
            constrname.append(cm1[i].getAttr('ConstrName'))
        dual = dict(zip(constrname,dual_value))
        gamma = [[0 for j in range(ni)] for i in range(ni)]
        delta = gamma.copy()
        epsilon = [0 for j in range(ni)]
        lamda = epsilon.copy()
        mu = epsilon.copy()
        for i in range(ni):
            for j in range(ni):
                gamma_name = ''.join(['gamma[',str(max_k),',',str(i),',',str(j),']'])
                delta_name = ''.join(['delta[',str(max_k),',',str(i),',',str(j),']'])
                gamma[i][j] = dual[gamma_name]
                delta[i][j] = dual[delta_name]
        for n in range(ni):
            epsilon_name = ''.join(['epsilon[',str(max_k),',',str(n),']'])
            lamda_name = ''.join(['lamda[',str(max_k),',',str(n),']'])
            mu_name = ''.join(['mu[',str(max_k),',',str(n),']'])
            epsilon[n] = dual[epsilon_name]
            lamda[n] = dual[lamda_name]
            mu[n] = dual[mu_name]
        nu_name = ''.join(['nu[',str(max_k),']'])
        nu = dual[nu_name]
        if iteration == 1:
            omega = m1.addVar(vtype=GRB.CONTINUOUS, obj=1, name="omega")
        ''' Benders' cut
         omega >= sumsum(gamma_k'ij*y) + sum_j(-lamda*y) +  nu*sum_j((aj(k')-1)*y)
         + sumsum((1-aj(k'))*delta_ij)+sum_i(epsilon)+sum_j(lamda)
         + sum_j(1-aj(k'))*mu_j + p*nu
        '''
        gamma_y = []
        for i in range(ni):
            for j in range(ni):
                gamma_y.append(gamma[i][j])
        ajk_y = []
        for j in range(ni):
            ajk_y.append(sk[max_k][j]-1)
        aaa=0
        c_y = LinExpr(gamma_y,y.select()*ni) - LinExpr(lamda,y.select()) \
        + nu*LinExpr(ajk_y,y.select())
        constant_delta = 0
        for i in range(ni):
            constant_delta += quicksum([(1-sk[max_k][j])*delta[i][j] for j in range(ni)])
        constant = quicksum(epsilon) + quicksum(lamda) + constant_delta +\
        quicksum([(1-sk[max_k][j])*mu[j] for j in range(ni)])
        m.addConstr(omega >= c_y + constant)
        return m
    def update_sub(m1,value_y):
        m1=1
        return m1
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
            m = update_master(m,m1)
        filename= ''.join(['.\model\master(',str(iteration),').lp'])
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
        filename = ''.join(['.\model\sub(',str(iteration),').lp'])
        m1.write(filename)
        m1.optimize()
        # get dual variable value for gerating Benders cut
        cm1 = m1.getConstrs()
        num_c = len(cm1)
        dual_value = []
        constrname = []
        for i in range(num_c):
            dual_value.append(cm1[i].getAttr('Pi'))
            constrname.append(cm1[i].getAttr('ConstrName'))
        dual = dict(zip(constrname,dual_value))

        # update UB = ;

        iteration += 1 #
        #
        gap = (UB-LB)/LB
        #Output
        # for v in m.getVars():
        #      print('%s %g' % (v.varName, v.x))
        # print('Obj: %g' % m.objVal)
except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
