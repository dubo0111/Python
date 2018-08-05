# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 08:24:56 2018

2-stage recoverable p-center model:
    Benders' decomposition (relaxed dual sub model)
    failed: dual sub model cannot be solved to optimal by gurobi. constraints violated 1e-13
@author: DUBO
"""

import data_generator1 as dg
#INPUT Parameters:p, cost matrix, cost matrix of each scenarios, disruption scenarios
#p,cd = dg.ins_small()
#p,cd = dg.ins_big(5)
p,cd,cdk,sk = dg.ins_k(3,2,3) #(ni,nk,randomseed)
# !!!!! Make sure index match: cdk V.S. v_ij(k) [k][i][j]
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
        m1 = Model('sub model')
        # beta gamma delta epsilon lamda mu nu
        beta = m1.addVars(nk,ni,ub = 0,lb = -float('inf'), vtype=GRB.CONTINUOUS, name="beta")
        gamma = m1.addVars(nk,ni,ni,ub = 0,lb = -float('inf'), vtype=GRB.CONTINUOUS, name="gamma")
        delta = m1.addVars(nk,ni,ni,ub = 0,lb = -float('inf'), vtype=GRB.CONTINUOUS, name="delta")
        epsilon = m1.addVars(nk,ni,lb = 0,vtype=GRB.CONTINUOUS, name="epsilon")
        lamda = m1.addVars(nk,ni,ub = 0,lb = -float('inf'), vtype=GRB.CONTINUOUS, name="lambda")
        mu = m1.addVars(nk,ni,ub = 0,lb = -float('inf'), vtype=GRB.CONTINUOUS, name="mu")
        nu = m1.addVars(nk,vtype=GRB.CONTINUOUS, name="nu")
        # Qk are auxiliary variables, maximize every subproblem
        Qk = m1.addVars(nk,vtype=GRB.CONTINUOUS,obj = 1, name="Qk")
        # Set sub model objective to maximize
        m1.modelSense = GRB.MAXIMIZE
        #(1) Q(k) = sum_i(sum_j(y_j*gamma_kij))+sumsum_ij((1-a_kj)*delta_kij
        #           + sum_i(epsilon_ki) + sum_j((1-y_j)*lambda_kj) +
        #           + sum_j((1-a_kj)*mu_kj) + (p+sum_j(a_kj*y_j)-sum(y_j))*nu_k   forall k
        c_delta,c_lamda,c_mu,c_nu = c_constr1(value_y) # update coeff
        m1.addConstrs(
                 (Qk[k] == LinExpr(value_y*ni,gamma.select(k,'*','*')) + \
                 LinExpr(c_delta[k],delta.select(k,'*','*')) + \
                 epsilon.sum(k,'*') + LinExpr(c_lamda,lamda.select(k,'*')) + \
                 LinExpr(c_mu[k],mu.select(k,'*')) + c_nu[k]*nu[k] for k in range(nk)),
                 "Q(k)")
        #(2) -sum_i(gamma_kij)+lambda_kj+mu_kj+nu_k<=0  forall k,j
        m1.addConstrs(
                (-gamma.sum(k,'*',j)+lamda[k,j]+mu[k,j]+nu[k] <= 0 for k in range(nk) for j in range (ni)),
                "u")
        #(3) c_kij*d_ki*beta_ki+gamma_kij+delta_kij+epsilon_ki<=0  forall k,i,j
        m1.addConstrs(
                (cdk[k][i][j]*beta[k,i] + gamma[k,i,j] + delta[k,i,j] + epsilon[k,i] <= 0 \
                for k in range(nk) for i in range(ni) for j in range(ni)),
                "v")
        #(4) -sum_i(beta_i)<=1 forall k
        m1.addConstrs(
                (-beta.sum(k,'*') <= 1 for k in range(nk)),
                "L3")
        # m1.addConstr((beta[0,0] == -0.40801750478276144))
        # m1.addConstr((beta[0,2] == -0.5919824952172386))
        # m1.addConstr((gamma[0,0,0] == -927.3994470233015))
        # m1.addConstr((gamma[0,2,1] == -927.3994470233015))
        # m1.addConstr((delta[0,0,2] == -574.579597469973))
        # m1.addConstr((delta[0,2,2] == -1836.2115594007862))
        # m1.addConstr((epsilon[0,0] == 927.3994470233015))
        # m1.addConstr((epsilon[0,2] == 1836.2115594007862))
        # m1.addConstr((nu[0] == -927.3994470233015))
        return(m1)
    def update_master(m,subx):
        m=1
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
#        m1.write(filename)
#        m1.params.FeasibilityTol = 1e-6
#        m1.params.OptimalityTol = 1e-6
        #m1.params.MarkowitzTol = 0.0001
        #m1.params.Method = -1
        #m1.params.Aggregate = 0
#        m1.params.ScaleFlag = 3
#        m1.params.ObjScale = 100
        m1.optimize()
        #
        value_Q = []
        for k in range(nk):
            name = ''.join(['Qk[',str(k),']'])
            temp = m1.getVarByName(name)
            value_Q.append(temp.x)
        print(max(value_Q)*0.5+LB)

        # print dual variable value
        cm1 = m1.getConstrs()
        num_c = len(cm1)
        dual_value = []
        constrname = []
        for i in range(num_c):
            dual_value.append(cm1[i].getAttr('Pi'))
            constrname.append(cm1[i].getAttr('ConstrName'))
        dual = dict(zip(constrname,dual_value))
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
