# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 08:24:56 2018

2-stage recoverable p-center model:
    Benders' decomposition (relaxed model)

@author: DUBO
"""

import data_generator1 as dg
#INPUT Parameters:p, cost matrix
#p,cd = dg.ins_small()
#p,cd = dg.ins_big(5)
p,cd,cdk,sk = dg.ins_k(3,2) #(ni,nk,sumk)
# !!!!! Make sure index match: cdk VS. v_ij(k) [k][i][j]
from gurobipy import *
#import os

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
    L = m.addVar(vtype=GRB.CONTINUOUS,obj=a1,name="L")

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

    m.write('master.lp')

    # ---------- Sub problem ----------
    # Create sub model
    m1 = Model('Sub model')

    # beta gamma delta epsilon lamda mu nu
    beta = m1.addVars(nk,ni,ub = 0,vtype=GRB.CONTINUOUS, name="beta")
    gamma = m1.addVars(nk,ni,ni,ub = 0,vtype=GRB.CONTINUOUS, name="gamma")
    delta = m1.addVars(nk,ni,ni,ub = 0,vtype=GRB.CONTINUOUS, name="delta")
    epsilon = m1.addVars(nk,ni,vtype=GRB.CONTINUOUS, name="epsilon")
    lamda = m1.addVars(nk,ni,ub = 0,vtype=GRB.CONTINUOUS, name="lambda")
    mu = m1.addVars(nk,ni,ub = 0,vtype=GRB.CONTINUOUS, name="mu")
    nu = m1.addVars(nk,vtype=GRB.CONTINUOUS, name="nu")
    # Qk are auxiliary variables, minimizing every subproblem
    Qk = m1.addVars(nk,vtype=GRB.CONTINUOUS,obj = 1, name="Qk")

    # Set sub model objective to minimize
    m1.modelSense = GRB.MAXIMIZE
    #(1) Q(k) = sum_i(sum_j(y_j*gamma_kij))+sumsum_ij((1-a_kj)*delta_kij
    #           + sum_i(epsilon_ki) + sum_j((1-y_j)*lambda_kj) +
    #           + sum_j((1-a_kj)*mu_kj) + (p+sum_j(a_kj*y_j)-sum(y_j))*nu_k   forall k\
    ##
    value_y = [1,0,0]
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
    m1.addConstrs(
             (Qk[k] == LinExpr(value_y*ni,gamma.select(k,'*','*')) + \
             LinExpr(c_delta[k],delta.select(k,'*','*')) + \
             epsilon.sum(k,'*') + LinExpr(c_lamda,lamda.select(k,'*')) + \
             LinExpr(c_mu[k],mu.select(k,'*')) + c_nu[k]*nu[k] for k in range(nk)),
             "Q(k)")
    #(2) -gamma_kij+lambda_kj+mu_kj+nu_k<=0  forall k,i,j
    m1.addConstrs(
            (-gamma[k,i,j]+lamda[k,j]+mu[k,j]+nu[k] <= 0 for k in range(nk) for i in range(ni) for j in range (ni)),
            "u")
    #(3) c_kij*d_ki*beta_ki+gamma_kij+delta_kij+epsilon_ki<=0  forall k,i,j
    m1.addConstrs(
            (cdk[k][i][j]*beta[k,i] + gamma[k,i,j] + delta[k,i,j] +epsilon[k,i] <=0 for k in range(nk) for i in range(ni) for j in range(ni)),
            "v")
    #(4) -sum_i(beta_i)<=1 forall k
    m1.addConstrs(
            (-beta.sum(k,'*') <= 1 for k in range(nk)),
            "L3")
            
    # save model and optimize
    m1.write('sub.lp')
    #m.optimize()

    #Output
#    for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    print('Obj: %g' % m.objVal)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
