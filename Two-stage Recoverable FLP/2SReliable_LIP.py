# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 2018

@author: DUBO

p-center model
"""
import data_generator1 as dg
#INPUT Parameters:p, cost matrix
#p,cd = dg.ins_small()
#p,cd = dg.ins_big(5)
p,cd,cdk,sk = dg.ins_k(6,5) #(ni,nk,sumk)
# !!!!! Make sure index match: cdk VS. v_ij(k) [k][i][j]
from gurobipy import *

try:

    # Create a new model
    m = Model("p-center")

    # Number of nodes
    ni = len(cd)
    nk = len(cdk)

    # Number of scenarios
    #nk = 1
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
    # ---------- Sub problem ----------
    # v:allocations u:location L3,eta: auxiliary variable
    v = m.addVars(nk,ni,ni,vtype=GRB.BINARY, name="v")
    u = m.addVars(nk,ni,vtype=GRB.BINARY, name="u")
    L3 = m.addVars(nk,vtype=GRB.CONTINUOUS, name="L3")
    eta = m.addVar(vtype=GRB.CONTINUOUS, obj=a2, name="eta")

    #(5) eta >= L3(k) forall k
    m.addConstrs(
            (eta >= L3[k] for k in range(nk)),
            "eta>k")
    #(6) L3(k) >= c'd'v(k) forall i,k
    cdv = v.copy()
    for k in range(nk):
        for i in range(ni):
            for j in range(ni):
               cdv[k,i,j]=cdk[k][i][j]
    m.addConstrs(
            (v.prod(cdv,k,i,'*') <= L3[k] for k in range(nk) for i in range(ni)),
            "sumcdv<L3k")
    #(7) v(k) <= y + u(k) forall k,i,j
    m.addConstrs(
            (v[k,i,j] <= y[j] + u[k,j] for k in range(nk) for i in range(ni) for j in range(ni)),
            "v<y+u")
    #(8) v(k) <= 1 - a_j(k) forall k,i,j
    m.addConstrs(
            (v[k,i,j] <= 1 - sk[k][j] for k in range(nk) for i in range(ni) for j in range(ni)),
            "v<1-ak")
    #(9) sum(v) = 1 forall i,k
    m.addConstrs(
            (v.sum(k,i,'*') == 1 for k in range(nk) for i in range(ni)),
            "sumv")
    #(10) u(k) + y <= 1 forall k,j
    m.addConstrs(
            (u[k,j] + y[j] <= 1 for k in range(nk) for j in range(ni)),
            "u+y<1")
    #(11) u(k) + a_j(k) <= 1 forall k,j
    m.addConstrs(
            (u[k,j] + sk[k][j] <= 1 for k in range(nk) for j in range(ni)),
            "u+ak<1")
    #(12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
    m.addConstrs(
            (u.sum(k,'*') + y.sum() - LinExpr(sk[k],y.select()) == p for k in range(nk)),
            "2S-p")
    # save model and optimize
    m.write('2S_Recovarable_pcenter.lp')
    m.optimize()

    #Output
#    for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    print('Obj: %g' % m.objVal)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')