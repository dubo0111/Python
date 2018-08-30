# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 2018

2-stage recoverable p-center model:
    Linear Reformulation

@author: DUBO
"""
#%reset -f
import time
#import data_generator1 as dg
import data_generator0 as dg0
#INPUT Parameters:p, cost matrix
data = dg0.data_gen(100,2,2)
p,cd,cdk,sk = data.data()

from gurobipy import *
start_time = time.time()

try:

    # Create a new model
    m = Model("p-center")

    # Number of nodes
    ni = len(cd)
    nk = len(cdk)

    # Number of scenarios
    #nk = 1
    # weights of two stages
    a1 = 0.4
    a2 = 1 - a1

    # --------- Master problem ---------
    # Create variables
    # x:allocations y:location L:auxiliary variable
    x = m.addVars(ni,ni,vtype=GRB.CONTINUOUS, name="x")
    y = m.addVars(ni,vtype=GRB.BINARY, name="y")
    L = m.addVar(vtype=GRB.CONTINUOUS,obj=a1,name="L")

    # Set objective to minimize
    m.modelSense = GRB.MINIMIZE
#    m.params.OutputFlag = 0
    m.params.Presolve = 0
#    m.params.ScaleFlag = 3
#    m.params.NumericFocus = 3
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

    m.optimize()
    print('========================')
    print('L:',L.x)

    # ---------- Sub problem ----------
    # v:allocations u:location L3,eta: auxiliary variable
    v = m.addVars(nk,ni,ni,vtype=GRB.CONTINUOUS, name="v")
#    u = m.addVars(nk,ni,vtype=GRB.BINARY, name="u")
    u = m.addVars(nk,ni,vtype=GRB.CONTINUOUS, name="u")
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
#    m.write(".\model\LIP.lp")
    m.optimize()

    #Output
#    for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    value_L = m.getVarByName('L')
    value_eta = m.getVarByName('eta')
#    print('Obj: %g' % m.objVal)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
print('=================LIP SOLUTION==================')
print('1st stage:',L.x)
print('2nd stage:',eta.x)
print('obj:',a1*L.x+a2*eta.x)
print("--- %s seconds ---" % round((time.time() - start_time),2))

# value_u = []
# for k in range(nk):
#     for i in range(ni):
#         u_name = ''.join(['u[',str(k),',',str(i),']'])
#         value_u.append(m.getVarByName(u_name))
# for i in value_u:
#     if i.x != 0:
#         print(i.x)
#
#m.reset()
#try:
#    for k in range(nk):
#        for j in range(ni):
#            u_name=''.join(['u[',str(k),',',str(j),']'])
#            m.getVarByName(u_name).setAttr('vtype',GRB.BINARY)
#    m.optimize()
#    value_L = m.getVarByName('L')
#    value_eta = m.getVarByName('eta')
#except GurobiError as e:
#    print('Error code ' + str(e.errno) + ": " + str(e))
#except AttributeError:
#    print('Encountered an attribute error')
#print('=================LIP SOLUTION With MIP Subproblem==================')
#print('1st stage:',L.x)
#print('2nd stage:',eta.x)
#print('obj:',a1*L.x+a2*eta.x)
#print("--- %s seconds ---" % round((time.time() - start_time),2))
