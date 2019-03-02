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
from gurobipy import *


#INPUT Parameters:p, cost matrix
def one_stage(ni,sk,rnd,a1):
    data = dg0.data_gen(ni,sk,rnd)
    cd,cdk,sk,_ = data.illustrative()
    p=2

    start_time = time.time()

    try:

        # Create a new model
        m = Model("p-center")

        # Number of nodes
        ni = len(cd)
        # nk = len(cdk)
        nk = ni

        # Number of scenarios
        #nk = 1
        # weights of two stages
#        a1 = 0.4
        a2 = 1 - a1

        # --------- Master problem ---------
        # Create variables
        # x:allocations y:location L:auxiliary variable
        x = m.addVars(ni,ni,vtype=GRB.BINARY, name="x")
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
        m.params.OutputFlag = 0
        m.optimize()
    #    print('========================')
    #    print('L:',L.x)

        m1 = Model("p-center")
        # ---------- Sub problem ----------
        # v:allocations u:location L3,eta: auxiliary variable
        v = m1.addVars(nk,ni,ni,vtype=GRB.BINARY, name="v")
        u = m1.addVars(nk,ni,vtype=GRB.BINARY, name="u")
        L3 = m1.addVars(nk,vtype=GRB.CONTINUOUS, obj=0.00001,name="L3")
        eta = m1.addVar(vtype=GRB.CONTINUOUS, obj=a2, name="eta")

        #(5) eta >= L3(k) forall k
        m1.addConstrs(
                (eta >= L3[k] for k in range(nk)),
                "eta>k")
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
                (v[k,i,j] <= y[j].x + u[k,j] for k in range(nk) for i in range(ni) for j in range(ni)),
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
                (u[k,j] + y[j].x <= 1 for k in range(nk) for j in range(ni)),
                "u+y<1")
        #(11) u(k) + a_j(k) <= 1 forall k,j
        m1.addConstrs(
                (u[k,j] + sk[k][j] <= 1 for k in range(nk) for j in range(ni)),
                "u+ak<1")
        #(12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
        sumy = 0
        for j in range(ni):
            sumy += y[j].x
        for k in range(nk):
            sumky = 0
            for j in range(ni):
                sumky += y[j].x * sk[k][j]
            m1.addConstr(
                    (u.sum(k,'*') + sumy - sumky == p),
                    "2S-p")
        # save model and optimize
    #    m.write(".\model\LIP.lp")
        m1.params.OutputFlag = 0
        m1.optimize()

        #Output
    #    for v in m.getVars():
    #         print('%s %g' % (v.varName, v.x))
        value_L = m.getVarByName('L')
        value_eta = m1.getVarByName('eta')
    #    print('Obj: %g' % m.objVal)

    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
    print('=================LIP SOLUTION==================')
    print('1st stage:',L.x)
    print('2nd stage:',eta.x)
    print('obj:',a1*L.x+a2*eta.x)

    # print("--- %s seconds ---" % round((time.time() - start_time),2))

    value_y = []
    for j in range(ni):
        y_name = ''.join(['y[',str(j),']'])
        value_y.append(m.getVarByName(y_name).x)
    value_x = [[0 for j in range(ni)] for i in range(ni)]
    for i in range(ni):
        for j in range(ni):
            x_name = ''.join(['x[',str(i),',',str(j),']'])
            value_x[i][j] = m.getVarByName(x_name).x
    value_u =  [[0 for j in range(ni)] for k in range(nk)]
    for k in range(nk):
        for i in range(ni):
            u_name = ''.join(['u[',str(k),',',str(i),']'])
            value_u[k][i] = m1.getVarByName(u_name).x
    value_v =[[[0 for j in range(ni)] for i in range(ni)] for k in range(nk)]
    for k in range(nk):
        for i in range(ni):
            for j in range(ni):
                v_name = ''.join(['v[',str(k),',',str(i),',',str(j),']'])
                value_v[k][i][j]= m1.getVarByName(v_name).x
    value_L3 = []
    for k in range(nk):
        L3_name = ''.join(['L3[',str(k),']'])
        value_L3.append(m1.getVarByName(L3_name).x)

    return L.x, eta.x,cd,value_y,value_x,value_u,value_v,value_L3
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
