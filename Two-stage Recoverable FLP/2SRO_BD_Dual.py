# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 08:24:56 2018

2-stage recoverable p-center model:
    Benders' decomposition (relaxed dual sub model)
    failed: dual sub model cannot be solved to optimal by gurobi. constraints violated 1e-13
@author: DUBO
"""

import data_generator1 as dg
# debug:(3,1,3):primal dual value error;
# debug:(3,1,2):dual problem error
p,cd,cdk,sk = dg.ins_k(10,10,1) #(ni,nk,randomseed)

from gurobipy import *
ni = len(cd)
nk = len(cdk)
# weights of two stages
a1 = 0.5
a2 = 1 - a1

def cal_DualObj():
    #Q(k) = sum_i(sum_j(y_j*gamma_kij))+sumsum_ij((1-a_kj)*delta_kij
    #           + sum_i(epsilon_ki) + sum_j((1-y_j)*lambda_kj) +
    #           + sum_j((1-a_kj)*mu_kj) + (p+sum_j(a_kj*y_j)-sum(y_j))*nu_k
    sum_gamma = 0
    sum_delta = 0
    sum_epsilon = sum(epsilon)
    sum_lamda = 0
    sum_mu = 0
    sum_nu = p+sum([sk[max_k][j]*value_y[j] for j in range(ni)])+sum(value_y)
    for i in range(ni):
        for j in range(ni):
            sum_gamma += gamma[i][j]*value_y[j]
            sum_delta += delta[i][j]*(1-sk[max_k][j])
    for n in range(ni):
        sum_lamda += 1-value_y[n]*lamda[n]
        sum_mu += 1-sk[max_k][n]
    dual_obj1 = sum_gamma+sum_delta+sum_epsilon+sum_lamda+sum_mu+sum_nu
    dual_obj2 = sum_gamma+sum_delta+sum_epsilon
    return dual_obj1,dual_obj2
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

# --------- Master problem ---------
# Create master model
m = Model("master model")
m.params.OutputFlag = 0
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
def sub_dual(value_y):
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
    # m1.addConstrs(
    #         (Qk[k] == LinExpr(value_y*ni,gamma.select(k,'*','*')) + \
    #         LinExpr(c_delta[k],delta.select(k,'*','*')) + \
    #         epsilon.sum(k,'*') for k in range(nk)),
    #         "Q(k)")
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
    return m1
def sub_model(value_y):
    # Create relaxed sub model
    m2 = Model('sub model')
    
    # ---------- Sub problem ----------
    # v:allocations u:location L3,eta: auxiliary variable
    v = m2.addVars(nk,ni,ni,lb=0,ub=1, vtype=GRB.CONTINUOUS, name="v")
    u = m2.addVars(nk,ni,lb=0,ub=1,vtype=GRB.CONTINUOUS, name="u")
    L3 = m2.addVars(nk,vtype=GRB.CONTINUOUS, name="L3")
    eta = m2.addVar(vtype=GRB.CONTINUOUS, obj=1, name="eta")
    m2.modelSense = GRB.MINIMIZE
    #(5) eta == sum(L3(k)) forall k
    m2.addConstr(
            (eta == L3.sum()),
            "eta=sumL")
    #(6) L3(k) >= c'd'v(k) forall i,k
    cdv = v.copy()
    for k in range(nk):
        for i in range(ni):
            for j in range(ni):
               cdv[k,i,j]=cdk[k][i][j]
    m2.addConstrs(
            (v.prod(cdv,k,i,'*') <= L3[k] for k in range(nk) for i in range(ni)),
            "beta")
    #(7) v(k) <= y + u(k) forall k,i,j
    m2.addConstrs(
            (v[k,i,j] <= value_y[j] + u[k,j] for k in range(nk) for i in range(ni) for j in range(ni)),
            "gamma")
    #(8) v(k) <= 1 - a_j(k) forall k,i,j
    m2.addConstrs(
            (v[k,i,j] <= 1 - sk[k][j] for k in range(nk) for i in range(ni) for j in range(ni)),
            "delta")
    #(9) sum(v) = 1 forall i,k
    m2.addConstrs(
            (v.sum(k,i,'*') == 1 for k in range(nk) for i in range(ni)),
            "epsilon")
    #(10) u(k) + y <= 1 forall k,j
    m2.addConstrs(
            (u[k,j] + value_y[j] <= 1 for k in range(nk) for j in range(ni)),
            "lamda")
    #(11) u(k) + a_j(k) <= 1 forall k,j
    m2.addConstrs(
            (u[k,j] + sk[k][j] <= 1 for k in range(nk) for j in range(ni)),
            "mu")
    #(12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
    ky = [0 for k in range(nk)]
    for k in range(nk):
        ky[k] = sum([sk[k][j]*value_y[j] for j in range(ni)])
    m2.addConstrs(
            (u.sum(k,'*') + sum(value_y) - ky[k] == p for k in range(nk)),
            "nu")
    return m2
#=========================================
m.optimize()
# extract value_y
value_y = []
for j in range(ni):
    y_name = ''.join(['y[',str(j),']'])
    y_temp = m.getVarByName(y_name)
    value_y.append(y_temp.x)
m1=sub_dual(value_y)
#m1.params.FeasibilityTol = 1e-9
#m1.params.OptimalityTol = 1e-9
#m1.params.MarkowitzTol = 0.0001
#m1.params.IntFeasTol = 1e-5

m1.params.ScaleFlag = 3
#m1.params.ObjScale = 100
m1.params.NumericFocus = 3
#m1.params.NormAdjust = 3
m1.params.InfUnbdInfo = 0 #1
m1.params.Quad = -1 #1
#m1.params.Sifting = -1 #2
#m1.params.SiftMethod = -1 # 2
#m1.params.SimplexPricing = -1 #3

#m1.params.Method = -1
#m1.params.AggFill = 0
#m1.params.Aggregate = 0
#m1.params.DualReductions = 1 #0
#m1.params.PreDual = 2 #2
m1.params.Presolve = 0 #2
m1.write('dual.lp')
m1.params.OutputFlag = 0
m1.optimize()
m2=sub_model(value_y)
m2.write('primal.lp')
m2.params.FeasibilityTol = 1e-9
m2.params.OptimalityTol = 1e-9
#m2.params.MarkowitzTol = 0.0001
#m2.params.IntFeasTol = 1e-5

m2.params.ScaleFlag = 3
#m2.params.ObjScale = 0
m2.params.NumericFocus = 3
#m2.params.NormAdjust = -1
#m2.params.InfUnbdInfo = 0 #1
#m2.params.Quad = -1 #1
#m2.params.Sifting = -1 #2
#m2.params.SiftMethod = -1 # 2
#m2.params.SimplexPricing = -1 #3

#m2.params.Method = -1
#m2.params.AggFill = -1  #inf
#m2.params.Aggregate = 1 #1
#m2.params.DualReductions = 0 #0
#m2.params.PreDual = -1 #2
m2.params.Presolve = 0 #2
m2.params.OutputFlag = 0
m2.optimize()

# max_k for dual
value_Qk = []
for k in range(nk):
    Qk_name = ''.join(['Qk[',str(k),']'])
    Qk_temp = m1.getVarByName(Qk_name)
    value_Qk.append(Qk_temp.x)
# maximum L3 and its index (worst k)
max_Qk = max([[v,i] for i,v in enumerate(value_Qk)])
max_k_Qk = max_Qk[1]
#max_k for primal
value_L3 = []
for k in range(nk):
    L3_name = ''.join(['L3[',str(k),']'])
    L3_temp = m2.getVarByName(L3_name)
    value_L3.append(L3_temp.x)
# maximum L3 and its index (worst k)
max_L3 = max([[v,i] for i,v in enumerate(value_L3)])
max_k = max_L3[1]
# get dual
cm1 = m2.getConstrs()
num_c = len(cm1)
dual_value = []
constrname = []
for i in range(num_c):
    dual_value.append(cm1[i].getAttr('Pi'))
    constrname.append(cm1[i].getAttr('ConstrName'))
dual = dict(zip(constrname,dual_value))
#print(dual)
gamma = [[0 for j in range(ni)] for i in range(ni)]
delta = [[0 for j in range(ni)] for i in range(ni)]
epsilon = [0 for j in range(ni)]
lamda = [0 for j in range(ni)]
mu = [0 for j in range(ni)]
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
# Benders' cut
# omega >= sumsum(gamma_k'ij*y) + sum_j(-lamda*y) +  nu*sum_j((aj(k')-1)*y)
# + sumsum((1-aj(k'))*delta_ij)+sum_i(epsilon)+sum_j(lamda)
# + sum_j(1-aj(k'))*mu_j + p*nu
#gamma_y = []
#for i in range(ni):
#    for j in range(ni):
#        gamma_y.append(gamma[i][j])
#ajk_y = []
#for j in range(ni):
#    ajk_y.append(sk[max_k][j]-1)
#c_y = LinExpr(gamma_y,y.select()*ni) - LinExpr(lamda,y.select()) \
#+ nu*LinExpr(ajk_y,y.select())
#constant_delta = 0
#for i in range(ni):
#    constant_delta += quicksum([(1-sk[max_k][j])*delta[i][j] for j in range(ni)])
#constant = quicksum(epsilon) + quicksum(lamda) + constant_delta +\
#quicksum([(1-sk[max_k][j])*mu[j] for j in range(ni)]) + p*nu
#constr_y = c_y + constant
#print(constr_y)
#
print('=============================')
print('dual:',max_Qk[0],'k=',max_Qk[1])
print('=============================')
print('primal:',max_L3[0],'k=',max_L3[1])
print('=============================')
dualobj1,dualobj2 = cal_DualObj()
print('primal-dual1:',dualobj1)
print('primal-dual2:',dualobj2)
