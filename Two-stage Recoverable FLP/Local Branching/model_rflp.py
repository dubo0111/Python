# class for:
# bulid or update specific models
# get variables,dual variables
# VERSION: 1.for local branching
#          2.for iterative experiments: added reset function
import copy
import math
import time
import numpy as np
from gurobipy import *


class rflp:
    # Attributes
    master_model = Model()
    sub_model = Model()
    sub_dual = Model()
    sub_cover = Model()  # COVER
    value_y = []
    value_omega = 0
    max_k = 0
    max_Lk = []
    constr_y = LinExpr()
    integer_cut = LinExpr()
    UB = float('inf')
    LB = -float('inf')
    y = []
    omega = []
    L=[]
    value_L=[]
    soft = []
    Qk = []
    beta = 0
    gamma = 0
    delta = 0
    epsilon = 0
    lamda = 0
    mu = 0
    nu = 0
    add_cut_scen = []
    iteration = 0
    gap = float('inf')
    int_gap = float('inf')
    dual = 1
    intSP = 0
    error = 0
    Benders_cut = 0
    callback_time = 0
    callback_num = 0
    Benders_cut = 0
    Integer_cut = 0
    warm = 0
    convergence = []
    # cuts adding
    violation = []
    freq = []
    viol_freq = []
    lift = 0
    zero_half = 0
    # Result output
    num_cut = 0  # num of lazycut added
    opt = 0
    # Input
    p = 0
    ni = 0
    nk = 0
    a1 = 0
    a2 = 0
    cd = []
    cdk = []
    sk = []
    # COVER
    ne = 0
    cd_cover = 0
    bije = 0
    # Local Branching
    LB_terminate = 0
    LB_branch = 0
    LB_cuts = []
    LB_cuts_y = []
    LB_root = 0 #
    LB_value_y = []
    LB_omega = []
    tabu = []
    bestbound = 0
    vn_end = 0
    pr_end = 0
    # Avoid solving identical relaxed dual-sp and sp
    save_max_Lk_DualLP = [] # for benders cut
    save_max_Lk_SP = [] # for integer cut
    save_y = []
    save_y_int = []

    def __init__(self, p, ni, nk, a1, a2, cd, cdk, sk):
        self.reset()
        self.p = p
        self.ni = ni
        self.nk = nk
        self.a1 = a1
        self.a2 = a2
        self.cd = cd
        self.cdk = cdk
        self.sk = sk
        self.value_y = [0 for i in range(ni)]
        self.violation = [0 for k in range(nk)]
        self.freq = [0 for k in range(nk)]
        self.viol_freq = [0 for k in range(nk)]

    def reset(self):
        # Attributes
        self.master_model = Model()
        self.sub_model = Model()
        self.sub_dual = Model()
        self.sub_cover = Model()  # COVER
        self.value_y = []
        self.value_omega = 0
        self.max_k = 0
        self.max_Lk = []
        self.constr_y = LinExpr()
        self.integer_cut = LinExpr()
        self.UB = float('inf')
        self.LB = -float('inf')
        self.y = []
        self.L=[]
        self.value_L=[]
        self.omega = []
        self.soft = []
        self.Qk = []
        self.beta = 0
        self.gamma = 0
        self.delta = 0
        self.epsilon = 0
        self.lamda = 0
        self.mu = 0
        self.nu = 0
        self.add_cut_scen = []
        self.iteration = 0
        self.gap = float('inf')
        self.int_gap = float('inf')
        self.dual = 1
        self.intSP = 0
        self.error = 0
        self.Benders_cut = 0
        self.callback_time = 0
        self.callback_num = 0
        self.Benders_cut_call = 0
        self.Integer_cut = 0
        self.warm = 0
        self.convergence = []
        # cuts adding
        self.violation = []
        self.freq = []
        self.viol_freq = []
        self.lift = 0
        self.zero_half = 0
        # Result output
        self.num_cut = 0  # num of lazycut added
        self.opt = 0
        # Input
        self.p = 0
        self.ni = 0
        self.nk = 0
        self.a1 = 0
        self.a2 = 0
        self.cd = []
        self.cdk = []
        self.sk = []
        # COVER
        self.ne = 0
        self.cd_cover = 0
        self.bije = 0
        # Local Branching
        self.LB_terminate = 0
        self.LB_branch = 0
        self.LB_cuts = []
        self.LB_root = 0
        self.LB_value_y = []
        self.LB_omega = []
        self.tabu = []
        self.bestbound = 0
        self.vn_end = 0
        self.pr_end = 0

        self.save_max_Lk_DualLP = [] # for benders cut
        self.save_max_Lk_SP = [] # for integer cut
        self.save_y = []
        self.save_y_int = []


    def master(self, relax=0):
        # Create variables
        # x:allocations y:location L:auxiliary variable
        x = self.master_model.addVars(
            self.ni, self.ni, vtype=GRB.CONTINUOUS, name="x")
        if relax == 0:
            self.y = self.master_model.addVars(
                self.ni, vtype=GRB.BINARY, name="y")
        else:
            self.y = self.master_model.addVars(
                self.ni, vtype=GRB.CONTINUOUS, name="y")
        self.soft=self.master_model.addVar(
            vtype=GRB.BINARY, name = "soft")
        self.L = self.master_model.addVar(
            vtype=GRB.CONTINUOUS, obj=self.a1, name="L")
        self.omega = self.master_model.addVar(lb=0, ub=float(
            'inf'), vtype=GRB.CONTINUOUS, obj=self.a2, name="omega")

        # Set objective to minimize
        self.master_model.modelSense = GRB.MINIMIZE
        # (1) Maximum cost constraints (objective): L>sum(cdx) forall i
        cdx = x.copy()
        for i in range(self.ni):
            for j in range(self.ni):
                cdx[i, j] = self.cd[i][j]
        self.master_model.addConstrs(
            (x.prod(cdx, i, '*') <= self.L for i in range(self.ni)),
            "epigraph")
        # (2) Constraints sum(y)=p
        self.master_model.addConstr(
            (self.y.sum() == self.p),
            "p")
        # (3) x<=y forall i,j
        self.master_model.addConstrs(
            (x[i, j] <= self.y[j] for i in range(self.ni)
             for j in range(self.ni)),
            "x<y")
        # (4) sum(x)=1 forall i
        self.master_model.addConstrs(
            (x.sum(i, '*') == 1 for i in range(self.ni)),
            "sumx")
        self.master_model.update()  # build master problem model
        # return self.master_model

    def sub(self, callback=0):
        # ---------- Sub problem ----------
        if callback == 0:
            self.update_y()
        # v:allocations u:location L3,eta: auxiliary variable
        v = self.sub_model.addVars(
            self.nk, self.ni, self.ni, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="v")
        if self.intSP == 0:
            u = self.sub_model.addVars(
                self.nk, self.ni, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="u")
        elif self.intSP == 1:
            u = self.sub_model.addVars(
                self.nk, self.ni, vtype=GRB.BINARY, name="u")
        L3 = self.sub_model.addVars(self.nk, vtype=GRB.CONTINUOUS, name="L3")
        eta = self.sub_model.addVar(vtype=GRB.CONTINUOUS, obj=1, name="eta")
        self.sub_model.modelSense = GRB.MINIMIZE
        # (5) eta == sum(L3(k)) forall k
        self.sub_model.addConstr(
            (eta == quicksum(L3)),
            "eta=sumL")
        # (6) L3(k) >= c'd'v(k) forall i,k
        cdv = v.copy()
        for k in range(self.nk):
            for i in range(self.ni):
                for j in range(self.ni):
                    cdv[k, i, j] = self.cdk[k][i][j]
        self.sub_model.addConstrs(
            (v.prod(cdv, k, i, '*') <= L3[k]
             for k in range(self.nk) for i in range(self.ni)),
            "beta")
        # (7) v(k) <= y + u(k) forall k,i,j
        self.sub_model.addConstrs(
            (v[k, i, j] <= self.value_y[j] + u[k, j] for k in range(self.nk)
             for i in range(self.ni) for j in range(self.ni)),
            "gamma")
        # (8) v(k) <= 1 - a_j(k) forall k,i,j
        self.sub_model.addConstrs(
            (v[k, i, j] <= 1 - self.sk[k][j] for k in range(self.nk)
             for i in range(self.ni) for j in range(self.ni)),
            "delta")
        # (9) sum(v) = 1 forall i,k
        self.sub_model.addConstrs(
            (v.sum(k, i, '*') == 1 for k in range(self.nk)
             for i in range(self.ni)),
            "epsilon")
        # (10) u(k) + y <= 1 forall k,j
        self.sub_model.addConstrs(
            (u[k, j] + self.value_y[j] <= 1 for k in range(self.nk)
             for j in range(self.ni)),
            "lamda")
        # (11) u(k) + a_j(k) <= 1 forall k,j
        self.sub_model.addConstrs(
            (u[k, j] + self.sk[k][j] <= 1 for k in range(self.nk)
             for j in range(self.ni)),
            "mu")
        # (12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
        ky = [0 for k in range(self.nk)]
        for k in range(self.nk):
            ky[k] = sum([self.sk[k][j] * self.value_y[j]
                         for j in range(self.ni)])
        self.sub_model.addConstrs(
            (u.sum(k, '*') + sum(self.value_y) -
             ky[k] == self.p for k in range(self.nk)),
            "nu")
        # build subproblem model(get dual variables through Constr.getAttr(Pi))
        self.sub_model.update()
        # return(self.sub_model)

    def dual_sub(self, callback=0):
        if callback == 0:
            self.update_y()
        # Create dual sub model
        self.sub_dual = Model('dual sub model')
        # beta gamma delta epsilon lamda mu nu
        self.beta = self.sub_dual.addVars(
            self.nk, self.ni, ub=0, lb=-float('inf'), vtype=GRB.CONTINUOUS, name="beta")
        self.gamma = self.sub_dual.addVars(
            self.nk, self.ni, self.ni, ub=0, lb=-float('inf'), vtype=GRB.CONTINUOUS, name="gamma")
        self.delta = self.sub_dual.addVars(
            self.nk, self.ni, self.ni, ub=0, lb=-float('inf'), vtype=GRB.CONTINUOUS, name="delta")
        self.epsilon = self.sub_dual.addVars(
            self.nk, self.ni, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="epsilon")
        self.lamda = self.sub_dual.addVars(
            self.nk, self.ni, ub=0, lb=-float('inf'), vtype=GRB.CONTINUOUS, name="lamda")
        self.mu = self.sub_dual.addVars(
            self.nk, self.ni, ub=0, lb=-float('inf'), vtype=GRB.CONTINUOUS, name="mu")
        self.nu = self.sub_dual.addVars(
            self.nk, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="nu")
        # Qk are auxiliary variables, maximize every subproblem
        self.Qk = self.sub_dual.addVars(
            self.nk, vtype=GRB.CONTINUOUS, obj=1, name="Qk")
        # Set sub model objective to maximize
        self.sub_dual.modelSense = GRB.MAXIMIZE
        # (1)
        self.sub_dual_obj()
        # (2) -sum_i(gamma_kij)+lambda_kj+mu_kj+nu_k<=0  forall k,j
        self.sub_dual.addConstrs(
            (-self.gamma.sum(k, '*', j) + self.lamda[k, j] + self.mu[k, j] + self.nu[k]
             <= 0 for k in range(self.nk) for j in range(self.ni)),
            "u")
        # (3) c_kij*d_ki*beta_ki+gamma_kij+delta_kij+epsilon_ki<=0  forall k,i,j
        self.sub_dual.addConstrs(
            (self.cdk[k][i][j] * self.beta[k, i] + self.gamma[k, i, j] + self.delta[k, i, j] + self.epsilon[k, i] <= 0
             for k in range(self.nk) for i in range(self.ni) for j in range(self.ni)),
            "v")
        # (4) -sum_i(beta_i)<=1 forall k
        self.sub_dual.addConstrs(
            (-self.beta.sum(k, '*') <= 1 for k in range(self.nk)),
            "L3")
        self.sub_dual.update()  # build dual subproblem model

    def cover_pre(self):  # COVER
        # preprocessing cdk
        cd_cover = [[] for k in range(self.nk)]
        self.ne = [0 for k in range(self.nk)]
        for k in range(self.nk):
            cd0 = list(itertools.chain.from_iterable(
                self.cdk[k]))  # combine lists
            cd_cover[k] = sorted(set(cd0))  # sort without duplicates
            self.ne[k] = len(cd_cover[k])
        # Build all b_ije(k)
        b = [[[[0 for e in range(self.ne[k])] for j in range(self.ni)]
              for i in range(self.ni)] for k in range(self.nk)]
        for k in range(self.nk):  # !four layers
            for i in range(self.ni):
                for j in range(self.ni):
                    for e in range(self.ne[k]):
                        if self.cdk[k][i][j] <= cd_cover[k][e]:
                            b[k][i][j][e] = 1
        self.cd_cover = cd_cover #
        self.bije = b

    def cover_sub(self):  # COVER
        # self.sub_dual = Model("p-center-cover")
        # Create variables
        # z:ordered cost, y:location
        z = [[0 for e in range(self.ne[k])] for k in range(self.nk)]
        for k in range(self.nk):
            for e in range(self.ne[k]):
                z_name = ''.join(['z[', str(k), ',', str(e), ']'])
                z[k][e] = self.sub_cover.addVar(vtype=GRB.BINARY, name=z_name)
        y = self.sub_cover.addVars(self.nk, self.ni, vtype=GRB.BINARY, name="y")
        L = self.sub_cover.addVars(
            self.nk, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="L")  # ??
        # Set objective to minimize
        self.sub_cover.modelSense = GRB.MINIMIZE
        # Minimize: \sum L
        self.sub_cover.setObjective(quicksum(L))
        # (00) Q = max(L(k)) # ???Will obj affect solution speed: Obj:max Q & max sum(L(k))
        # (0) L(k) = \sum_e cd_cover[k][e]z[k][e]  \forall k
        for k in range(self.nk):
            sum_cdz = 0
            for e in range(self.ne[k]):
                sum_cdz += self.cd_cover[k][e] * z[k][e]
            self.sub_cover.addConstr(
                (L[k] == sum_cdz),
                "L(k)=rho*z")
        # (1) \sum_j b_ije[k]*y_j >= z_e \forall nk,i,e
        for k in range(self.nk): # !
            for i in range(self.ni):
                for e in range(self.ne[k]):
                    sum_ay = 0
                    for j in range(self.ni):
                        sum_ay += self.bije[k][i][j][e] * y[k,j]
                    self.sub_cover.addConstr(
                        (sum_ay >= z[k][e]),
                        'ay>e' + str(k) + str(i) + str(e))
        # (2) \sum y_j[k] = p \forall k
        self.sub_cover.update() # ?
        self.sub_cover.addConstrs((y.sum(k,'*') == self.p for k in range(self.nk)),'sump')
        # (3) \sum_e z_e[k] = 1 \forall k
        self.sub_cover.addConstrs((quicksum(z[k]) == 1 for k in range(self.nk)),'sumz')
        self.sub_cover.update() #
        # fix varibale: sk
        # for k in range(self.nk):
        #     for j in range(self.ni):
        #         if self.sk[k][j] == 1:
        #             y_name = ''.join(['y[',str(k),',', str(j), ']'])
        #             self.sub_cover.getVarByName(y_name).lb = 0
        #             self.sub_cover.getVarByName(y_name).ub = 0

    def cover_bound(self): # COVER
        # Function for a simple p-center problem
        def p_center(cd,p=1,LP=0):
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
            m.params.OutputFlag = 0
            m.optimize()
            return m.objVal,m.Runtime
        # For each iteration,
        # Input: value_y (variable) & sk(fixed)
        # Find UB1
        # Exclude unavailable nodes and Find 'farthest' nodes (for each k)
        y_set = [set() for k in range(self.nk)] # Store open facilies
        y_now = [0 for k in range(self.nk)]
        cd_matrix = [[] for k in range(self.nk)]
        cd_matrix_1 = [[] for k in range(self.nk)]
        x = [[[0 for j in range(self.ni)] for i in range(self.ni)] for k in range(self.nk)]
        # t1 = time.time()
        for k in range(self.nk):
            cd_matrix[k] = np.array(self.cdk[k])
            cd_matrix_1[k] = np.array(self.cdk[k])
            for j in range(self.ni):
                y_name = ''.join(['y[',str(k),',', str(j), ']']) # Find and fix y
                if self.value_y[j]-self.sk[k][j] > 0:
                    self.sub_cover.getVarByName(y_name).lb = 1
                    self.sub_cover.getVarByName(y_name).ub = 1
                    y_set[k].update({j})
                elif self.sk[k][j]==1:
                    self.sub_cover.getVarByName(y_name).lb = 0
                    self.sub_cover.getVarByName(y_name).ub = 0
                else:
                    self.sub_cover.getVarByName(y_name).lb = 0
                    self.sub_cover.getVarByName(y_name).ub = 1
                if self.sk[k][j] == 1:
                    cd_matrix[k][:,j] = 0
            while len(y_set[k]) < self.p:
                y_curr = list(y_set[k])
                cd_temp = cd_matrix[k][y_curr,:] #
                (a,b) = np.unravel_index(cd_temp.argmax(), cd_temp.shape)
                cd_matrix[k][y_curr[a],b] = 0
                y_set[k].update({b})
            y_now[k] = [0 for j in range(self.ni)]
            for x in y_set[k]:
                y_now[k][x] = 1
            # Find closest facility for each demand
            for j in range(self.ni):
                if y_now[k][j] == 0:
                    cd_matrix_1[k][:,j] = 1e8 # set j column to Big M
            facility = np.argmin(cd_matrix_1[k], axis=1) # index of maximum value in each row
            cost = cd_matrix_1[k][np.arange(cd_matrix_1[k].shape[0]),facility] # advanced index returning maximum value of each row
            UB1 = max(cost) # Upper Bound UB1
            # Find UB2
            cluster = [[] for n in range(self.p)]
            y_idx = list(y_set[k])
            for i in range(self.p):
                cluster[i] = np.argwhere(facility == y_idx[i]).ravel()#.tolist()
            cd_array = np.array(self.cdk[k])
            L_cluster = [0 for i in range(self.p)]
            for i in range(self.p):
                L_cluster[i],_ = p_center(cd_array[cluster[i][:,None],cluster[i]]) #
            UB2 = max(L_cluster)
            # LB
            LB2,_ = p_center(self.cdk[k],self.p,1)
            # print('Bounds:',UB1)
            # print(UB2)
            # print(LB2)
            # Use UB & LB to restrict 'z' in self.sub_cover; Input self.cd_cover
            # z = [[0 for e in range(self.ne[k])] for k in range(self.nk)]
            for e in range(self.ne[k]):
                z_name = ''.join(['z[', str(k), ',', str(e), ']'])
                if self.cd_cover[k][e] < LB2 and self.cd_cover[k][e] > UB2:
                    # self.sub_cover.getVarByName(z_name).lb = 0
                    self.sub_cover.getVarByName(z_name).ub = 0
                else:
                    # self.sub_cover.getVarByName(z_name).lb = 0
                    self.sub_cover.getVarByName(z_name).ub = 1
        # t2= time.time()-t1
        # print('bound time')
        #print(round(t2,2))

    def sub_dual_obj(self):  # Caculate objective function
        def c_constr1():
            c_delta = [[0 for i in range(self.ni * self.ni)]
                       for k in range(self.nk)]
            for k in range(self.nk):
                c_delta[k] = [(1 - self.sk[k][j])
                              for j in range(self.ni)] * self.ni
            c_lamda = [(1 - self.value_y[j]) for j in range(self.ni)]
            c_mu = [[0 for i in range(self.ni)] for k in range(self.nk)]
            for k in range(self.nk):
                c_mu[k] = [1 - self.sk[k][j] for j in range(self.ni)]
            c_nu = []
            for k in range(self.nk):
                c_nu.append(self.p + sum([self.sk[k][j] * self.value_y[j]
                                          for j in range(self.ni)]) - sum(self.value_y))
            return c_delta, c_lamda, c_mu, c_nu
        # (1) Q(k) = sum_i(sum_j(y_j*gamma_kij))+sumsum_ij((1-a_kj)*delta_kij
        #           + sum_i(epsilon_ki) + sum_j((1-y_j)*lambda_kj) +
        #           + sum_j((1-a_kj)*mu_kj) + (p+sum_j(a_kj*y_j)-sum(y_j))*nu_k   forall k
        c_delta, c_lamda, c_mu, c_nu = c_constr1()  # update coeff
        self.sub_dual.addConstrs(
            (self.Qk[k] == LinExpr(self.value_y * self.ni, self.gamma.select(k, '*', '*')) +
             LinExpr(c_delta[k], self.delta.select(k, '*', '*')) +
             self.epsilon.sum(k, '*') + LinExpr(c_lamda, self.lamda.select(k, '*')) +
             LinExpr(c_mu[k], self.mu.select(k, '*')) + c_nu[k] * self.nu[k] for k in range(self.nk)),
            "Q(k)")

    def update_master(self):
        self.update_cut()
        self.master_model.addConstr(self.omega >= self.constr_y)
        # self.master_model.update()

    def update_sub(self, callback=0):
        if callback == 0:
            self.update_y()
        for k in range(self.nk):
            for i in range(self.ni):
                for j in range(self.ni):
                    gamma_name = ''.join(
                        ['gamma[', str(k), ',', str(i), ',', str(j), ']'])
                    self.sub_model.getConstrByName(
                        gamma_name).rhs = self.value_y[j]
        for k in range(self.nk):
            for j in range(self.ni):
                lamda_name = ''.join(['lamda[', str(k), ',', str(j), ']'])
                self.sub_model.getConstrByName(
                    lamda_name).rhs = 1 - self.value_y[j]
        ky = [0 for k in range(self.nk)]
        for k in range(self.nk):
            ky[k] = sum([self.sk[k][j] * self.value_y[j]
                         for j in range(self.ni)])
        for k in range(self.nk):
            nu_name = ''.join(['nu[', str(k), ']'])
            self.sub_model.getConstrByName(
                nu_name).rhs = self.p + ky[k] - sum(self.value_y)
        self.sub_model.update()
        # return self.sub_model

    def update_sub_dual(self, callback=0):
        if callback == 0:
            self.update_y()
        for k in range(self.nk):
            constr_name = ''.join(['Q(k)[', str(k), ']'])
            self.sub_dual.remove(self.sub_dual.getConstrByName(constr_name))
        self.sub_dual_obj()
        self.sub_dual.update()

    def update_y(self):
        self.value_y = []
        for j in range(self.ni):
            y_name = ''.join(['y[', str(j), ']'])
            y_temp = self.master_model.getVarByName(y_name)
            self.value_y.append(y_temp.x)
        # update value of y for subproblem in each iteration
        self.value_y = [round(x) for x in self.value_y]

    def update_cut(self, numk=None, lift=0,save = 0):
        gamma = []
        delta = []
        epsilon = []
        lamda = []
        mu = []
        beta = []
        nu = []
        if numk == None:
            numk = self.max_k
        if save == 0:
            gamma = [[0 for j in range(self.ni)] for i in range(self.ni)]
            delta = [[0 for j in range(self.ni)] for i in range(self.ni)]
            epsilon = [0 for j in range(self.ni)]
            lamda = [0 for j in range(self.ni)]
            mu = [0 for j in range(self.ni)]
            beta = [0 for j in range(self.ni)]
            if self.dual == 0:  # for primal sub problem
                cm1 = self.sub_model.getConstrs()
                num_c = len(cm1)
                dual_value = []
                constrname = []
                for i in range(num_c):
                    dual_value.append(cm1[i].getAttr('Pi'))
                    constrname.append(cm1[i].getAttr('ConstrName'))
                dual = dict(zip(constrname, dual_value))
                for i in range(self.ni):
                    for j in range(self.ni):
                        gamma_name = ''.join(
                            ['gamma[', str(numk), ',', str(i), ',', str(j), ']'])
                        delta_name = ''.join(
                            ['delta[', str(numk), ',', str(i), ',', str(j), ']'])
                        gamma[i][j] = dual[gamma_name]
                        delta[i][j] = dual[delta_name]
                for n in range(self.ni):
                    epsilon_name = ''.join(
                        ['epsilon[', str(numk), ',', str(n), ']'])
                    lamda_name = ''.join(['lamda[', str(numk), ',', str(n), ']'])
                    mu_name = ''.join(['mu[', str(numk), ',', str(n), ']'])
                    epsilon[n] = dual[epsilon_name]
                    lamda[n] = dual[lamda_name]
                    mu[n] = dual[mu_name]
                nu_name = ''.join(['nu[', str(numk), ']'])
                nu = dual[nu_name]
            elif self.dual == 1:  # for dual sub problem
                for i in range(self.ni):
                    for j in range(self.ni):
                        gamma_name = ''.join(
                            ['gamma[', str(numk), ',', str(i), ',', str(j), ']'])
                        delta_name = ''.join(
                            ['delta[', str(numk), ',', str(i), ',', str(j), ']'])
                        gamma[i][j] = self.sub_dual.getVarByName(gamma_name).x
                        delta[i][j] = self.sub_dual.getVarByName(delta_name).x
                for n in range(self.ni):
                    beta_name = ''.join(['beta[', str(numk), ',', str(n), ']'])
                    epsilon_name = ''.join(
                        ['epsilon[', str(numk), ',', str(n), ']'])
                    lamda_name = ''.join(['lamda[', str(numk), ',', str(n), ']'])
                    mu_name = ''.join(['mu[', str(numk), ',', str(n), ']'])
                    beta[n] = self.sub_dual.getVarByName(beta_name).x
                    epsilon[n] = self.sub_dual.getVarByName(epsilon_name).x
                    lamda[n] = self.sub_dual.getVarByName(lamda_name).x
                    mu[n] = self.sub_dual.getVarByName(mu_name).x
                nu_name = ''.join(['nu[', str(numk), ']'])
                nu = self.sub_dual.getVarByName(nu_name).x
        else:
            gamma = save[0][numk]
            delta = save[1][numk]
            epsilon = save[2][numk]
            lamda = save[3][numk]
            mu = save[4][numk]
            beta = save[5][numk]
            nu = save[6][numk]
        # dual_lifting
        # if lift == 1:
            # sum_gamma = [0 for i in range(self.ni)]
            # for j in range(self.ni):
            #     for i in range(self.ni):
            #         sum_gamma[j] += gamma[i][j]
            # for j in range(self.ni):
            #     if self.value_y[j] == 1:
            #         new_lamda = -nu-mu[j]+sum_gamma[j]
            #         if lamda[j] < new_lamda:
            #            # print('lift lamda')
            #             lamda[j] = new_lamda
            # for i in range(self.ni):
            #     for j in range(self.ni):
            #         if self.value_y[j] == 0:
            #             new_gamma = -epsilon[i] - delta[i][j] - \
            #                 self.cdk[numk][i][j] * beta[i]
            #             if gamma[i][j] < new_gamma:
            #                 # print('lift gamma')
            #                 gamma[i][j] = new_gamma
        # Benders' cut
        # omega >= sumsum(gamma_k'ij*y) + sum_j(-lamda*y) +  nu*sum_j((aj(k')-1)*y)
        # + sumsum((1-aj(k'))*delta_ij)+sum_i(epsilon)+sum_j(lamda)
        # + sum_j(1-aj(k'))*mu_j + p*nu
        gamma_y = []
        for i in range(self.ni):
            for j in range(self.ni):
                gamma_y.append(gamma[i][j])
        ajk_y = []
        for j in range(self.ni):
            ajk_y.append(self.sk[numk][j] - 1)
        c_y = LinExpr(gamma_y, self.y.select() * self.ni) - LinExpr(lamda, self.y.select()) \
            + nu * LinExpr(ajk_y, self.y.select())
        constant_delta = 0
        for i in range(self.ni):
            constant_delta += quicksum([(1 - self.sk[numk][j])
                                        * delta[i][j] for j in range(self.ni)])
        constant = quicksum(epsilon) + quicksum(lamda) + constant_delta +\
            quicksum([(1 - self.sk[numk][j]) * mu[j]
                      for j in range(self.ni)]) + self.p * nu
        # update cut to be added to master problem in each iteration
        self.constr_y = c_y + constant

    def get_subdual_vals(self):
        gamma = [[[0 for j in range(self.ni)] for i in range(self.ni)] for k in range(self.nk)]
        delta = [[[0 for j in range(self.ni)] for i in range(self.ni)] for k in range(self.nk)]
        epsilon = [[0 for j in range(self.ni)] for k in range(self.nk)]
        lamda = [[0 for j in range(self.ni)] for k in range(self.nk)]
        mu =[[0 for j in range(self.ni)] for k in range(self.nk)]
        beta = [[0 for j in range(self.ni)] for k in range(self.nk)]
        nu = [0 for k in range(self.nk)]
        for numk in range(self.nk):
            for i in range(self.ni):
                for j in range(self.ni):
                    gamma_name = ''.join(
                        ['gamma[', str(numk), ',', str(i), ',', str(j), ']'])
                    delta_name = ''.join(
                        ['delta[', str(numk), ',', str(i), ',', str(j), ']'])
                    gamma[numk][i][j] = self.sub_dual.getVarByName(gamma_name).x
                    delta[numk][i][j] = self.sub_dual.getVarByName(delta_name).x
            for n in range(self.ni):
                beta_name = ''.join(['beta[', str(numk), ',', str(n), ']'])
                epsilon_name = ''.join(
                    ['epsilon[', str(numk), ',', str(n), ']'])
                lamda_name = ''.join(['lamda[', str(numk), ',', str(n), ']'])
                mu_name = ''.join(['mu[', str(numk), ',', str(n), ']'])
                beta[numk][n] = self.sub_dual.getVarByName(beta_name).x
                epsilon[numk][n] = self.sub_dual.getVarByName(epsilon_name).x
                lamda[numk][n] = self.sub_dual.getVarByName(lamda_name).x
                mu[numk][n] = self.sub_dual.getVarByName(mu_name).x
            nu_name = ''.join(['nu[', str(numk), ']'])
            nu[numk] = self.sub_dual.getVarByName(nu_name).x
        save_subdual_vals = [gamma,delta,epsilon,lamda,mu,beta,nu]
        return save_subdual_vals

    def update_integer_cut(self,cover=0,save=0):
        # simplified cut (use no dual information)
        sum_c_y = LinExpr()
        for i in range(self.ni):
            if self.value_y[i] == 1:
                sum_c_y += self.y[i] - 1
            elif self.value_y[i] == 0:
                sum_c_y += -self.y[i]
        sum_c_y += 1
        if save!=0 and cover==0:
            self.max_Lk = save
        elif cover == 0 and save == 0:
            self.max_Lk = self.worst_scenario(1)
        else:
            self.max_Lk = self.worst_scenario(1,1)
        self.integer_cut = self.max_Lk[0] * sum_c_y
        self.LB_cuts.append(self.integer_cut)
        # print('integer cut')

    def gap_calculation(self, MIP_SP=0, Check_optimal=0,cover=0):
        if MIP_SP == 1:  # in mycallback
            # vals = self.master_model.cbGetSolution(self.master_model._vars)
            # value_L = vals(-2)
            if cover == 0:
                self.max_Lk = self.worst_scenario(1)
            else:
                self.max_Lk = self.worst_scenario(1,1)
            self.int_gap = self.max_Lk[0] - self.value_omega #
           # print('self.max_Lk:',self.max_Lk[0])
        else:
            # extract L
            var_L = self.master_model.getVarByName('L')
            value_L = var_L.x
            # update LB = objective value of master problem
            obj_master = self.master_model.getObjective()
            self.LB = obj_master.getValue()
            if Check_optimal == 0:
                self.max_Lk = self.worst_scenario()
            elif Check_optimal == 1:
                self.max_Lk = self.worst_scenario(1)
            # update UB
            self.UB = min([self.UB, self.a1 * value_L + self.a2 * self.max_Lk[0]])
            self.gap = (self.UB - self.LB) / self.LB
        # return self.gap

    def worst_scenario(self, MIP_SP=0, cover=0):
        if cover == 0:
            if self.dual == 0 or MIP_SP == 1:
                value_L3 = []
                for k in range(self.nk):
                    L3_name = ''.join(['L3[', str(k), ']'])
                    L3_temp = self.sub_model.getVarByName(L3_name)
                    value_L3.append(L3_temp.x)
                # maximum L3 and its index (worst k)
                self.max_Lk = max([[v, i] for i, v in enumerate(value_L3)])
            elif self.dual == 1 and MIP_SP != 1:
                value_Qk = []
                for k in range(self.nk):
                    Qk_name = ''.join(['Qk[', str(k), ']'])
                    Qk_temp = self.sub_dual.getVarByName(Qk_name)
                    value_Qk.append(Qk_temp.x)
                # maximum Qk and its index (worst k)
                self.max_Lk = max([[v, i] for i, v in enumerate(value_Qk)])
            self.max_k = self.max_Lk[1]
        elif cover == 1:
            value_L = []
            for k in range(self.nk):
                L_name = ''.join(['L[', str(k), ']'])
                value_L.append(self.sub_cover.getVarByName(L_name).x)
            self.max_Lk = max([[v, i] for i, v in enumerate(value_L)])
            self.max_k = self.max_Lk[1]
        return self.max_Lk

    def update_status(self):
        self.add_cut_scen.append(self.max_k)
        self.iteration += 1
        # print('==========================================')
        # print('Current iteration:', str(self.iteration))
        # print('gap = ', str(self.gap))
        # print('int_gap = ',str(self.int_gap))
        # print('Cuts added form scenario:', str(self.add_cut_scen[0:-1]))

    def params_tuneup(self, accu=0):
        # References:
        #m1.params.ScaleFlag = 3
        #m1.params.ObjScale = 100
        #m1.params.NumericFocus = 3
        #m1.params.NormAdjust = 3
        # m1.params.InfUnbdInfo = 0 #1
        # m1.params.Quad = -1 #1
        # m1.params.Sifting = -1 #2
        # m1.params.SiftMethod = -1 # 2
        # m1.params.SimplexPricing = -1 #3
        #m1.params.Method = -1
        #m1.params.AggFill = 0
        #m1.params.Aggregate = 0
        # m1.params.DualReductions = 1 #0
        # m1.params.PreDual = 2 #2
        #m1.params.Presolve = 0
        # tune parameters to avoid numerical issues for subproblem
        # wrong optimal solutions appear for both sub&dual_sub
        self.master_model.params.OutputFlag = 1
        self.sub_model.params.OutputFlag = 0
        self.sub_dual.params.OutputFlag = 0
        self.sub_cover.params.OutputFlag = 0
        # self.master_model.params.Cuts = 0
        # self.master_model.params.Presolve = 0
        # self.master_model.params.Heuristics = 0
        if accu == 1:
            # self.master_model.params.Cuts = 0
            # self.sub_model.params.Cuts = 0
            # self.sub_dual.params.Cuts = 0

            self.master_model.params.Presolve = 0
            self.master_model.params.ScaleFlag = 0
            self.master_model.params.NumericFocus = 0

            self.sub_model.params.Presolve = 0
            self.sub_model.params.ScaleFlag = 0
            self.sub_model.params.NumericFocus = 0

            self.sub_dual.params.Presolve = 0
            self.sub_dual.params.ScaleFlag = 0
            self.sub_dual.params.NumericFocus = 0

    def error_check(self):
        if self.gap <= -0.1:
            self.error = 1
            print('# WARNING: gap is a negative value', '\n',
                  '   Wrong dual problem solution.')
        return self.error

    def update_multiple_scenario(self,SP_Qk,save=0):
        # print('===============sorting==================')
        # extract omega and Q(k)
        # Verify cut improvement #LB
        violation_now = [
            x - self.value_omega for x in SP_Qk]

        # for k in range(self.nk):
        #     if violation_now[k] > 0:
        #         self.violation[k] += violation_now[k]
        #         self.freq[k] += 1
        #         self.viol_freq[k] = self.violation[k] / self.freq[k]
        # print(self.viol_freq)
        # rank = sorted(range(len(self.viol_freq)), reverse=True, key=self.viol_freq.__getitem__)
        rank = sorted(range(len(violation_now)), reverse=True,
                      key=violation_now.__getitem__)
        # print(rank[0])
        # odd_1 = 0
        # odd_2 = 0
        # for n in range(round(self.nk / 4 + 1)): #!!!fewer cuts
        for n in range(round(self.nk)):
            # for n in range(1):
            if violation_now[rank[n]] > 0:
                if save == 0:
                    self.update_cut(rank[n], self.lift)
                else:
                    self.update_cut(rank[n], self.lift,save)
                if self.zero_half == 0:
                    self.master_model.cbLazy(self.omega >= self.constr_y)
                    self.LB_cuts.append(self.constr_y)
                    # self.LB_cuts_y.append(self.value_y)
                    self.Benders_cut += 1
                # elif self.zero_half == 1:
                #     if n == 0:
                #         self.master_model.cbLazy(self.omega >= self.constr_y)
                #         coeff_strong, constant_strong = self.reformulate_linexpr(
                #             self.constr_y)
                #         constr_strong = LinExpr(
                #             coeff_strong, self.y.select()) + constant_strong
                #         if [k for k in coeff_strong if k % 2] != []:
                #             odd_1 = 1
                #         # print(constr_strong)
                #     else:
                #         # zero-half
                #         coeff_constr_y, constant_y = self.reformulate_linexpr(
                #             self.constr_y)
                #         # if self.master_model.cbGet(GRB.Callback.MIPSOL_NODCNT) == 0:
                #         if [k for k in coeff_constr_y if k % 2] != []:
                #             odd_2 = 1
                #         if odd_1 == 1 and odd_2 == 1:
                #             coeff_sum, constant_sum = self.reformulate_linexpr(
                #                 self.constr_y + constr_strong)
                #             coeff_sum = [math.ceil(0.5 * x) for x in coeff_sum]
                #             constant_sum = math.ceil(0.5 * constant_sum)
                #             constr_sum = LinExpr(
                #                 coeff_sum, self.y.select()) + constant_sum
                #             self.master_model.cbLazy(self.omega >= constr_sum)
                #             # print(const r_sum)
                #         else:
                #             # print(self.constr_y)
                #             self.master_model.cbLazy(
                #                 self.omega >= self.constr_y)
            else:
                break  # Callback: update scenario violation and frequency⁠

    def reformulate_linexpr(self, formulation):
        coeff_y = [0 for j in range(self.ni)]
        for n in range(formulation.size()):
            y_index = int(formulation.getVar(n).VarName[2])
            coeff_y[y_index] += formulation.getCoeff(n)
        constant_y = formulation.getConstant()
        return coeff_y, constant_y

    def warm_start(self, warm=0):
        # self.master_model.optimize()
        # self.update_sub_dual(0)
        # self.sub_dual.optimize()
        # self.gap_calculation()
        # self.master_model.addConstr(
        #     self.a1 * self.master_model.getVars()[-1] + self.a1 * self.omega >= self.LB)
        if warm == 0:
            # chase the carrot??
            self.master(1)
            self.params_tuneup()
            self.master_model.optimize()
            self.update_y()
            # interior point:(not sure)
            y_in = [self.p / self.ni for j in range(self.ni)]
            y_optimal = self.value_y
            bound_no_impro = 0
            last_bound = GRB.INFINITY
            inter = 0.1  # step length
            y_step = 0.1  # step length
            while bound_no_impro < 5:
                self.value_y = [inter * a +
                                (1 - inter) * b for a, b in zip(y_optimal, y_in)]
                # self.value_y = round_y(self.value_y)
                self.update_sub_dual(1)
                self.sub_dual.optimize()
                self.update_cut()
                self.master_model.addConstr(self.omega >= self.constr_y)
                self.master_model.optimize()
                self.update_y()
                y_optimal = self.value_y
                y_in = [y_step * a + (1 - y_step) * b for a,
                        b in zip(y_in, y_optimal)]
                if last_bound == self.master_model.getObjective().getValue():
                    bound_no_impro += 1
                last_bound = self.master_model.getObjective().getValue()
            # reset y to binary
            for j in range(self.ni):
                y_name = ''.join(['y[', str(j), ']'])
                self.master_model.getVarByName(y_name).vtype = 'B'
                # print(self.master_model.getVarByName(y_name).vtype)
            self.warm = 'over'
        else:
            self.master()
            self.params_tuneup()
            self.warm = 'over'

    def round_y(self, value_y):  # ?
        new_y = [0 for j in range(self.ni)]
        rank = sorted(range(len(value_y)), reverse=True,
                      key=value_y.__getitem__)
        for j in range(self.p):
            new_y[rank[j]] = 1
        return new_y

    def add_LB(self, Branching_record,neighbourhood,initial=0): # add Hanmming constraints in master_model
        # Input: y(var); value_y
        # Delta(y_new,y_now) = [value_y - y]
        delta_y = 0
        bst_value_y = Branching_record[1][-3 - self.ni:-3]
        for j in range(self.ni):
            if bst_value_y[j] == 1:
                delta_y += 1 - self.y[j]
            else:
                delta_y += self.y[j]
        if initial == 1:
            self.master_model.addConstr(delta_y >= 2)
        self.master_model.addConstr(delta_y <= neighbourhood)
        self.master_model.update() #

    def add_master_bound(self,UB=0,LB=0):
        self.master_model.addConstr(self.a1*self.L+self.a2*self.omega <= UB)

    def set_initial(self,value):
        Vars = self.master_model.getVars()
        for n in range(len(Vars)):
            Vars[n].Start = value[n]
        # for j in range(self.ni):
        #     y_name = ''.join(['y[', str(j), ']'])
        #     self.master_model.getVarByName(y_name).Start = value[j]

    def add_proximity(self,Branching_record,impro = 0,soft = 0):
        bigM = 1e5 # soft
        delta_y = 0
        bst_value_y = Branching_record[1][-3 - self.ni:-3]
        UB = Branching_record[0]
        LB = self.bestbound
        rhs = (UB + LB)/2
        soft_rhs = rhs+(UB-rhs)/2
        for j in range(self.ni):
            if bst_value_y[j] == 1:
                delta_y += 1 - self.y[j]
            else:
                delta_y += self.y[j]
        self.master_model.setObjective(delta_y+self.soft*bigM) # set obj
        self.master_model.addConstr(delta_y>=2)
        self.master_model.addConstr(
            self.a1*self.L+self.a2*self.omega <= rhs+self.soft*((UB-rhs)/2))
        # print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        # print('rhs = ',rhs,' soft rhs= ',rhs+((UB-rhs)/2))
        return rhs,soft_rhs

    def remove_proximity(self):
        self.master_model.setObjective(self.a1*self.L+self.a2*self.omega)

    def record_best_sol(self,Branching_record,start_time):
        best_incumbent = []
        improve = 0
        if self.master_model.Objval < Branching_record[0]:
            improve = 1
            Vars = self.master_model.getVars()
            for n in Vars:
                best_incumbent.append(n.x)
            Branching_record = [self.master_model.Objval,best_incumbent,time.time()-start_time]
        # else:
        #     Branching_record = Branching_record
        return Branching_record,improve
