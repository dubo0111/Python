# class for:
# bulid or update specific model
# get variables,dual variables
from gurobipy import *
class rflp:
    # Attributes
    master_model = Model()
    sub_model = Model()
    value_y = []
    max_k = 0
    constr_y = LinExpr()
    UB = float('inf')
    LB = -float('inf')
    y = []
    omega = []
    add_cut_scen = []
    iteration = 0
    gap = 0
    # Input
    p = 0
    ni = 0
    nk = 0
    a1= 0
    a2 = 0
    cd = []
    cdk = []
    sk = []
    def __init__(self,p,ni,nk,a1,a2,cd,cdk,sk):
        self.p = p
        self.ni = ni
        self.nk = nk
        self.a1 = a1
        self.a2 = a2
        self.cd = cd
        self.cdk = cdk
        self.sk = sk
        self.value_y = [0 for i in range(ni)]
    def master(self):
        # Create variables
        # x:allocations y:location L:auxiliary variable
        x = self.master_model.addVars(self.ni,self.ni,vtype=GRB.BINARY, name="x")
        self.y = self.master_model.addVars(self.ni,vtype=GRB.BINARY, name="y")
        L = self.master_model.addVar(vtype=GRB.CONTINUOUS,obj = self.a1, name="L")
        self.omega = self.master_model.addVar(lb=0,ub=float('inf'),vtype=GRB.CONTINUOUS,obj=self.a2, name="omega")
        # Set objective to minimize
        self.master_model.modelSense = GRB.MINIMIZE
        # (1) Maximum cost constraints (objective): L>sum(cdx) forall i
        cdx = x.copy()
        for i in range(self.ni):
            for j in range(self.ni):
               cdx[i,j]=self.cd[i][j]
        self.master_model.addConstrs(
                (x.prod(cdx,i,'*') <= L for i in range(self.ni)),
                "epigraph")
        # (2) Constraints sum(y)=p
        self.master_model.addConstr(
                (self.y.sum() == self.p),
                "p")
        # (3) x<=y forall i,j
        self.master_model.addConstrs(
                (x[i,j] <= self.y[j] for i in range(self.ni) for j in range(self.ni)),
                "x<y")
        # (4) sum(x)=1 forall i
        self.master_model.addConstrs(
                (x.sum(i,'*') == 1 for i in range(self.ni)),
                "sumx")
        self.master_model.update()
        #return self.master_model

    def sub(self,callback = 0):
        # ---------- Sub problem ----------
        if callback == 0:
            self.update_y()
        # v:allocations u:location L3,eta: auxiliary variable
        v = self.sub_model.addVars(self.nk,self.ni,self.ni,lb=0,ub=1, vtype=GRB.CONTINUOUS, name="v")
        u = self.sub_model.addVars(self.nk,self.ni,lb=0,ub=1,vtype=GRB.CONTINUOUS, name="u")
        L3 = self.sub_model.addVars(self.nk,vtype=GRB.CONTINUOUS, name="L3")
        eta = self.sub_model.addVar(vtype=GRB.CONTINUOUS, obj=1, name="eta")
        self.sub_model.modelSense = GRB.MINIMIZE
        #(5) eta == sum(L3(k)) forall k
        self.sub_model.addConstr(
                (eta == L3.sum()),
                "eta=sumL")
        #(6) L3(k) >= c'd'v(k) forall i,k
        cdv = v.copy()
        for k in range(self.nk):
            for i in range(self.ni):
                for j in range(self.ni):
                   cdv[k,i,j]=self.cdk[k][i][j]
        self.sub_model.addConstrs(
                (v.prod(cdv,k,i,'*') <= L3[k] for k in range(self.nk) for i in range(self.ni)),
                "beta")
        #(7) v(k) <= y + u(k) forall k,i,j
        self.sub_model.addConstrs(
                (v[k,i,j] <= self.value_y[j] + u[k,j] for k in range(self.nk) for i in range(self.ni) for j in range(self.ni)),
                "gamma")
        #(8) v(k) <= 1 - a_j(k) forall k,i,j
        self.sub_model.addConstrs(
                (v[k,i,j] <= 1 - self.sk[k][j] for k in range(self.nk) for i in range(self.ni) for j in range(self.ni)),
                "delta")
        #(9) sum(v) = 1 forall i,k
        self.sub_model.addConstrs(
                (v.sum(k,i,'*') == 1 for k in range(self.nk) for i in range(self.ni)),
                "epsilon")
        #(10) u(k) + y <= 1 forall k,j
        self.sub_model.addConstrs(
                (u[k,j] + self.value_y[j] <= 1 for k in range(self.nk) for j in range(self.ni)),
                "lamda")
        #(11) u(k) + a_j(k) <= 1 forall k,j
        self.sub_model.addConstrs(
                (u[k,j] + self.sk[k][j] <= 1 for k in range(self.nk) for j in range(self.ni)),
                "mu")
        #(12) sum(u(k)) + sum(y) - sum(a_j(k)*y) = p forall k (or <=)
        ky = [0 for k in range(self.nk)]
        for k in range(self.nk):
            ky[k] = sum([self.sk[k][j]*self.value_y[j] for j in range(self.ni)])
        self.sub_model.addConstrs(
                (u.sum(k,'*') + sum(self.value_y) - ky[k] == self.p for k in range(self.nk)),
                "nu")
        self.sub_model.update()
        #return(self.sub_model)

    def update_master(self):
        self.update_cut()
        #self.master_model.getVarByName('omega').Obj = self.a2
        self.master_model.addConstr(self.omega >= self.constr_y)#2054.9917037337914)#
        self.master_model.update()

    def update_sub(self,callback = 0):
        if callback == 0:
            self.update_y()
        for k in range(self.nk):
            for i in range(self.ni):
                for j in range(self.ni):
                    gamma_name = ''.join(['gamma[',str(k),',',str(i),',',str(j),']'])
                    self.sub_model.getConstrByName(gamma_name).rhs = self.value_y[j]
        for k in range(self.nk):
            for j in range(self.ni):
                lamda_name = ''.join(['lamda[',str(k),',',str(j),']'])
                self.sub_model.getConstrByName(lamda_name).rhs = 1 - self.value_y[j]
        ky = [0 for k in range(self.nk)]
        for k in range(self.nk):
            ky[k] = sum([self.sk[k][j]*self.value_y[j] for j in range(self.ni)])
        for k in range(self.nk):
            nu_name = ''.join(['nu[',str(k),']'])
            self.sub_model.getConstrByName(nu_name).rhs = self.p + ky[k] - sum(self.value_y)
        self.sub_model.update()
        #return self.sub_model

    def update_y(self):
        self.value_y = []
        for j in range(self.ni):
            y_name = ''.join(['y[',str(j),']'])
            y_temp = self.master_model.getVarByName(y_name)
            self.value_y.append(y_temp.x)

    def update_cut(self):
        cm1 = self.sub_model.getConstrs()
        num_c = len(cm1)
        dual_value = []
        constrname = []
        for i in range(num_c):
            dual_value.append(cm1[i].getAttr('Pi'))
            constrname.append(cm1[i].getAttr('ConstrName'))
        dual = dict(zip(constrname,dual_value))
        gamma = [[0 for j in range(self.ni)] for i in range(self.ni)]
        delta = [[0 for j in range(self.ni)] for i in range(self.ni)]
        epsilon = [0 for j in range(self.ni)]
        lamda = [0 for j in range(self.ni)]
        mu = [0 for j in range(self.ni)]
        for i in range(self.ni):
            for j in range(self.ni):
                gamma_name = ''.join(['gamma[',str(self.max_k),',',str(i),',',str(j),']'])
                delta_name = ''.join(['delta[',str(self.max_k),',',str(i),',',str(j),']'])
                gamma[i][j] = dual[gamma_name]
                delta[i][j] = dual[delta_name]
        for n in range(self.ni):
            epsilon_name = ''.join(['epsilon[',str(self.max_k),',',str(n),']'])
            lamda_name = ''.join(['lamda[',str(self.max_k),',',str(n),']'])
            mu_name = ''.join(['mu[',str(self.max_k),',',str(n),']'])
            epsilon[n] = dual[epsilon_name]
            lamda[n] = dual[lamda_name]
            mu[n] = dual[mu_name]
        nu_name = ''.join(['nu[',str(self.max_k),']'])
        nu = dual[nu_name]
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
            ajk_y.append(self.sk[self.max_k][j]-1)
        c_y = LinExpr(gamma_y,self.y.select()*self.ni) - LinExpr(lamda,self.y.select()) \
        + nu*LinExpr(ajk_y,self.y.select())
        constant_delta = 0
        for i in range(self.ni):
            constant_delta += quicksum([(1-self.sk[self.max_k][j])*delta[i][j] for j in range(self.ni)])
        constant = quicksum(epsilon) + quicksum(lamda) + constant_delta +\
        quicksum([(1-self.sk[self.max_k][j])*mu[j] for j in range(self.ni)]) + self.p*nu
        self.constr_y = c_y + constant

    def gap_calculation(self):
        # extract L
        var_L = self.master_model.getVarByName('L')
        value_L = var_L.x
        # update LB = objective value of master problem
        obj_master = self.master_model.getObjective()
        self.LB = obj_master.getValue()
        max_L3 = self.worst_scenario()
        # update UB
        self.UB = min([self.UB,0.5*value_L + 0.5*max_L3[0]])
        self.gap = (self.UB-self.LB)/self.LB
        return self.gap

    def worst_scenario(self):
        value_L3 = []
        for k in range(self.nk):
            L3_name = ''.join(['L3[',str(k),']'])
            L3_temp = self.sub_model.getVarByName(L3_name)
            value_L3.append(L3_temp.x)
        # maximum L3 and its index (worst k)
        max_L3 = max([[v,i] for i,v in enumerate(value_L3)])
        self.max_k = max_L3[1]
        return max_L3

    def update_status(self):
        self.add_cut_scen.append(self.max_k)
        self.iteration += 1
        print('==========================================')
        print('Current iteration:',str(self.iteration))
        print('gap = ',str(self.gap))
        print('Cuts added form scenario:',str(self.add_cut_scen[0:-1]))
