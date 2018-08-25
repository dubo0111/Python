# class for:
# bulid or update specific model
# get variables,dual variables
import copy
import math
from gurobipy import *


class rflp:
    # Attributes
    master_model = Model()
    sub_model = Model()
    sub_dual = Model()
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
    # cuts adding
    violation = []
    freq = []
    viol_freq = []
    lift = 0
    zero_half = 0
    # Input
    p = 0
    ni = 0
    nk = 0
    a1 = 0
    a2 = 0
    cd = []
    cdk = []
    sk = []

    def __init__(self, p, ni, nk, a1, a2, cd, cdk, sk):
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
    # build master problem model

    def master(self):
        # Create variables
        # x:allocations y:location L:auxiliary variable
        x = self.master_model.addVars(
            self.ni, self.ni, vtype=GRB.CONTINUOUS, name="x")
        self.y = self.master_model.addVars(self.ni, vtype=GRB.BINARY, name="y")
        L = self.master_model.addVar(
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
            (x.prod(cdx, i, '*') <= L for i in range(self.ni)),
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
        self.master_model.update()
        # return self.master_model
    # build subproblem model
    # (get dual variables through Constr.getAttr(Pi))

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
            (eta == L3.sum()),
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
        self.sub_model.update()
        # return(self.sub_model)
    # build dual subproblem model

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
        self.sub_dual.update()
    #

    def sub_dual_obj(self):
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
    #

    def update_master(self):
        self.update_cut()
        self.master_model.addConstr(self.omega >= self.constr_y)
        # self.master_model.update()
    #

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
    #

    def update_sub_dual(self, callback=0):
        if callback == 0:
            self.update_y()
        for k in range(self.nk):
            constr_name = ''.join(['Q(k)[', str(k), ']'])
            self.sub_dual.remove(self.sub_dual.getConstrByName(constr_name))
        self.sub_dual_obj()
        self.sub_dual.update()
    # update value of y for subproblem in each iteration

    def update_y(self):
        self.value_y = []
        for j in range(self.ni):
            y_name = ''.join(['y[', str(j), ']'])
            y_temp = self.master_model.getVarByName(y_name)
            self.value_y.append(y_temp.x)
    # update cut to be added to master problem in each iteration

    def update_cut(self, numk=None, lift=0):
        gamma = [[0 for j in range(self.ni)] for i in range(self.ni)]
        delta = [[0 for j in range(self.ni)] for i in range(self.ni)]
        epsilon = [0 for j in range(self.ni)]
        lamda = [0 for j in range(self.ni)]
        mu = [0 for j in range(self.ni)]
        beta = [0 for j in range(self.ni)]
        if numk == None:
            numk = self.max_k
        if self.dual == 0:
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
        elif self.dual == 1:
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
        # dual_lifting
        if lift != 0:
            #             sum_gamma = [0 for i in range(self.ni)]
            #             for j in range(self.ni):
            #                 for i in range(self.ni):
            #                     sum_gamma[j] += gamma[i][j]
            #             for j in range(self.ni):
            #                 if self.value_y[j] == 1:
            #                     new_lamda = -nu-mu[j]+sum_gamma[j]
            #                     if lamda[j] != new_lamda:
            # #                        print('lift lamda')
            #                         lamda[j] = new_lamda
            for i in range(self.ni):
                for j in range(self.ni):
                    if self.value_y[j] == 0:
                        new_gamma = -epsilon[i] - delta[i][j] - \
                            self.cdk[numk][i][j] * beta[i]
                        if gamma[i][j] != new_gamma:
                            #                            print('lift gamma')
                            gamma[i][j] = new_gamma
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
        self.constr_y = c_y + constant
    #

    def update_integer_cut(self):
        sum_c_y = LinExpr()
        for i in range(self.ni):
            if self.value_y[i] == 1:
                sum_c_y += self.y[i] - 1
            elif self.value_y[i] == 1:
                sum_c_y += -self.y[i]
        sum_c_y += 1
        max_Lk = self.worst_scenario(1)
        self.integer_cut = max_Lk[0] * sum_c_y
    #

    def gap_calculation(self, MIP_SP=0, Check_optimal=0):
        if MIP_SP == 1:  # in mycallback
            vals = self.master_model.cbGetSolution(self.master_model._vars)
            #value_L = vals(-2)
            max_Lk = self.worst_scenario(1)
            self.int_gap = max_Lk[0] - self.value_omega
           # print('max_Lk:',max_Lk[0])
        else:
            # extract L
            var_L = self.master_model.getVarByName('L')
            value_L = var_L.x
            # update LB = objective value of master problem
            obj_master = self.master_model.getObjective()
            self.LB = obj_master.getValue()
            if Check_optimal == 0:
                max_Lk = self.worst_scenario()
            elif Check_optimal == 1:
                max_Lk = self.worst_scenario(1)
            # update UB
            self.UB = min([self.UB, self.a1 * value_L + self.a2 * max_Lk[0]])
            self.gap = (self.UB - self.LB) / self.LB
        # return self.gap
    #

    def worst_scenario(self, MIP_SP=0):
        if self.dual == 0 or MIP_SP == 1:
            value_L3 = []
            for k in range(self.nk):
                L3_name = ''.join(['L3[', str(k), ']'])
                L3_temp = self.sub_model.getVarByName(L3_name)
                value_L3.append(L3_temp.x)
            # maximum L3 and its index (worst k)
            max_Lk = max([[v, i] for i, v in enumerate(value_L3)])
        elif self.dual == 1 and MIP_SP != 1:
            value_Qk = []
            for k in range(self.nk):
                Qk_name = ''.join(['Qk[', str(k), ']'])
                Qk_temp = self.sub_dual.getVarByName(Qk_name)
                value_Qk.append(Qk_temp.x)
            # maximum Qk and its index (worst k)
            max_Lk = max([[v, i] for i, v in enumerate(value_Qk)])
        self.max_k = max_Lk[1]
        return max_Lk
    #

    def update_status(self):
        self.add_cut_scen.append(self.max_k)
        self.iteration += 1
        print('==========================================')
        print('Current iteration:', str(self.iteration))
        print('gap = ', str(self.gap))
        print('Cuts added form scenario:', str(self.add_cut_scen[0:-1]))
    # tune parameters to avoid numerical issues for subproblem
    # wrong optimal solutions appear for both sub&dual_sub

    def params_tuneup(self):
        # self.master_model.params.OutputFlag = 0
        #        self.master_model.params.Presolve = 0
        #        self.master_model.params.ScaleFlag = 3
        #        self.master_model.params.NumericFocus = 3
        self.master_model.params.PreCrush = 1

        self.sub_model.params.OutputFlag = 0
        # self.sub_model.params.Presolve = 0
        # self.sub_model.params.ScaleFlag = 3
        # self.sub_model.params.NumericFocus = 3

        self.sub_dual.params.OutputFlag = 0
        # self.sub_dual.params.Presolve = 0
        # self.sub_dual.params.ScaleFlag = 3
        # self.sub_dual.params.NumericFocus = 3
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
    #

    def error_check(self):
        if self.gap <= -0.1:
            self.error = 1
            print('# WARNING: gap is a negative value', '\n',
                  '   Wrong dual problem solution.')
        return self.error
    # Callback: update scenario violation and frequency⁠

    def update_scenario_sorting(self):
        # extract omega and Q(k)
        violation_now = [
            x.x - self.value_omega for x in self.sub_dual.getVars()[-self.nk:]]
        for k in range(self.nk):
            if violation_now[k] > 0:
                self.violation[k] += violation_now[k]
                self.freq[k] += 1
                self.viol_freq[k] = self.violation[k] / self.freq[k]
        #print(self.viol_freq)
        rank = sorted(range(len(self.viol_freq)), reverse=True, key=self.viol_freq.__getitem__)
        rank = sorted(range(len(violation_now)), reverse=True,
                      key=violation_now.__getitem__)
        #print(rank[0])
        for n in range(round(self.nk / 4)):
            #        for n in range(1):
            if violation_now[rank[n]] > 0:
                self.update_cut(rank[n], self.lift)
                if self.zero_half == 0:
                    self.master_model.cbLazy(self.omega >= self.constr_y)
                elif self.zero_half == 1:
                    if n == 0:
                        self.master_model.cbLazy(self.omega >= self.constr_y)
                        constr_strong = self.constr_y.copy()
                        coeff_strong, constant_strong = self.reformulate_linexpr(
                            constr_strong)
                    else:
                        # zero-half
                        coeff_constr_y, constant_y = self.reformulate_linexpr(
                            self.constr_y)
                        coeff_sum = [0.5 * (x + y)
                                     for x, y in zip(coeff_strong, coeff_constr_y)]
                        coeff_zero_half = [math.ceil(x) for x in coeff_sum]
                        constant_zero_half = math.ceil(
                            0.5 * (constant_strong + constant_y))
                        constr_y_zero_half = LinExpr(
                            coeff_zero_half, self.y.select()) + constant_zero_half
                        #print(constr_y_zero_half)
                        self.master_model.cbLazy(self.omega >= constr_y_zero_half)
            else:
                break

    def reformulate_linexpr(self, formulation):
        coeff_y = [0 for j in range(self.ni)]
        for n in range(formulation.size()):
            y_index = int(formulation.getVar(n).VarName[2])
            coeff_y[y_index] += formulation.getCoeff(n)
        constant_y = formulation.getConstant()
        return coeff_y, constant_y
        # def dual_lifting(self):
        #     for j in range(ni):
        #         if self.value_y[j] = 0
        # def zero_half(self):
