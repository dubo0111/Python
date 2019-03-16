'''
Benders' Decomposition:
 Branch and cut
 Multiple scenario generation
 Improved Integer cut generation
 Variable Neighbourhood Branching (After Root Nodes) Version 2
 Proximity search

Du Bo
'''
import model_rflp as mr
from gurobipy import *
import time
def bra_cut(p,cd,cdk,sk,a1, tl_total, tl_node,tl_pr_node,tl_pr_total,branch_step):
    convergence = []
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    try:
        def mycallback(model, where): # callback fuction: benders cut & integer cut
            if where ==GRB.Callback.MIP:
                if TSRFLP.LB_terminate == 1 and TSRFLP.LB_branch == 0:
                    if TSRFLP.vn_end == 1:
                        if time.time()-pr_node_time>=tl_pr_node:
                            model.terminate()
                        if time.time()-pr_time>=tl_pr_total:
                            model.terminate()
                            TSRFLP.pr_end = 1
                    else:
                        if time.time()-LB_time>=tl_node:
                            model.terminate()
                        if time.time()-vn_time>=tl_total:
                            model.terminate()
                            TSRFLP.vn_end = 1
                if TSRFLP.LB_branch == 1:
                    objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
                    objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
                    if objbst < 1e10:
                        convergence.append([objbst,objbnd,time.time()-start_time])
                if time.time() - start_time >= 1000: # Stop criteria
                    model.terminate()
            if where == GRB.Callback.MIPSOL:
                nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
                objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-3 - ni:-3]
                if TSRFLP.warm == 'over':
                    TSRFLP.value_y = [round(x) for x in TSRFLP.value_y] # make sure y are binary
                TSRFLP.value_omega = vals[-1]
                if nodecnt > 0 and TSRFLP.LB_terminate == 0: # LB right after root node
                    TSRFLP.LB_terminate = 1
                    TSRFLP.bestbound = objbnd
                    model.terminate()
                # if TSRFLP.value_y not in TSRFLP.tabu: # save best incumbent
                #     TSRFLP.tabu.append(TSRFLP.value_y)
                TSRFLP.update_sub_dual(callback=1)
                TSRFLP.sub_dual.optimize()
                max_Lk = TSRFLP.worst_scenario()
                if max_Lk[0] - TSRFLP.value_omega >=1e-4: # ----benders cut----
                    TSRFLP.update_multiple_scenario()
                else: # ----integer cut----
                    TSRFLP.update_sub(callback=1)
                    time_sub = time.time()
                    TSRFLP.sub_model.optimize()
                    TSRFLP.worst_scenario(1) # calculate max L3
                    TSRFLP.gap_calculation(1) # calculate int_gap
                    if TSRFLP.int_gap >= 1e-4:
                        # cut incumbent solution
                        TSRFLP.update_integer_cut()
                        model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
            if where == GRB.Callback.MESSAGE: # Record lazy constraints
                if TSRFLP.LB_branch == 1:
                    msg = model.cbGet(GRB.Callback.MSG_STRING)
                    cutname = 'Lazy constraints'
                    if cutname in msg:
                        TSRFLP.num_cut += int(msg[20:-1])

        TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk) # instantiate class
        # setting algorithm environment
        TSRFLP.dual = 1
        TSRFLP.intSP = 1.0
        TSRFLP.lift = 0
        TSRFLP.zero_half = 0
        gap = 1 #
        # initailization
        TSRFLP.dual_sub(callback=1)
        TSRFLP.sub(callback=1)
        TSRFLP.warm_start(1) # 1:no warm start 0: warm start
        TSRFLP.params_tuneup()
        start_time = time.time() # set initail time
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        TSRFLP.master_model.optimize(mycallback) # terminate after root node
        vn_time = time.time()
        # pr_time = 0
        # pr_node_time = 0
        Branching_record = [1e6,[]] # best objval; best solution;
        # Branching
        while TSRFLP.vn_end == 0: #
            LB_time = time.time() # time Limits for one neighbourhood
            TSRFLP.add_LB(branch_step)
            TSRFLP.master_model.optimize(mycallback)
            best_incumbent = []
            if TSRFLP.master_model.Objval < Branching_record[0]:
                Vars = TSRFLP.master_model.getVars()
                for n in Vars:
                    best_incumbent.append(n.x)
                Branching_record = [TSRFLP.master_model.Objval,best_incumbent]
            TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-1])
            TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-2])
        # Proximity search
        pr_time = time.time()
        pr_gap = 1
        while TSRFLP.pr_end == 0 and pr_gap > 0.2:
            pr_node_time = time.time()
            rhs,soft_rhs=TSRFLP.add_proximity(Branching_record)
            TSRFLP.master_model.optimize(mycallback)
            if TSRFLP.master_model.Status in [2,11]: # optimal or interrupted
                if TSRFLP.master_model.ObjVal < 1e10: # optimal or feasible
                    best_incumbent = []
                    obj_now = TSRFLP.a1*TSRFLP.L.x+TSRFLP.a2*TSRFLP.omega.x
                    if obj_now < Branching_record[0]:
                        # print(obj_now)
                        Vars = TSRFLP.master_model.getVars()
                        for n in Vars:
                            best_incumbent.append(n.x)
                        Branching_record = [obj_now,best_incumbent]
                    if abs(soft_rhs-obj_now) < 0.01:
                        TSRFLP.bestbound = rhs
                else: # cannot find feasible solution
                    TSRFLP.pr_end == 1 # stop
            elif TSRFLP.master_model.Status in [3,4,5]: #infeasible
                TSRFLP.bestbound = Branching_record[0]-(Branching_record[0]-TSRFLP.bestbound)/4
            pr_gap = (Branching_record[0]-TSRFLP.bestbound)/(1+Branching_record[0])
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('gap: ',pr_gap,' UB= ',Branching_record[0],' LB= ',TSRFLP.bestbound)
            TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-1])
        TSRFLP.remove_proximity()
        TSRFLP.LB_branch = 1
        for x in TSRFLP.LB_cuts:
            TSRFLP.master_model.addConstr(TSRFLP.omega >= x)
            TSRFLP.master_model.getConstrs()[-1].Lazy = 1
        if Branching_record[1] != []:
            TSRFLP.set_initial(Branching_record[1])
        TSRFLP.master_model.optimize(mycallback) # final optimization
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
    runtime =  round((time.time() - start_time), 2)
    if TSRFLP.master_model.Status == 2:
        TSRFLP.opt = 1
    objval = round(TSRFLP.master_model.Objval,2)
    gap= TSRFLP.master_model.MIPGap
    if abs(gap)<=1e-5:
        gap = 0
    var_y=[]
    for j in range(TSRFLP.ni):
        y_name = ''.join(['y[', str(j), ']'])
        y_temp = TSRFLP.master_model.getVarByName(y_name)
        var_y.append(y_temp.x)
    convergence = [*zip(*convergence)]
    return var_y,runtime,TSRFLP.num_cut,TSRFLP.opt,objval,gap,convergence,len(TSRFLP.LB_cuts)
