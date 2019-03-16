'''
Benders' Decomposition:
 Branch and cut
 Multiple scenario generation
 Improved Integer cut generation
 Local Branching (After Root Nodes)
Backup
Du Bo
'''
import model_rflp as mr
from gurobipy import *
import time
def bra_cut(p,cd,cdk,sk,a1,tl_total,branch_step):
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
                    if time.time() - start_time >= tl_total: # time Limits
                        model.terminate()
                if TSRFLP.LB_branch == 1:
                    objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
                    objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
                    convergence.append([objbst,objbnd,time.time()-start_time])
                if time.time() - start_time >= 1000: # Stop criteria
                    model.terminate()
            if where == GRB.Callback.MIPSOL:
                # Status output
                nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-3 - ni:-3]
                if TSRFLP.warm == 'over':
                    TSRFLP.value_y = [round(x) for x in TSRFLP.value_y] # make sure y are binary
                TSRFLP.value_omega = vals[-1]
                # Local Branching
                if nodecnt > 0 and TSRFLP.LB_terminate == 0: # LB right after root node
                    TSRFLP.LB_terminate = 1
                    model.terminate()
                TSRFLP.update_sub_dual(callback=1)
                time_subdual = time.time()
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
                # Message callback
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
        #
        gap = 1 #
        TSRFLP.dual_sub(callback=1)
        TSRFLP.sub(callback=1)
        TSRFLP.warm_start(1) # 1:no warm start 0: warm start
        TSRFLP.params_tuneup()
        start_time = time.time()
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        TSRFLP.master_model.optimize(mycallback) # terminate after root node
        Branching_record = [1e6,[]]
        Branching_record,better_sol = TSRFLP.record_best_sol(Branching_record)
        LB_total_time = time.time()
        TSRFLP.add_LB(branch_step)
        while time.time() - LB_total_time < 2000: # tl_total
            TSRFLP.master_model.optimize(mycallback)
            Branching_record,better_sol = TSRFLP.record_best_sol(Branching_record)
            if TSRFLP.master_model.Status in [2]:
                branch_step += 2 #
                TSRFLP.master_model.getConstrs()[-2].rhs = branch_step-2 #
                TSRFLP.master_model.getConstrs()[-1].rhs = branch_step
                print('+++++++++++++++++branch_step: ',branch_step)
            if branch_step > TSRFLP.p*2 or TSRFLP.master_model.Status in [11]:
                break
        TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-1])
        TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-2])
        TSRFLP.LB_branch = 1
        # add all lazy cuts
        # print('Cuts to be added:   ',len(TSRFLP.LB_cuts))
        for x in TSRFLP.LB_cuts:
            TSRFLP.master_model.addConstr(TSRFLP.omega >= x)
            TSRFLP.master_model.getConstrs()[-1].Lazy = 1
        TSRFLP.set_initial(Branching_record[1])
        TSRFLP.master_model.optimize(mycallback) # final optimization

    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
        # TSRFLP.error_check()
    # OUTPUT
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
    return var_y,runtime,TSRFLP.num_cut,TSRFLP.opt,objval,gap,convergence
