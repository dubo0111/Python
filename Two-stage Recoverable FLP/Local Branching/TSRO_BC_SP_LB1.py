'''
Benders' Decomposition:
 Branch and cut
 Multiple scenario generation
 Improved Integer cut generation
 Variable Neighbourhood Branching (After Root Nodes) Version 2


Du Bo
'''
import model_rflp as mr
from gurobipy import *
import time
import copy


def bra_cut(Time_Limit, p, cd, cdk, sk, a1, tl_total, tl_node, branch_step):
    convergence = []
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    try:
        def mycallback(model, where):  # callback fuction: benders cut & integer cut
            if where == GRB.Callback.MIP:
                if TSRFLP.LB_terminate == 1 and TSRFLP.LB_branch == 0:
                    if time.time() - LB_time >= tl_node:
                        model.terminate()
                    if time.time() - VN_time >= tl_total:
                        model.terminate()
                        TSRFLP.vn_end = 1
                # if TSRFLP.LB_branch == 1:
                objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
                objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
                if objbst < 1e10:
                    if TSRFLP.LB_terminate == 1 and TSRFLP.vn_end == 0:
                        # print('22222222222222222')
                        if convergence != [] and objbst <= convergence[-1][0]:
                            convergence.append(
                                [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                            convergence.append(
                                [objbst, TSRFLP.bestbound, time.time() - start_time])
                        if convergence[-1][0] == 742.5:
                            aaa =1
                    if TSRFLP.vn_end == 1 and TSRFLP.LB_branch == 1:
                        # print('3333333333333333')
                        if convergence != [] and objbnd >= TSRFLP.bestbound and objbst <= convergence[-1][0]:
                            convergence.append(
                                [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                            convergence.append(
                                [objbst, objbnd, time.time() - start_time])
                        if convergence[-1][0] == 742.5:
                            aaa =1
                    if TSRFLP.LB_terminate == 0:
                        # print('111111111111111')
                        if convergence != []:
                            convergence.append(
                                [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                        convergence.append(
                            [objbst, objbnd, time.time() - start_time])
                        if convergence[-1][0] == 742.5:
                            aaa =1
                nodecnt = model.cbGet(GRB.Callback.MIP_NODCNT)
                if time.time() - start_time >= Time_Limit and nodecnt > 0:  # Stop criteria
                    model.terminate()
            if where == GRB.Callback.MIPSOL:
                nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
                objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-3 - ni:-3]
                if TSRFLP.warm == 'over':
                    # make sure y are binary
                    TSRFLP.value_y = [round(x) for x in TSRFLP.value_y]
                TSRFLP.value_omega = vals[-1]
                if nodecnt > 200 and TSRFLP.LB_terminate == 0:  # LB right after root node
                    TSRFLP.LB_terminate = 1
                    TSRFLP.bestbound = copy.deepcopy(objbnd)
                    model.terminate()
                if TSRFLP.value_y not in TSRFLP.save_y and TSRFLP.value_y not in TSRFLP.save_y_int:
                    TSRFLP.update_sub_dual(callback=1)
                    TSRFLP.sub_dual.optimize()
                    max_Lk = TSRFLP.worst_scenario()
                    SP_Qk = [i.x for i in TSRFLP.sub_dual.getVars()
                                                                  [-TSRFLP.nk:]]
                    save_sub = TSRFLP.get_subdual_vals()
                    TSRFLP.save_max_Lk_DualLP.append(
                        [TSRFLP.value_y, TSRFLP.max_Lk, SP_Qk, save_sub])
                    TSRFLP.save_y.append(TSRFLP.value_y)
                    if max_Lk[0] - TSRFLP.value_omega >= 1e-4:  # ----benders cut----
                        TSRFLP.update_multiple_scenario(SP_Qk)
                    else:  # ----integer cut----
                        TSRFLP.update_sub(callback=1)
                        TSRFLP.sub_model.optimize()
                        TSRFLP.worst_scenario(1)  # calculate max L3
                        TSRFLP.save_max_Lk_SP.append(
                            [TSRFLP.value_y, TSRFLP.max_Lk])
                        TSRFLP.save_y_int.append(TSRFLP.value_y)
                        TSRFLP.gap_calculation(1)  # calculate int_gap
                        if TSRFLP.int_gap >= 1e-4:
                            TSRFLP.update_integer_cut()
                            model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
                else:
                    save_index = [(i, x.index(TSRFLP.value_y)) for i, x in enumerate(
                        TSRFLP.save_max_Lk_DualLP) if TSRFLP.value_y in x]
                    # ----benders cut----
                    if TSRFLP.save_max_Lk_DualLP[save_index[0][0]][1][0] - TSRFLP.value_omega >= 1e-4:
                        TSRFLP.update_multiple_scenario(TSRFLP.save_max_Lk_DualLP[save_index[0][0]][2],
                                                        TSRFLP.save_max_Lk_DualLP[save_index[0][0]][3])
                    else:
                        save_index = [(i, x.index(TSRFLP.value_y)) for i, x in enumerate(
                            TSRFLP.save_max_Lk_SP) if TSRFLP.value_y in x]
                        if save_index != []:
                            if TSRFLP.save_max_Lk_SP[save_index[0][0]][1][0] - TSRFLP.value_omega >= 1e-4:
                                TSRFLP.update_integer_cut(
                                    0, TSRFLP.save_max_Lk_SP[save_index[0][0]][1])
                                model.cbLazy(TSRFLP.omega >=
                                             TSRFLP.integer_cut)
                        else:
                            TSRFLP.update_sub(callback=1)
                            TSRFLP.sub_model.optimize()
                            TSRFLP.worst_scenario(1)  # calculate max L3
                            TSRFLP.save_max_Lk_SP.append(
                                [TSRFLP.value_y, TSRFLP.max_Lk])
                            TSRFLP.save_y_int.append(TSRFLP.value_y)
                            TSRFLP.gap_calculation(1)  # calculate int_gap
                            if TSRFLP.int_gap >= 1e-4:
                                TSRFLP.update_integer_cut()
                                model.cbLazy(TSRFLP.omega >=
                                             TSRFLP.integer_cut)
            if where == GRB.Callback.MESSAGE:  # Record lazy constraints
                # Message callback
                if TSRFLP.LB_branch == 1:
                    msg = model.cbGet(GRB.Callback.MSG_STRING)
                    cutname = 'Lazy constraints'
                    if cutname in msg:
                        TSRFLP.num_cut += int(msg[20:-1])

        TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)  # instantiate class
        # setting algorithm environment
        TSRFLP.dual = 1
        TSRFLP.intSP = 1.0
        TSRFLP.lift = 0
        TSRFLP.zero_half = 0
        #
        gap = 1
        TSRFLP.dual_sub(callback=1)
        TSRFLP.sub(callback=1)
        TSRFLP.warm_start(1)  # 1:no warm start 0: warm start
        TSRFLP.params_tuneup()
        start_time = time.time()  # set initail time
        # TSRFLP.master_model._lastnode = -GRB.INFINITY
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        TSRFLP.master_model.optimize(mycallback)  # terminate after root node
        if TSRFLP.master_model.status == 2:
            TSRFLP.LB_terminate = 1
        TSRFLP.master_model.addConstr(
            TSRFLP.a1 * TSRFLP.L + TSRFLP.a2 * TSRFLP.omega >= TSRFLP.bestbound)
        TSRFLP.Branching_record = [1e6, []]
        TSRFLP.Branching_record, better_sol, convergence, _ = TSRFLP.record_best_sol(
            TSRFLP.Branching_record, start_time, convergence)
        TSRFLP.add_LB(TSRFLP.Branching_record, branch_step, 1)
        LB_cut = 2
        VN_time = time.time()
        while TSRFLP.vn_end == 0 and tl_total != 0:
            LB_time = time.time()  # time Limits for one neighbourhood
            TSRFLP.master_model.optimize(mycallback)
            if TSRFLP.vn_end == 1:
                break
            if TSRFLP.master_model.status in [3,4,5]:
                if branch_step <= TSRFLP.p*2-2:
                    if branch_step == 2:
                        TSRFLP.reverse_LB(TSRFLP.Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(TSRFLP.Branching_record,branch_step+2) #
                        LB_cut += 1
                        # print('*********************','++first reverse++',branch_step,TSRFLP.Branching_record[0],'*********************')
                    else:
                        TSRFLP.reverse_LB(TSRFLP.Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(TSRFLP.Branching_record,branch_step+2) #
                        # print('*********************','++reverse++',branch_step,'incumbent',TSRFLP.Branching_record[0],'*********************')
                else:
                    TSRFLP.vn_end = 1
                    break
            if TSRFLP.master_model.status in [2,11]:
                TSRFLP.Branching_record,better_sol,convergence,Reverse_record = TSRFLP.record_best_sol(
                                                        TSRFLP.Branching_record,start_time,convergence)
                if better_sol == 1:
                    TSRFLP.reverse_LB(Reverse_record,branch_step,better_sol=1)
                    branch_step = 2
                    TSRFLP.add_LB(TSRFLP.Branching_record,branch_step) #
                    LB_cut += 1
                    # print('*********************','++better_sol++',branch_step,'incumbent',TSRFLP.Branching_record[0],'last one',Reverse_record[0],'*********************')
                else:
                    if branch_step == 2:
                        TSRFLP.reverse_LB(TSRFLP.Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(TSRFLP.Branching_record,branch_step+2) #
                        LB_cut += 1
                        # print('*********************','++first reverse++',branch_step,'incumbent',TSRFLP.Branching_record[0],'*********************')
                    else:
                        TSRFLP.reverse_LB(TSRFLP.Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(TSRFLP.Branching_record,branch_step+2) #
                        # print('*********************','++reverse++',branch_step,'incumbent',TSRFLP.Branching_record[0],'*********************')
        # extract heuristics solution
        for n in range(len(convergence)):
            if abs(TSRFLP.Branching_record[0] - convergence[n][0])<=1e-5:
                Heu_sol = [round(TSRFLP.Branching_record[0],2),round(convergence[n][2],2)]
                break
        for n in range(LB_cut):
            TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-n-1])
        TSRFLP.LB_branch = 1
#        for x in TSRFLP.LB_cuts:
#            TSRFLP.master_model.addConstr(TSRFLP.omega >= x)
#            TSRFLP.master_model.getConstrs()[-1].Lazy = 1
        for x in TSRFLP.LB_cuts:
            TSRFLP.master_model.addConstr(TSRFLP.omega >= x)
        TSRFLP.master_model.update()
        for n in range(len(TSRFLP.LB_cuts)):
            TSRFLP.master_model.getConstrs()[-1-n].Lazy = 1
        TSRFLP.master_model.update()
        if TSRFLP.Branching_record[1] != []:
            TSRFLP.set_initial(TSRFLP.Branching_record[1])
        TSRFLP.master_model.addConstr(TSRFLP.a1*TSRFLP.L+TSRFLP.a2*TSRFLP.omega >= TSRFLP.bestbound)
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
    if objval < 1e10:
        for j in range(TSRFLP.ni):
            y_name = ''.join(['y[', str(j), ']'])
            y_temp = TSRFLP.master_model.getVarByName(y_name)
            var_y.append(y_temp.x)
    convergence = [*zip(*convergence)]
    gap=round(gap,2)
    Heu_sol.append(round(((Heu_sol[0]-TSRFLP.master_model.Objval)/(1+Heu_sol[0])),2))
    return var_y,runtime,TSRFLP.num_cut,TSRFLP.opt,objval,gap,convergence,len(TSRFLP.LB_cuts),Heu_sol
