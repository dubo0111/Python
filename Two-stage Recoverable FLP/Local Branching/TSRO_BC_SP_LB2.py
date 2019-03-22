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


def bra_cut(p, cd, cdk, sk, a1, tl_total, tl_node, tl_pr_node, tl_pr_total, branch_step, stop_gap, pr_terminate, pr_step):
    convergence = []
    Heu_sol = []
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    try:
        def mycallback(model, where):  # callback fuction: benders cut & integer cut
            if where == GRB.Callback.MIP:
                if TSRFLP.LB_terminate == 1 and TSRFLP.LB_branch == 0:
                    if TSRFLP.vn_end == 1:
                        if time.time() - pr_node_time >= tl_pr_node:
                            model.terminate()
                        if time.time() - pr_time >= tl_pr_total:
                            # model.terminate()
                            TSRFLP.pr_end = 1
                    else:
                        if time.time() - LB_time >= tl_node:
                            model.terminate()
                        if time.time() - vn_time >= tl_total:
                            model.terminate() #
                            TSRFLP.vn_end = 1
                objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
                objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
                if objbst < 1e10:
                    if TSRFLP.LB_terminate == 1 and TSRFLP.vn_end == 0:
                        # print('22222222222222222')
                        if convergence != []:
                            convergence.append(
                                [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                        convergence.append(
                            [objbst, TSRFLP.bestbound, time.time() - start_time])
                    # if TSRFLP.vn_end == 1 and TSRFLP.pr_end == 0:
                    #     print('3333333333333333')
                    #     print(TSRFLP.a1*TSRFLP.value_L+TSRFLP.a2*TSRFLP.value_omega)
                    #     if convergence[-1][0] == 328.49999999999983:
                    #         aaa = 1
                    #     obj_now=TSRFLP.a1*TSRFLP.value_L+TSRFLP.a2*TSRFLP.value_omega
                    #     if obj_now <= convergence[-1][0]:
                    #         convergence.append([convergence[-1][0],convergence[-1][1],time.time()-start_time])
                    #         convergence.append([obj_now,TSRFLP.bestbound,time.time()-start_time])
                    if TSRFLP.pr_end == 1 and TSRFLP.LB_branch == 1:
                        # print('444444444444444444')
                        if convergence != [] and objbnd >= TSRFLP.bestbound and objbst <= convergence[-1][0]:
                            convergence.append(
                                [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                            convergence.append(
                                [objbst, objbnd, time.time() - start_time])
                    if TSRFLP.LB_terminate == 0:
                        # print('111111111111111')
                        if convergence != []:
                            convergence.append(
                                [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                        convergence.append(
                            [objbst, objbnd, time.time() - start_time])
                nodecnt = model.cbGet(GRB.Callback.MIP_NODCNT)
                if time.time() - start_time >= 1000 and nodecnt > 0:  # Stop criteria
                    model.terminate()
            if where == GRB.Callback.MIPSOL:
                nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
                objbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)

                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-3 - ni:-3]
                if TSRFLP.warm == 'over':
                    # make sure y are binary
                    TSRFLP.value_y = [round(x) for x in TSRFLP.value_y]
                TSRFLP.value_omega = vals[-1]
                TSRFLP.value_L = vals[-2]  # for recording objval in prox
                if nodecnt > 100 and TSRFLP.LB_terminate == 0:  # LB right after root node
                    TSRFLP.LB_terminate = 1
                    TSRFLP.bestbound = objbnd
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
                            if TSRFLP.pr_end == 0 and TSRFLP.vn_end == 1:
                                if objbst < pr_terminate:  # terninate at hard incumbent
                                    TSRFLP.pr_end = 1
#                                    model.terminate()
                                obj_now = TSRFLP.a1 * TSRFLP.value_L + TSRFLP.a2 * TSRFLP.value_omega
                                if obj_now <= convergence[-1][0]:
                                    convergence.append(
                                        [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                                    convergence.append(
                                        [obj_now, TSRFLP.bestbound, time.time() - start_time])
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
                            else:
                                if TSRFLP.pr_end == 0 and TSRFLP.vn_end == 1:
                                    if objbst < pr_terminate:  # terninate at hard/soft incumbent
                                        TSRFLP.pr_end = 1
#                                        model.terminate()
                                    obj_now = TSRFLP.a1 * TSRFLP.value_L + TSRFLP.a2 * TSRFLP.value_omega
                                    if obj_now <= convergence[-1][0]:
                                        convergence.append(
                                            [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                                        convergence.append(
                                            [obj_now, TSRFLP.bestbound, time.time() - start_time])

            if where == GRB.Callback.MESSAGE:  # Record lazy constraints
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
        gap = 1
        # initailization
        TSRFLP.dual_sub(callback=1)
        TSRFLP.sub(callback=1)
        TSRFLP.warm_start(1)  # 1:no warm start 0: warm start
        TSRFLP.params_tuneup()
        start_time = time.time()  # set initail time
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        TSRFLP.master_model.optimize(mycallback)  # terminate after root node
        TSRFLP.master_model.addConstr(
            TSRFLP.a1 * TSRFLP.L + TSRFLP.a2 * TSRFLP.omega >= TSRFLP.bestbound)
        rootval = TSRFLP.master_model.objval
        Branching_record = [1e6, []]
        Branching_record, better_sol, convergence,_ = TSRFLP.record_best_sol(
            Branching_record, start_time, convergence)
        TSRFLP.add_LB(Branching_record, branch_step, 1)
        LB_cut = 2
        vn_time = time.time()
        LB_time = 0
        # Branching
        while TSRFLP.vn_end == 0 and tl_total != 0:
            LB_time = time.time()  # time Limits for one neighbourhood
            TSRFLP.master_model.optimize(mycallback)
            if TSRFLP.master_model.status in [3,4,5]:
                if branch_step <= TSRFLP.p*2-2:
                    if branch_step == 2:
                        TSRFLP.reverse_LB(Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(Branching_record,branch_step+2) #
                        LB_cut += 1
                        print('*********************','++first reverse++',branch_step,Branching_record[0],'*********************')
                    else:
                        TSRFLP.reverse_LB(Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(Branching_record,branch_step+2) #
                        print('*********************','++reverse++',branch_step,'incumbent',Branching_record[0],'*********************')
                else:
                    TSRFLP.vn_end = 1
                    break
            if TSRFLP.master_model.status in [2,11]:
                Branching_record,better_sol,convergence,Reverse_record = TSRFLP.record_best_sol(
                                                        Branching_record,start_time,convergence)
                if better_sol == 1:
                    TSRFLP.reverse_LB(Reverse_record,branch_step,better_sol=1)
                    branch_step = 2
                    TSRFLP.add_LB(Branching_record,branch_step) #
                    LB_cut += 1
                    print('*********************','++better_sol++',branch_step,'incumbent',Branching_record[0],'last one',Reverse_record[0],'*********************')
                else:
                    if branch_step == 2:
                        TSRFLP.reverse_LB(Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(Branching_record,branch_step+2) #
                        LB_cut += 1
                        print('*********************','++first reverse++',branch_step,'incumbent',Branching_record[0],'*********************')
                    else:
                        TSRFLP.reverse_LB(Branching_record,branch_step)
                        branch_step += 2
                        TSRFLP.add_LB(Branching_record,branch_step+2) #
                        print('*********************','++reverse++',branch_step,'incumbent',Branching_record[0],'*********************')
        for n in range(LB_cut):
            TSRFLP.master_model.remove(
                TSRFLP.master_model.getConstrs()[-n - 1])
        if tl_total == 0:
            TSRFLP.vn_end = 1
        # Proximity search
        pr_time = time.time()
        pr_gap = 1
        while TSRFLP.pr_end == 0 and tl_pr_total != 0:
            pr_node_time = time.time()
            rhs, soft_rhs = TSRFLP.add_proximity(Branching_record, pr_step)
            TSRFLP.master_model.optimize(mycallback)
            if TSRFLP.master_model.Status in [2, 11]:  # optimal or interrupted
                if TSRFLP.master_model.ObjVal < 1e10:  # optimal or feasible
                    best_incumbent = []
                    obj_now = TSRFLP.a1 * TSRFLP.L.x + TSRFLP.a2 * TSRFLP.omega.x
                    if obj_now < Branching_record[0]:
                        Vars = TSRFLP.master_model.getVars()
                        for n in Vars:
                            best_incumbent.append(n.x)
                        Branching_record = [
                            obj_now, best_incumbent, time.time() - start_time]
                    if abs(soft_rhs - obj_now) < abs(rhs-obj_now) and TSRFLP.master_model.Status ==2:
                        TSRFLP.bestbound = rhs
                    convergence.append(
                        [convergence[-1][0], convergence[-1][1], time.time() - start_time])
                    convergence.append(
                        [Branching_record[0], TSRFLP.bestbound, time.time() - start_time])
                else:  # cannot find feasible solution
                    if TSRFLP.master_model.ObjBound > 1e5:
                        TSRFLP.bestbound = rhs
                    TSRFLP.pr_end = 1  # stop
            if TSRFLP.master_model.Status in [3, 4, 5]:  # infeasible
                TSRFLP.bestbound = soft_rhs
            pr_gap = (
                Branching_record[0] - TSRFLP.bestbound) / (1 + Branching_record[0])
            if pr_gap <= stop_gap:
                TSRFLP.pr_end = 1
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('gap: ', pr_gap, ' UB= ',
                  Branching_record[0], ' LB= ', TSRFLP.bestbound)
            TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-1])
            TSRFLP.master_model.remove(TSRFLP.master_model.getConstrs()[-2])
        TSRFLP.remove_proximity()
        if tl_pr_total == 0:
            TSRFLP.pr_end = 1
        for n in range(len(convergence)):
            if abs(Branching_record[0] - convergence[n][0]) <= 1e-5:
                Heu_sol = [round(Branching_record[0], 2),
                           round(convergence[n][2], 2)]
                break
        TSRFLP.LB_branch = 1
        for x in TSRFLP.LB_cuts:
            TSRFLP.master_model.addConstr(TSRFLP.omega >= x)
        TSRFLP.master_model.update()
        for n in range(len(TSRFLP.LB_cuts)):
            TSRFLP.master_model.getConstrs()[-1-n].Lazy = 1
        TSRFLP.master_model.update()
        if Branching_record[1] != []:
            TSRFLP.set_initial(Branching_record[1])
        TSRFLP.master_model.addConstr(
            TSRFLP.a1 * TSRFLP.L + TSRFLP.a2 * TSRFLP.omega >= TSRFLP.bestbound)
        TSRFLP.master_model.optimize(mycallback)  # final optimization
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
    runtime = round((time.time() - start_time), 2)
    if TSRFLP.master_model.Status == 2:
        TSRFLP.opt = 1
    objval = round(TSRFLP.master_model.Objval, 2)
    gap = TSRFLP.master_model.MIPGap
    if abs(gap) <= 1e-5:
        gap = 0
    var_y = []
    if objval < 1e10: # prevent infeasible
        for j in range(TSRFLP.ni):
            y_name = ''.join(['y[', str(j), ']'])
            y_temp = TSRFLP.master_model.getVarByName(y_name)
            var_y.append(y_temp.x) #
    convergence = [*zip(*convergence)]
    gap = round(gap, 2)
    rootval = round(rootval, 2)
    Heu_sol.append(
        round((Heu_sol[0] - TSRFLP.master_model.Objval) / (1 + Heu_sol[0]), 2))
    return var_y, runtime, TSRFLP.num_cut, TSRFLP.opt, objval, gap, convergence, len(TSRFLP.LB_cuts), Heu_sol, rootval
