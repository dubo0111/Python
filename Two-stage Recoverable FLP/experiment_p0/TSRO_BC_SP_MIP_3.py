# Branch and cut
# multiple scenario generation
import model_rflp as mr
#import data_generator1 as dg
#p, cd, cdk, sk = dg.ins_k(20, 100, 40)  # (ni,nk,randomseed*)
# import data_generator0 as dg0
# data = dg0.data_gen(20,20,2)
# p,cd,cdk,sk = data.data()
from gurobipy import *
import time
def bra_cut(p,cd,cdk,sk,a1):
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    try:
        # if gap doesn't change for * iterations, add fractional cut
        def mycallback(model, where):
            # time1 = time.time()
            if where == GRB.Callback.MIPSOL:
                # status
                # nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
                # obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
                # solcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
                # objbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                # objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                # gap_mipsol = abs(objbst - objbnd)/(1.0 + abs(objbst))
                # print('**** New solution at node %d, obj %g, sol %d, '
                #       'gap = %g ****' % (nodecnt, obj, solcnt, gap_mipsol))
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-2 - ni:-2]
                TSRFLP.value_omega = vals[-1]
                TSRFLP.update_sub_dual(callback=1)
                # time_subdual = time.time()
                TSRFLP.sub_dual.optimize()
                # print("SUB_callback--- %s seconds ---" % round((time.time() - time_subdual), 2))
                TSRFLP.update_multiple_scenario()
                # print("callback--- %s seconds ---" % round((time.time() - time1), 2))
            if where == GRB.Callback.MESSAGE:
                # Message callback
                msg = model.cbGet(GRB.Callback.MSG_STRING)
                cutname = 'Lazy constraints'
                if cutname in msg:
                    TSRFLP.num_cut += float(msg[-2])

        def mycallback_int(model,where):
            if where == GRB.Callback.MIPSOL:
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-2 - ni:-2]
                # optimality cut (Benders' cut)
                TSRFLP.value_omega = vals[-1]
                TSRFLP.update_sub_dual(callback=1)
                TSRFLP.sub_dual.optimize()
                TSRFLP.worst_scenario()
                TSRFLP.update_cut()
                model.cbLazy(TSRFLP.omega >= TSRFLP.constr_y)
                # integer l-shaped cut
                TSRFLP.update_sub(callback=1)
                TSRFLP.sub_model.optimize()
                TSRFLP.worst_scenario(1)
                TSRFLP.gap_calculation(1)
                print('---------gap:',TSRFLP.int_gap)
                if abs(TSRFLP.int_gap) >= 1e-4:
                   TSRFLP.update_integer_cut()
                   # cut incumbent solution
                   model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
            if where == GRB.Callback.MESSAGE:
                # Message callback
                msg = model.cbGet(GRB.Callback.MSG_STRING)
                cutname = 'Lazy constraints'
                if cutname in msg:
                    TSRFLP.num_cut += int(msg[-2])

        TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)
        TSRFLP.dual = 1
        TSRFLP.intSP = 1
        TSRFLP.lift = 0
        TSRFLP.zero_half = 0
        # --------------------
        gap = 1
        stop = 1e-5
        # build = time.time()
        TSRFLP.master()
        # print("BUILDING MASTER--- %s seconds ---" % round((time.time() - build), 2))
        # build = time.time()
        TSRFLP.dual_sub(callback=1)
        # print("BUILDING SUBDUAL--- %s seconds ---" % round((time.time() - build), 2))
        # build = time.time()
        TSRFLP.sub(callback=1)
        # print("BUILDING SUB--- %s seconds ---" % round((time.time() - build), 2))
        TSRFLP.params_tuneup()
        # set initail time
        start_time = time.time()
        TSRFLP.master_model._lastnode = -GRB.INFINITY
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        TSRFLP.master_model.Params.TimeLimit = 2000 # 2000 seconds
        # warm start
        TSRFLP.master_model.optimize()
        TSRFLP.update_sub_dual(0)
        TSRFLP.sub_dual.optimize()
        TSRFLP.gap_calculation()
        TSRFLP.master_model.addConstr(TSRFLP.a1*TSRFLP.master_model.getVars()[-1]+TSRFLP.a1*TSRFLP.omega >= TSRFLP.LB)

        TSRFLP.master_model.optimize(mycallback)

        # start integer cut
        TSRFLP.update_sub(callback=0)
        TSRFLP.sub_model.optimize()
        TSRFLP.worst_scenario(1)
        TSRFLP.gap_calculation(0,1)
        if abs(TSRFLP.gap) > 1e-4:
            # print('=========== Start Integer L-shaped cut =========')
            TSRFLP.master_model.reset()
            TSRFLP.master_model.optimize(mycallback_int)
        print('Optimal solution found: %g' % TSRFLP.master_model.objVal)
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
        # TSRFLP.error_check()
    print("--- %s seconds ---" % round((time.time() - start_time), 2))
    # OUTPUT
    runtime =  round((time.time() - start_time), 2)
    if TSRFLP.master_model.Status == 2:
        TSRFLP.opt = 1
    objval = round(TSRFLP.master_model.Objval,2)
    gap= TSRFLP.master_model.MIPGap
    if abs(gap)<=1e-5:
        gap = 0
    return runtime,TSRFLP.num_cut,TSRFLP.opt,objval,gap
