'''
Benders' Decomposition:
 Branch and cut
 Multiple scenario generation
 Improved Integer cut generation
 Test: Callback frequency and times
Du Bo
'''
import model_rflp as mr
from gurobipy import *
import time
def bra_cut(p,cd,cdk,sk,a1):
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    try:
        # callback_time = 0
        # callback_num = 0
        # Benders_cut_call = 0
        # Integer_cut = 0
        def mycallback(model, where):
            if where == GRB.Callback.MIPSOL:
                TSRFLP.callback_num += 1
                print('callback_num: ',TSRFLP.callback_num)
                time1 = time.time()
                # status
                # nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
                # obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
                # solcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
                # objbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                # objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                # gap_mipsol = abs(objbst - objbnd)/(1.0 + abs(objbst))
                # #print('**** New solution at node %d, obj %g, sol %d, '
                #       'gap = %g ****' % (nodecnt, obj, solcnt, gap_mipsol))
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-2 - ni:-2]
                TSRFLP.value_y = [round(x) for x in TSRFLP.value_y] # make sure y are binary
                TSRFLP.value_omega = vals[-1]
                TSRFLP.update_sub_dual(callback=1)
                time_subdual = time.time()
                TSRFLP.sub_dual.optimize()
                max_Lk = TSRFLP.worst_scenario()
                print("DUAL_SUB_callback--- %s seconds ---" % round((time.time() - time_subdual), 2))
                if max_Lk[0] - TSRFLP.value_omega >=1e-4:
                    TSRFLP.update_multiple_scenario()
                    TSRFLP.Benders_cut_call += 1
                    print('Benders_cut: ',TSRFLP.Benders_cut_call)
                    print('Benders_total:', TSRFLP.Benders_cut)

                # ------- integer cut --------
                else:
                    TSRFLP.update_sub(callback=1)
                    time_sub = time.time()
                    TSRFLP.sub_model.optimize()
                    print("PRIMAL_SUB_callback--- %s seconds ---" % round((time.time() - time_sub), 2))
                    TSRFLP.worst_scenario(1) # calculate max L3
                    TSRFLP.gap_calculation(1) # calculate int_gap
                    # print('----Integer gap:',TSRFLP.int_gap)
                    if TSRFLP.int_gap >= 1e-4:
                        # cut incumbent solution
                        TSRFLP.update_integer_cut()
                        model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
                        TSRFLP.Integer_cut += 1
                        print('Integer_cut: ',TSRFLP.Integer_cut)
                TSRFLP.callback_time += time.time() - time1
                print('total callback time: ',TSRFLP.callback_time)
                print("callback--- %s seconds ---" % round((time.time() - time1), 2))
            if where == GRB.Callback.MESSAGE:
                # Message callback
                msg = model.cbGet(GRB.Callback.MSG_STRING)
                cutname = 'Lazy constraints'
                if cutname in msg:
                    TSRFLP.num_cut += int(msg[20:-1])
            # print(time.time() - start_time)


            if time.time() - start_time >= 2000:
                model.terminate()

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
        # #print("BUILDING MASTER--- %s seconds ---" % round((time.time() - build), 2))
        # build = time.time()
        TSRFLP.dual_sub(callback=1)
        # #print("BUILDING SUBDUAL--- %s seconds ---" % round((time.time() - build), 2))
        # build = time.time()
        TSRFLP.sub(callback=1)
        # #print("BUILDING SUB--- %s seconds ---" % round((time.time() - build), 2))
        TSRFLP.params_tuneup()
        # set initail time
        start_time = time.time()
        TSRFLP.master_model._lastnode = -GRB.INFINITY
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        #
#        TSRFLP.warm_start()
        TSRFLP.master_model.optimize(mycallback)

    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
        # TSRFLP.error_check()
    #print("--- %s seconds ---" % round((time.time() - start_time), 2))
    # OUTPUT
    runtime =  round((time.time() - start_time), 2)
    if TSRFLP.master_model.Status == 2:
        TSRFLP.opt = 1
    objval = round(TSRFLP.master_model.Objval,2)
    gap= TSRFLP.master_model.MIPGap
    if abs(gap)<=1e-5:
        gap = 0

    Benders_total = TSRFLP.Benders_cut

    var_y=[]
    for j in range(TSRFLP.ni):
        y_name = ''.join(['y[', str(j), ']'])
        y_temp = TSRFLP.master_model.getVarByName(y_name)
        var_y.append(y_temp.x)

    return var_y,runtime,TSRFLP.num_cut,TSRFLP.opt,objval,gap
