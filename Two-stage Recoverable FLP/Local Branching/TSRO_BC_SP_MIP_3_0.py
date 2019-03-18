'''
Benders' Decomposition:
 Branch and cut
 Multiple scenario generation
 !Improved Integer cut generation

Du Bo
'''
import model_rflp as mr
from gurobipy import *
import time
def bra_cut(p,cd,cdk,sk,a1):
    convergence = []
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    try:
        def mycallback(model, where): # callback fuction: benders cut & integer cut
            # time1 = time.time()
            if where ==GRB.Callback.MIP:
                objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
                objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
                if objbst < 1e10:
                    if convergence != []:
                        if objbst < convergence[-1][0] :
                            convergence.append([convergence[-1][0],convergence[-1][1],time.time()-start_time])
                            convergence.append([objbst,objbnd,time.time()-start_time])
                        else:
                            convergence.append([objbst,objbnd,time.time()-start_time])
                    else:
                        convergence.append([objbst,objbnd,time.time()-start_time])
                if time.time() - start_time >= 1000:
                    model.terminate()
            if where == GRB.Callback.MIPSOL:
                # objbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                # objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                # if objbst < 1e10 and objbnd>0:
                #     if convergence != []:
                #         if objbst < convergence[-1][0] :
                #             convergence.append([convergence[-1][0],convergence[-1][1],time.time()-start_time])
                #             convergence.append([objbst,objbnd,time.time()-start_time])
                #         else:
                #             convergence.append([objbst,objbnd,time.time()-start_time])
                #     else:
                #         convergence.append([objbst,objbnd,time.time()-start_time])
                vals = model.cbGetSolution(model._vars)
                TSRFLP.value_y = vals[-3 - ni:-3]
                if TSRFLP.warm == 'over':
                    TSRFLP.value_y = [round(x) for x in TSRFLP.value_y] # make sure y are binary
                TSRFLP.value_omega = vals[-1]
                if TSRFLP.value_y not in TSRFLP.save_y and TSRFLP.value_y not in TSRFLP.save_y_int:
                    TSRFLP.update_sub_dual(callback=1)
                    TSRFLP.sub_dual.optimize()
                    max_Lk = TSRFLP.worst_scenario()
                    SP_Qk = [i.x for i in TSRFLP.sub_dual.getVars()[-TSRFLP.nk:]]
                    save_sub = TSRFLP.get_subdual_vals()
                    TSRFLP.save_max_Lk_DualLP.append([TSRFLP.value_y,TSRFLP.max_Lk,SP_Qk,save_sub])
                    TSRFLP.save_y.append(TSRFLP.value_y)
                    if max_Lk[0] - TSRFLP.value_omega >=1e-4: # ----benders cut----
                        TSRFLP.update_multiple_scenario(SP_Qk)
                    else: # ----integer cut----
                        TSRFLP.update_sub(callback=1)
                        TSRFLP.sub_model.optimize()
                        TSRFLP.worst_scenario(1) # calculate max L3
                        TSRFLP.save_max_Lk_SP.append([TSRFLP.value_y,TSRFLP.max_Lk])
                        TSRFLP.save_y_int.append(TSRFLP.value_y)
                        TSRFLP.gap_calculation(1) # calculate int_gap
                        if TSRFLP.int_gap >= 1e-4:
                            TSRFLP.update_integer_cut()
                            model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
                else:
                    save_index = [(i, x.index(TSRFLP.value_y)) for i, x in enumerate(TSRFLP.save_max_Lk_DualLP) if TSRFLP.value_y in x]
                    if TSRFLP.save_max_Lk_DualLP[save_index[0][0]][1][0] - TSRFLP.value_omega >=1e-4: # ----benders cut----
                        TSRFLP.update_multiple_scenario(TSRFLP.save_max_Lk_DualLP[save_index[0][0]][2],
                                                        TSRFLP.save_max_Lk_DualLP[save_index[0][0]][3])
                    else:
                        save_index = [(i, x.index(TSRFLP.value_y)) for i, x in enumerate(TSRFLP.save_max_Lk_SP) if TSRFLP.value_y in x]
                        if save_index != []:
                            if TSRFLP.save_max_Lk_SP[save_index[0][0]][1][0]-TSRFLP.value_omega >= 1e-4:
                                TSRFLP.update_integer_cut(0,TSRFLP.save_max_Lk_SP[save_index[0][0]][1])
                                model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
                        else:
                            TSRFLP.update_sub(callback=1)
                            TSRFLP.sub_model.optimize()
                            TSRFLP.worst_scenario(1) # calculate max L3
                            TSRFLP.save_max_Lk_SP.append([TSRFLP.value_y,TSRFLP.max_Lk])
                            TSRFLP.save_y_int.append(TSRFLP.value_y)
                            TSRFLP.gap_calculation(1) # calculate int_gap
                            if TSRFLP.int_gap >= 1e-4:
                                TSRFLP.update_integer_cut()
                                model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
#                    print('sadasdadsasdadasdadsasdadasdadsasd')
                    # indices = [i for i, x in enumerate(TSRFLP.LB_cuts_y) if x == TSRFLP.value_y]
                    # for n in indices:
                    #     model.cbLazy(TSRFLP.omega >= TSRFLP.LB_cuts[n])
            if where == GRB.Callback.MESSAGE: # lazy constraints
                # Message callback
                msg = model.cbGet(GRB.Callback.MSG_STRING)
                cutname = 'Lazy constraints'
                if cutname in msg:
                    TSRFLP.num_cut += int(msg[20:-1])
            # print(time.time() - start_time)


        TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk) # instantiate class
        # setting algorithm environment
        TSRFLP.dual = 1
        TSRFLP.intSP = 1.0
        TSRFLP.lift = 0
        TSRFLP.zero_half = 0
        #
        gap = 1 #
        # stop = 1e-5

        # build = time.time()
        # TSRFLP.master() # needed when .warm_start is turned off
        # print("BUILDING MASTER--- %s seconds ---" % round((time.time() - build), 2))

        # build = time.time()
        TSRFLP.dual_sub(callback=1)
        # print("BUILDING SUBDUAL--- %s seconds ---" % round((time.time() - build), 2))

        # build = time.time()
        TSRFLP.sub(callback=1)
        # print("BUILDING SUB--- %s seconds ---" % round((time.time() - build), 2))

        # warm_t = time.time()
        TSRFLP.warm_start(1)
        # print("Warm time %s seconds" % round((time.time() - warm_t), 2))

        TSRFLP.params_tuneup()
        # set initail time
        start_time = time.time()

        # TSRFLP.master_model._lastnode = -GRB.INFINITY
        TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
        TSRFLP.master_model.Params.lazyConstraints = 1
        #
        TSRFLP.master_model.optimize(mycallback)
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')
        # TSRFLP.error_check()
    # print("--- %s seconds ---" % round((time.time() - start_time), 2))
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
