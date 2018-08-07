# test class
import model_rflp as mr
import data_generator1 as dg
p,cd,cdk,sk = dg.ins_k(10,10,3) #(ni,nk,randomseed*)
from gurobipy import *
import time
# Number of nodes
ni = len(cd)
nk = len(cdk)
# weights of two stages
a1 = 0.5
a2 = 1 - a1
try:
    def mycallback(model, where):
        if where == GRB.Callback.MIPSOL:
            vals = model.cbGetSolution(model._vars)
            TSRFLP.value_y = vals[-2-ni:-2]
            TSRFLP.update_sub(callback = 1)
            TSRFLP.sub_model.optimize()
            TSRFLP.worst_scenario()
            print('++++++++++++++++++++++++++++++')
            TSRFLP.update_cut()
            model.cbLazy(TSRFLP.omega >= TSRFLP.constr_y)

    start_time = time.time()
    TSRFLP = mr.rflp(p,ni,nk,a1,a2,cd,cdk,sk)
    # ----------Benders' Decompisition----------
    iteration = 0
    gap = 1
    stop = 1e-5
    add_cut_scen = []
    TSRFLP.master()
    TSRFLP.sub(callback =1)
    TSRFLP.master_model.params.OutputFlag = 1
    TSRFLP.sub_model.params.OutputFlag = 0
    TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
    TSRFLP.master_model.Params.lazyConstraints = 1
    TSRFLP.master_model.optimize(mycallback)

    # while gap >= stop:
    #     if iteration != 0:
    #         TSRFLP.update_master()
    #     TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
    #     TSRFLP.master_model.optimize(mycallback)
    #     if iteration == 0:
    #         TSRFLP.sub()
    #     else:
    #         TSRFLP.update_sub()
    #     TSRFLP.sub_model.optimize()
    #     gap = TSRFLP.gap_calculation()
    #     TSRFLP.update_status()
    #     if gap <= stop:
    #         print('OPTIMAL SOLUTION FOUND !')
    #         print('Optimal Objective Value = ',str(TSRFLP.UB))
    #     iteration += 1
    #     if iteration >= 20:
    #         break
except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
print("--- %s seconds ---" % round((time.time() - start_time),2))
