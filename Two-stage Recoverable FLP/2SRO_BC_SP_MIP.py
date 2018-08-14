 # Branch and cut
# %reset -f
import model_rflp as mr
import data_generator1 as dg
p, cd, cdk, sk = dg.ins_k(100, 100, 40)  # (ni,nk,randomseed*)
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
        # if where == GRB.Callback.MIPNODE:
        #     if model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.Status.OPTIMAL:
        #         vals = model.cbGetNodeRel(model._vars)
        #         TSRFLP.value_y = vals[-2 - ni:-2]
        #         TSRFLP.update_sub_dual(callback=1)
        #         TSRFLP.sub_dual.optimize()
        #         TSRFLP.worst_scenario()
        #         TSRFLP.update_cut()
        #         model.cbCut(TSRFLP.omega >= TSRFLP.constr_y)
        #         print('++++')
        if where == GRB.Callback.MIPSOL:
            vals = model.cbGetSolution(model._vars)
            TSRFLP.value_y = vals[-2 - ni:-2]
            TSRFLP.value_omega = vals[-1]
#            print('omega:',TSRFLP.value_omega)
            TSRFLP.update_sub_dual(callback=1)
            TSRFLP.sub_dual.optimize()
            TSRFLP.worst_scenario()
            TSRFLP.update_cut()
            model.cbLazy(TSRFLP.omega >= TSRFLP.constr_y)
            # status
            nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
            obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
            solcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
            x = model.cbGetSolution(model._vars)
            print('**** New solution at node %d, obj %g, sol %d, '
                  'x[0] = %g ****' % (nodecnt, obj, solcnt, x[0]))
            # integer l-shaped cut
            TSRFLP.update_sub(callback=1)
            TSRFLP.sub_model.optimize()
            TSRFLP.worst_scenario(1)
            TSRFLP.gap_calculation(1)
            if abs(TSRFLP.int_gap) >= 1e-4:
                TSRFLP.update_integer_cut()
                # cut incumbent solution
                model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
    start_time = time.time()
    TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)
    TSRFLP.dual = 1
    TSRFLP.intSP = 1
    TSRFLP.params_tuneup()
    # ----------Benders' Decompisition----------
    iteration = 0
    gap = 1
    stop = 1e-5
    TSRFLP.master()
    TSRFLP.dual_sub(callback=1)
    TSRFLP.sub(callback=1)
    TSRFLP.master_model.params.OutputFlag = 1
    TSRFLP.sub_dual.params.OutputFlag = 0
    TSRFLP.sub_model.params.OutputFlag = 0
    TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
    TSRFLP.master_model.Params.lazyConstraints = 1
    TSRFLP.master_model.optimize(mycallback)
    print('Optimal solution found: %g' % TSRFLP.master_model.objVal)
except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
    # TSRFLP.error_check()
print("--- %s seconds ---" % round((time.time() - start_time), 2))
