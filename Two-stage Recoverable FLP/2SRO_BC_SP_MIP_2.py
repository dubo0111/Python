# Branch and cut
# %reset -f
import model_rflp as mr
#import data_generator1 as dg
#p, cd, cdk, sk = dg.ins_k(20, 100, 40)  # (ni,nk,randomseed*)
import data_generator0 as dg0
data = dg0.data_gen(30,100,21)
p,cd,cdk,sk = data.data()
from gurobipy import *
import time
# Number of nodes
ni = len(cd)
nk = len(cdk)
# weights of two stages
a1 = 0.5
a2 = 1 - a1
try:
    # if gap doesn't change for * iterations, add fractional cut
    gap_cnt = [-GRB.INFINITY,GRB.INFINITY]
    gap_iter = 0
    def mycallback(model, where):
        # if where == GRB.Callback.MIPNODE:
        #     nodecnt = model.cbGet(GRB.Callback.MIPNODE_NODCNT)
        if where == GRB.Callback.MIPSOL:
            vals = model.cbGetSolution(model._vars)
            TSRFLP.value_y = vals[-2 - ni:-2]
            TSRFLP.value_omega = vals[-1]
            TSRFLP.update_sub_dual(callback=1)
            TSRFLP.sub_dual.optimize()
            TSRFLP.worst_scenario()
            TSRFLP.update_cut()
            model.cbLazy(TSRFLP.omega >= TSRFLP.constr_y)
            # status
            nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
            obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
            solcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
            objbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
            objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
            gap_mipsol = abs(objbst - objbnd)/(1.0 + abs(objbst))
#            x = model.cbGetSolution(model._vars)
            print('**** New solution at node %d, obj %g, sol %d, '
                  'gap = %g ****' % (nodecnt, obj, solcnt, gap_mipsol))

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
#            TSRFLP.update_sub(callback=1)
#            TSRFLP.sub_model.optimize()
#            TSRFLP.worst_scenario(1)
#            TSRFLP.gap_calculation(1)
#            print('---------gap:',TSRFLP.int_gap)
#            if abs(TSRFLP.int_gap) >= 1e-4:
#                TSRFLP.update_integer_cut()
#                # cut incumbent solution
#                model.cbLazy(TSRFLP.omega >= TSRFLP.integer_cut)
    start_time = time.time()
    TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)
    TSRFLP.dual = 1
    TSRFLP.intSP = 1
    # ----------Benders' Decompisition----------
    iteration = 0
    gap = 1
    stop = 1e-5
    TSRFLP.master()
    TSRFLP.dual_sub(callback=1)
    TSRFLP.sub(callback=1)
    TSRFLP.params_tuneup()
    TSRFLP.master_model._lastnode = -GRB.INFINITY
    TSRFLP.master_model._vars = TSRFLP.master_model.getVars()
    TSRFLP.master_model.Params.lazyConstraints = 1
    TSRFLP.master_model.optimize(mycallback)

    TSRFLP.update_sub(callback=0)
    TSRFLP.sub_model.optimize()
    TSRFLP.worst_scenario(1)
    TSRFLP.gap_calculation(0,1)
    if abs(TSRFLP.gap) > 1e-4:
        print('============= Start Integer L-shaped cut ================')
#        x = []
#        for i in TSRFLP.master_model.getVars():
#            x.append(i.x)
        TSRFLP.master_model.reset()
#        for i in range(ni):
#            TSRFLP.master_model.getVars()[-ni-2+i].start = x[-ni-2+i]
        TSRFLP.master_model.optimize(mycallback_int)
    print('Optimal solution found: %g' % TSRFLP.master_model.objVal)
except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
    # TSRFLP.error_check()
print("--- %s seconds ---" % round((time.time() - start_time), 2))
