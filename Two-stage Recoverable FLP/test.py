# test class
import model_rflp as mr
import data_generator1 as dg
p,cd,cdk,sk = dg.ins_k(5,100) #(ni,nk,randomseed)
from gurobipy import *
# Number of nodes
ni = len(cd)
nk = len(cdk)
# weights of two stages
a1 = 0.5
a2 = 1 - a1

TSRFLP = mr.rflp(p,ni,nk,a1,a2,cd,cdk,sk)
# ----------Benders' Decompisition----------
iteration = 0
gap = 1
stop = 1e-5
add_cut_scen = []
TSRFLP.master()
TSRFLP.master_model.params.OutputFlag = 0
TSRFLP.sub_model.params.OutputFlag = 0
while gap >= stop: # stop criteria
    if iteration != 0:
        TSRFLP.update_master()
    TSRFLP.master_model.optimize()
    if iteration == 0:
        TSRFLP.sub()
    else:
        TSRFLP.update_sub()
    TSRFLP.sub_model.optimize()

    gap = TSRFLP.gap_calculation()
    #
    add_cut_scen.append(TSRFLP.max_k)
    print('==========================================')
    print('Current iteration:',str(iteration))
    print('gap = ',str(gap))
    print('Cuts added form scenario:',str(add_cut_scen[0:-1]))
    if gap <= stop:
        print('OPTIMAL SOLUTION FOUND !')
        print('Optimal Objective Value = ',str(TSRFLP.UB))
    # update iteration
    iteration += 1
    if iteration >= 20:
        break
