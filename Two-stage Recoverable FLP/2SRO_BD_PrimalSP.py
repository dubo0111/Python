# test class
import model_rflp as mr
#import data_generator1 as dg
#p, cd, cdk, sk = dg.ins_k(10, 100)  # (ni,nk,randomseed*)
import data_generator0 as dg0
data = dg0.data_gen(6,30,1)
p,cd,cdk,sk = data.data()

from gurobipy import *
import time
# Number of nodes
ni = len(cd)
nk = len(cdk)
# weights of two stages
a1 = 0.5
a2 = 1 - a1

start_time = time.time()
TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)
# primal sub problem
TSRFLP.dual = 0
TSRFLP.params_tuneup()
# ----------Benders' Decompisition----------
iteration = 0
gap = 1
stop = 1e-5
TSRFLP.master()
TSRFLP.master_model.params.OutputFlag = 0
TSRFLP.sub_model.params.OutputFlag = 0
while TSRFLP.gap >= stop:
    if iteration != 0:
        TSRFLP.update_master()
    #filename = ''.join(['.\model\master(',str(iteration),').lp'])
    # TSRFLP.master_model.write(filename)
    TSRFLP.master_model.optimize()
    if iteration == 0:
        TSRFLP.sub()
    else:
        TSRFLP.update_sub()
    #filename = ''.join(['.\model\sub(',str(iteration),').lp'])
    # TSRFLP.sub_model.write(filename)
    TSRFLP.sub_model.optimize()
    TSRFLP.gap_calculation()
    TSRFLP.update_status()
    TSRFLP.error_check()
    if TSRFLP.gap <= stop:
        print('OPTIMAL SOLUTION FOUND !')
        print('Optimal Objective Value = ', str(TSRFLP.UB))
    iteration += 1
    if iteration >= 50:
        break
print("--- %s seconds ---" % round((time.time() - start_time), 2))
