'''
Benders' Decomposition: old fashion
Du Bo
'''
import model_rflp as mr
#import data_generator1 as dg
import data_generator0 as dg0
data = dg0.data_gen(20,20,2)
p,cd,cdk,sk = data.data()

from gurobipy import *
import time
# Number of nodes
ni = len(cd)
nk = len(cdk)
# weights of two stages
a1 = 0.4 
a2 = 1 - a1
start_time = time.time()
TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)
#dual sub problem
TSRFLP.dual = 1
#TSRFLP.params_tuneup()
# ----------Benders' Decompisition----------
iteration = 0
gap = 1
stop = 1e-5
TSRFLP.master()
TSRFLP.master_model.params.OutputFlag = 0
TSRFLP.params_tuneup()
#TSRFLP.master_model.params.Presolve = 0
while abs(TSRFLP.gap) >= stop:
    if iteration != 0:
        TSRFLP.update_master()
    time_master= time.time()
    TSRFLP.master_model.optimize()
    print("---time_master--- %s seconds ---" % round((time.time() - time_master), 2))
    L0 = TSRFLP.master_model.getVarByName('L').x
    print('L0: ',L0)
    if iteration == 0:
        TSRFLP.dual_sub()
    else:
        TSRFLP.update_sub_dual()
    TSRFLP.sub_dual.params.OutputFlag = 0
    time_sub = time.time()
    TSRFLP.sub_dual.optimize()
    print("---time_sub--- %s seconds ---" % round((time.time() - time_sub), 2))
    TSRFLP.gap_calculation()
    TSRFLP.update_status()
    TSRFLP.error_check()
    if abs(TSRFLP.gap) <= stop:
        print('OPTIMAL SOLUTION FOUND !')
        print('Optimal Objective Value = ', str(TSRFLP.UB))
    iteration += 1
    if iteration >= 2000:
        break
#while abs()
print("--- %s seconds ---" % round((time.time() - start_time), 2))
