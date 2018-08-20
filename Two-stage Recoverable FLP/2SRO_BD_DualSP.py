# test class
#%reset -f
import model_rflp as mr
#import data_generator1 as dg
#p, cd, cdk, sk = dg.ins_k_alldiff(20,30)
import data_generator0 as dg0
data = dg0.data_gen(3,1,1)
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
#dual sub problem
TSRFLP.dual = 1
TSRFLP.params_tuneup()
# ----------Benders' Decompisition----------
iteration = 0
gap = 1
stop = 1e-5
TSRFLP.master()
TSRFLP.master_model.params.OutputFlag = 0
while abs(TSRFLP.gap) >= stop:
    if iteration != 0:
        TSRFLP.update_master()
    #filename = ''.join(['.\model\master(',str(iteration),').lp'])
    # TSRFLP.master_model.write(filename)
    TSRFLP.master_model.optimize()
    L0 = TSRFLP.master_model.getVarByName('L').x
    print('L0: ',L0)
    if iteration == 0:
        TSRFLP.dual_sub()
    else:
        TSRFLP.update_sub_dual()
    #filename = ''.join(['.\model\sub(',str(iteration),').lp'])
    # TSRFLP.sub_model.write(filename)
    TSRFLP.sub_dual.params.OutputFlag = 0
    TSRFLP.sub_dual.optimize()
    TSRFLP.gap_calculation()
    TSRFLP.update_status()
    TSRFLP.error_check()
    if abs(TSRFLP.gap) <= stop:
        print('OPTIMAL SOLUTION FOUND !')
        print('Optimal Objective Value = ', str(TSRFLP.UB))
    iteration += 1
    if iteration >= 20:
        break
print("--- %s seconds ---" % round((time.time() - start_time), 2))
