'''
Benders' Decomposition: old fashion
Du Bo
'''
import model_rflp as mr
import data_generator0 as dg0
from gurobipy import *
import time

def benders_deco(p,cd,cdk,sk,a1):
    convergence = []
    # Number of nodes
    ni = len(cd)
    nk = len(cdk)
    # weights of two stages
    a2 = 1 - a1
    start_time = time.time()
    TSRFLP = mr.rflp(p, ni, nk, a1, a2, cd, cdk, sk)
    #dual sub problem
    TSRFLP.dual = 1 # for dual sub problem
    TSRFLP.params_tuneup()
    TSRFLP.intSP = 1 #!!BINARY subproblem variable
    # ----------Benders' Decompisition----------
    iteration = 0
    stop = 1e-5
    TSRFLP.master()
    TSRFLP.master_model.params.OutputFlag = 0
    while abs(TSRFLP.gap) >= stop:
        if iteration != 0:
            TSRFLP.update_master()
        time_master= time.time()
        TSRFLP.master_model.optimize()
        if iteration == 0:
            TSRFLP.dual_sub()
        else:
            TSRFLP.update_sub_dual()
        TSRFLP.sub_dual.params.OutputFlag = 0
        time_sub = time.time()
        TSRFLP.sub_dual.optimize()
        TSRFLP.gap_calculation()
        TSRFLP.update_status()
        iteration += 1
        runtime = time.time() - start_time
        if convergence != []:
            convergence.append([convergence[-1][0],convergence[-1][1],time.time()-start_time])
        convergence.append([TSRFLP.UB,TSRFLP.LB,runtime])
        time_limit = 0
        if runtime >= 1000:
            time_limit = 1
            break
    # ------------Integer cut --------------
    if time_limit != 1:
        # Check optimality with primal SP first
        TSRFLP.UB = GRB.INFINITY #reset
        TSRFLP.sub()
        TSRFLP.sub_model.params.OutputFlag = 0
        TSRFLP.sub_model.optimize()
        TSRFLP.gap_calculation(0,1)
        if TSRFLP.gap > 1e-4:
            TSRFLP.update_integer_cut()
            TSRFLP.master_model.addConstr(TSRFLP.omega >= TSRFLP.integer_cut)
        # print('========= Start Integer L-shaped cut ==========')

    while abs(TSRFLP.gap) > 1e-4 and time_limit != 1:
        # Integer L-shaped cut
        TSRFLP.update_integer_cut()
        TSRFLP.master_model.addConstr(TSRFLP.omega >= TSRFLP.integer_cut)
        # Benders' cut
        TSRFLP.update_sub_dual()
        TSRFLP.sub_dual.optimize()
        TSRFLP.update_cut()
        TSRFLP.master_model.addConstr(TSRFLP.omega >= TSRFLP.constr_y)
        TSRFLP.master_model.optimize()
        TSRFLP.update_sub(callback=0)
        TSRFLP.sub_model.optimize()
        TSRFLP.gap_calculation(0,1) # calculate gap
        TSRFLP.update_status() #
        iteration += 1
        runtime = time.time() - start_time
        convergence.append([convergence[-1][0],convergence[-1][1],time.time()-start_time])
        convergence.append([TSRFLP.UB,TSRFLP.LB,runtime])
        if runtime >= 1000:
            break
    # print('Optimal Objective Value = ', str(TSRFLP.UB))
    # print("--- %s seconds ---" % runtime)
    if TSRFLP.int_gap == float('inf'):
        gap = TSRFLP.gap
    else:
        gap = TSRFLP.int_gap
    if abs(gap) <= 1e-5:
        TSRFLP.opt = 1
        gap = 0
    objval = round(TSRFLP.UB,2)
    convergence = [*zip(*convergence)]
    runtime = round(runtime,2)
    gap = round(gap,2)
    return runtime,iteration,TSRFLP.opt,objval,gap,convergence
