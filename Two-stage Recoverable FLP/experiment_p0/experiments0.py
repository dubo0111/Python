'''
Algorithm camparison
Time Limit = 2000s
'''
import numpy as np
import TSRO_BC_SP_MIP_3 as bc
import TSRO_BD_DualSP_INT as bd
import TSRO_LIP as lip
import data_generator0 as dg0
#import email_self as em

ex_N = [10,20,40,60]  # number of vertexes
ex_k = [10,30,60,100] # number of scenarios
ex_all = 100000000 # number of experiments for each combination
#ex_N = [5, 10]
#ex_k = [5, 10]
#ex_all = 1
a1 = 0.5
rnd_seed = 15  #starting random seed
##############
# bug record
# 20 20 17
# 10 10 15
##############

result = []
for n_N in ex_N:
    for n_k in ex_k:
        for n_e in range(ex_all):
            data = dg0.data_gen(n_N, n_k, rnd_seed)
            p, cd, cdk, sk = data.data()
            rnd_seed += 1
#            print(rnd_seed)
            y1,isinteger,t1, cut1, opt1, val1, gap1 = bc.bra_cut(p, cd, cdk, sk, a1)
#            t2, cut2, opt2, val2, gap2 = bd.benders_deco(p, cd, cdk, sk, a1)
            y3,t3, opt3, val3, gap3 = lip.LIP(p, cd, cdk, sk, a1)
            
            break
        break
    break
#            if isinteger == 1 and abs(val1-val3)>1e-5:
#                break
#        if isinteger == 1 and abs(val1-val3)>1e-5:
#            break
#    if isinteger == 1 and abs(val1-val3)>1e-5:
#        break

y_diff = []
for i in range(len(y1)):
    y_diff.append(y1[i]-y3[i])
print(y_diff)
                    

