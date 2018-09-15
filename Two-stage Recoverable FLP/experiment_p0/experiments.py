'''
Algorithm camparison
Time Limit = 2000s
'''
import numpy as np
import pandas as pd
import TSRO_BC_SP_MIP_3 as bc
import TSRO_BD_DualSP_INT as bd
import TSRO_LIP as lip
import data_generator0 as dg0
#import email_self as em

ex_N = [20,20,40,60]  # number of vertexes
ex_k = [20,30,60,100] # number of scenarios
ex_all = 10 # number of experiments for each combination
#ex_N = [5, 10]
#ex_k = [5, 10]
#ex_all = 1
a1 = 0.5
rnd_seed = 17  # starting random seed

try:
    result = []
    for n_N in ex_N:
        for n_k in ex_k:
            for n_e in range(ex_all):
                data = dg0.data_gen(n_N, n_k, rnd_seed)
                p, cd, cdk, sk = data.data()
                rnd_seed += 1
#                t1, cut1, opt1, val1, gap1 = bc.bra_cut(p, cd, cdk, sk, a1)
                t2, cut2, opt2, val2, gap2 = bd.benders_deco(p, cd, cdk, sk, a1)
#                t3, opt3, val3, gap3 = lip.LIP(p, cd, cdk, sk, a1)
#                result_i = [n_N, n_k, rnd_seed, t1, cut1, opt1, val1, gap1,
#                            t2, cut2, opt2, val2, gap2, t3, opt3, val3, gap3]
#                result.append(result_i)
#                result_pd = pd.DataFrame(result, columns=('|N|', '|k|', 'seed', 'BC:time', 'cuts', 'opt', 'objval',
#                                           'gap', 'BD:time', 'cuts', 'opt', 'objval', 'gap', 'LIP:time', 'opt', 'objval', 'gap'))
#                print(result_pd)
#                writer = pd.ExcelWriter('output.xlsx')
#                result_pd.to_excel(writer,'Sheet1')
                break
            break
        break
except:
    a=0
#    em.send(0)
else:
    a=0
#    em.send(1)
