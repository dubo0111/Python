'''
Algorithm camparison
Time Limit = 2000s
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import TSRO_BC_SP_MIP_3_0 as bc # B&C
import TSRO_BC_SP_LB as bc_LB # LB
import TSRO_BC_SP_LB1 as bc_VN # VNSB
import TSRO_BC_SP_LB0 as bc_LBRoot # LB each root nodes
import TSRO_BD_DualSP_INT as bd # BD
import TSRO_LIP as lip # LIP
import data_generator0 as dg0
#import email_self as em

ex_N = [20]  # number of vertexes
ex_k = [50,20,30,40,50] # number of scenarios
ex_all = 10 # number of experiments for each combination

##
a1 = 0.5
rnd_seed = 1 # starting random seed


try:
    result = []
    for n_N in ex_N:
        for n_k in ex_k:
            for n_e in range(ex_all):
                data = dg0.data_gen(n_N, n_k, rnd_seed)
                p, cd, cdk, sk = data.data()
                rnd_seed += 1
                # while rnd_seed in bug:
                #     rnd_seed += 1
                # print(rnd_seed)
                ###
#                y1,t1, cut1, opt1, val1, gap1, conv3 = bc.bra_cut(p, cd, cdk, sk, a1) # branch and cut
#                y10,t10, cut10, opt10, val10, gap10, conv3 = bc_LB.bra_cut(p, cd, cdk, sk, a1) # LB
#                y11,t11, cut11, opt11, val11, gap11,conv3 = bc_LBRoot.bra_cut(p, cd, cdk, sk, a1) # local branching
#                y12,t12, cut12, opt12, val12, gap12,conv3 = bc_VN.bra_cut(p, cd, cdk, sk, a1) # variable neighbourhood branching
#                t2, cut2, opt2, val2, gap2,conv2 = bd.benders_deco(p, cd, cdk, sk, a1)
                y3, t3, opt3, val3, gap3, conv3 = lip.LIP(p, cd, cdk, sk, a1)

                plt.plot(conv3[2],conv3[0],conv3[2],conv3[1])
                plt.show()

                ###
                # result_i = [n_N, n_k, rnd_seed, t1, cut1, opt1, val1, gap1,t2, cut2, opt2, val2, gap2, t3, opt3, val3, gap3]
                # result.append(result_i)
                # result_pd = pd.DataFrame(result, columns=('|N|', '|k|', 'seed', 'BC:time', 'cuts', 'opt', 'objval',
                #                            'gap', 'BD:time', 'cuts', 'opt', 'objval', 'gap', 'LIP:time', 'opt', 'objval', 'gap'))
                # print(result_i)
                # writer = pd.ExcelWriter('output.xlsx')
                # result_pd.to_excel(writer,'Sheet1')
                break
            break
        break
except:
    a=0
#    em.send(0)
else:
    a=0
#    em.send(1)
