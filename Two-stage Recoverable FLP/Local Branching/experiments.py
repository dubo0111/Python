'''
Algorithm camparison
Time Limit = 1000s
 '''
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib.pyplot as plt
import TSRO_BC_SP_MIP_3_0 as bc # B&C
import TSRO_BC_SP_LB as bc_LB # LB
import TSRO_BC_SP_LB1 as bc_VN # VNSB
import TSRO_BC_SP_LB2 as bc_PR # LB each root nodes
import TSRO_BD_DualSP_INT as bd # BD
import TSRO_LIP as lip # LIP
import data_generator0 as dg0
# import email_self as em

# experiments parameters
ex_N = [10,20,30]  # number of vertexes
#ex_k = [10,20,30,40,50,100,200,500]
ex_k = [10,20,50,100,200,500] # number of scenarios
ex_all = 5 # number of experiments for each combination
# problem parameters
a1 = 0.5
rnd_seed = 121# starting random seed
# algorithm parameters (LB)
tl_node = 1000
#tl_total = 10
tl_pr_node = 1000
#tl_pr_total = 20
branch_step = 2
pr_gap = 0.05
pr_terminate = 1e5 # hard
# pr_terminate = 1e6 #soft

try:
    result = []
    result_sum = []
    bb = np.empty((0,23),float)
    for n_N in ex_N:
        for n_k in ex_k:
            tl_total = n_N*n_k/100
            tl_pr_total = n_N*n_k/100
            for n_e in range(ex_all):
                p = round(n_N / 3) #
                sumk = int(round(p*0.5)) #
#                branch_step = round(p/2) #
                data = dg0.data_gen(n_N, n_k, sumk, rnd_seed)
#                data = dg0.data_gen(n_N, n_k, sumk)
                _, cd, cdk, sk = data.data()
                rnd_seed += 1
                # while rnd_seed in bug:
                #     rnd_seed += 1
                # print(rnd_seed)
                [ t3, gap3, t2, cut2, gap2, t1, cut1, gap1, t12, cut12, gap12,t13, cut13, gap13,rootval] = [0 for i in range(15)]
                Heu_sol12=[0,0,0]
                Heu_sol13=[0,0,0]

                
                y1,t1, cut1, opt1, val1, gap1, conv1 = bc.bra_cut(
                     p, cd, cdk, sk, a1) # branch and cut
                y12,t12, cut12, opt12, val12, gap12, conv12, pool12,Heu_sol12 = bc_VN.bra_cut(
                         p,cd, cdk, sk, a1, tl_total, tl_node,branch_step) # variable neighbourhood branching
                y13,t13, cut13, opt13, val13, gap13, conv13, pool13,Heu_sol13,rootval = bc_PR.bra_cut(
                     p, cd, cdk, sk, a1, tl_total, tl_node,tl_pr_node,tl_pr_total,branch_step,pr_gap,pr_terminate) # Proximity
                t2, cut2, opt2, val2, gap2,conv2 = bd.benders_deco(
                        p, cd, cdk, sk, a1)
                y3, t3, opt3, val3, gap3, conv3 = lip.LIP(
                        p, cd, cdk, sk, a1)
                
#                plt.plot(conv12[2],conv12[0],conv12[2],conv12[1])
#                plt.show()

                result_i = [[n_N, n_k, t3, gap3, t2, cut2, gap2, t1, cut1, gap1, t12, cut12, gap12,Heu_sol12[0],Heu_sol12[1],Heu_sol12[2],t13, cut13, gap13,Heu_sol13[0],Heu_sol13[1],Heu_sol13[2],rootval]]
                result.append(result_i[0])
                result_sum.append(result_i[0])
                result_pd = pd.DataFrame(result, columns=('|N|', '|k|', 'LIP:time', 'gap', 'BD:time', 'cuts', 'gap',
                                        'BC:time', 'cuts', 'gap', 'LB:time', 'cuts', 'gap','Heu_sol','time', 'opt_gap',
                                        'PR:time', 'cuts', 'gap','Heu_sol','time', 'opt_gap' ,'ROOT_val'))
                print(result_i)
                writer = pd.ExcelWriter('output.xlsx')
                result_pd.to_excel(writer,'Sheet1')
                writer.save()
#                break
#            break
#        break
except:
    a=0
    # em.send(0)
else:
    a=0
    # em.send(1)
