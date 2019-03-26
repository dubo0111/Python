'''
Algorithm camparison

 '''
import numpy as np
import pandas as pd
import pickle
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)
import matplotlib.pyplot as plt
import TSRO_BC_SP_MIP_3_0 as bc  # B&C
# import TSRO_BC_SP_LB as bc_LB # LB
import TSRO_BC_SP_LB1 as bc_VN  # LB
import TSRO_BC_SP_LB2 as bc_PR  # hybrid of LB & PR with time Limits
import TSRO_BC_SP_LB3 as bc_Hybrid  # hybrid of LB & PR when no improvement
import TSRO_BD_DualSP_INT as bd  # BD
import TSRO_LIP as lip  # LIP
import data_generator0 as dg0
import primal_integral as pi
# import email_self as em

# experiments parameters
ex_N = [10, 20, 30, 40, 50]  # number of vertexes
#ex_k = [10,20,30,40,50,100,200,500]
ex_k = [10, 20, 50, 100, 200, 500]  # number of scenarios
ex_all = 5  # number of experiments for each combination
# problem parameters
Time_Limit = 500
a1 = 0.5
rnd_seed = 100 # starting random seed
# algorithm parameters
# LB
tl_node = 60
tl_total = 300
# PR solo (13)
tl_node1 = 100
tl_total1 = 0
tl_pr_node1 = 1000
tl_pr_total1 = 1000
# LB+PR (14) : abort when no improve
tl_node2 = 60
tl_total2 = 200
tl_pr_node2 = 1000
tl_pr_total2 = 200
# REMARK: tl_pr_total
# PR parameters
branch_step = 2
pr_gap = 0.05
pr_terminate = 1e5  # hard
# pr_terminate = 1e6 # soft
pr_step = [2 / 3, 1 / 2]  # hard,soft

timeline = [i for i in range(10,510,10)]

try:
    result = []
    result_sum = []
    pi_final = []
    bb = np.empty((0, 26), float)
    for n_N in ex_N:
        for n_k in ex_k:
#            if n_N == 20:
#                tl_total = n_N*n_k/100+5 #
#            else:
#                tl_total = 200
            # tl_total1 = n_N*n_k/100+5
            # tl_pr_total = n_N*n_k/100+5
            for n_e in range(ex_all):
                p = round(n_N / 3)
                sumk = int(round(p * 0.5))
#                branch_step = round(p/2) #
                data = dg0.data_gen(n_N, n_k, sumk, rnd_seed)
#                data = dg0.data_gen(n_N, n_k, sumk)
                _, cd, cdk, sk = data.data()
                rnd_seed += 1
                [t1, cut1, gap1, t12, cut12, gap12, t13,
                    cut13, gap13, t14, cut14, gap14, rootval] = [0 for i in range(13)]
                Heu_sol1 = [0, 0]
                Heu_sol12 = [0, 0, 0]
                Heu_sol13 = [0, 0, 0]
                Heu_sol14 = [0, 0, 0]
                conv1=[]

                if n_N >= 40:
                    if n_k == 100 and n_e==4:
                        y1,t1, cut1, opt1, val1, gap1, conv1, Heu_sol1 = bc.bra_cut(
                             Time_Limit,p, cd, cdk, sk, a1) # branch and cut
                        y12, t12, cut12, opt12, val12, gap12, conv12, pool12, Heu_sol12 = bc_VN.bra_cut(
                             Time_Limit, p, cd, cdk, sk, a1, tl_total, tl_node, branch_step)
                        y13, t13, cut13, opt13, val13, gap13, conv13, pool13, Heu_sol13, rootval = bc_PR.bra_cut(
                             Time_Limit, p, cd, cdk, sk, a1, tl_total1, tl_node1, tl_pr_node1, tl_pr_total1, branch_step, pr_gap, pr_terminate, pr_step)
                        y14, t14, cut14, opt14, val14, gap14, conv14, pool14, Heu_sol14, rootval = bc_Hybrid.bra_cut(
                             Time_Limit, p, cd, cdk, sk, a1, tl_total2, tl_node2, tl_pr_node2, tl_pr_total2, branch_step, pr_gap, pr_terminate, pr_step)
                if n_N >= 40:
                    if n_k >=200:
                        y1,t1, cut1, opt1, val1, gap1, conv1, Heu_sol1 = bc.bra_cut(
                             Time_Limit,p, cd, cdk, sk, a1) # branch and cut
                        y12, t12, cut12, opt12, val12, gap12, conv12, pool12, Heu_sol12 = bc_VN.bra_cut(
                             Time_Limit, p, cd, cdk, sk, a1, tl_total, tl_node, branch_step)
                        y13, t13, cut13, opt13, val13, gap13, conv13, pool13, Heu_sol13, rootval = bc_PR.bra_cut(
                             Time_Limit, p, cd, cdk, sk, a1, tl_total1, tl_node1, tl_pr_node1, tl_pr_total1, branch_step, pr_gap, pr_terminate, pr_step)
                        y14, t14, cut14, opt14, val14, gap14, conv14, pool14, Heu_sol14, rootval = bc_Hybrid.bra_cut(
                             Time_Limit, p, cd, cdk, sk, a1, tl_total2, tl_node2, tl_pr_node2, tl_pr_total2, branch_step, pr_gap, pr_terminate, pr_step)
                if n_N >= 50:
                    y1,t1, cut1, opt1, val1, gap1, conv1, Heu_sol1 = bc.bra_cut(
                         Time_Limit,p, cd, cdk, sk, a1) # branch and cut
                    y12, t12, cut12, opt12, val12, gap12, conv12, pool12, Heu_sol12 = bc_VN.bra_cut(
                         Time_Limit, p, cd, cdk, sk, a1, tl_total, tl_node, branch_step)
                    y13, t13, cut13, opt13, val13, gap13, conv13, pool13, Heu_sol13, rootval = bc_PR.bra_cut(
                         Time_Limit, p, cd, cdk, sk, a1, tl_total1, tl_node1, tl_pr_node1, tl_pr_total1, branch_step, pr_gap, pr_terminate, pr_step)
                    y14, t14, cut14, opt14, val14, gap14, conv14, pool14, Heu_sol14, rootval = bc_Hybrid.bra_cut(
                         Time_Limit, p, cd, cdk, sk, a1, tl_total2, tl_node2, tl_pr_node2, tl_pr_total2, branch_step, pr_gap, pr_terminate, pr_step)

                                                # compute Primal Integral
                    # plt.plot(conv1[2],conv1[0],conv1[2],conv1[1])
                    # plt.plot(conv12[2], conv12[0], conv12[2], conv12[1])
                    # plt.plot(conv13[2], conv13[0], conv13[2], conv13[1])
                    # plt.plot(conv14[2], conv14[0], conv14[2], conv14[1])
                    # plt.show()
                if conv1 != []:
                    pi_all,best=pi.primal_integral(conv1,conv12,conv13,conv14,timeline)
                    pi_fname = 'data/'+str(n_N)+'_'+str(n_k)+'_'+str(n_e)+'_'+str(rnd_seed)
                    conv_name = 'data/'+str(n_N)+'_'+str(n_k)+'_'+str(n_e)+'_'+str(rnd_seed)+'conver'
                    with open(pi_fname, 'wb') as f:
                       pickle.dump(pi_all, f)
                    pi_final.append(pi_all)
                    with open(conv_name, 'wb') as f1:
                       pickle.dump([conv1,conv12,conv13,conv14], f1)

                result_i = [[n_N, n_k, t1, cut1, gap1, Heu_sol1[0], Heu_sol1[1], t12, cut12, gap12,
                             Heu_sol12[0], Heu_sol12[1], Heu_sol12[2], t13, cut13, gap13, Heu_sol13[0], Heu_sol13[1], Heu_sol13[2],
                             t14, cut14, gap14, Heu_sol14[0], Heu_sol14[1], Heu_sol14[2], rootval]]
                result.append(result_i[0])
                result_sum.append(result_i[0])
                result_pd = pd.DataFrame(result, columns=('|N|', '|k|',
                                                  'BC:time', 'cuts', 'gap', 'objval', 'time',
                                                  'LB:time', 'cuts', 'gap', 'Heu_sol', 'time', 'opt_gap',
                                                  'PR:time', 'cuts', 'gap', 'Heu_sol', 'time', 'opt_gap',
                                                  'HY:time', 'cuts', 'gap', 'Heu_sol', 'time', 'opt_gap',
                                                  'ROOT_val'))
                print(result_i)
                writer = pd.ExcelWriter('output.xlsx')
                result_pd.to_excel(writer, 'Sheet1')
                writer.save()
        #         break
        #     break
        # break
            aa = np.array(result_sum).sum(axis=0) / ex_all
            bb = np.vstack((bb, aa))
            result_sum = []
            Final = pd.DataFrame(bb, columns=('|N|', '|k|',
                                              'BC:time', 'cuts', 'gap', 'objval', 'time',
                                              'LB:time', 'cuts', 'gap', 'Heu_sol', 'time', 'opt_gap',
                                              'PR:time', 'cuts', 'gap', 'Heu_sol', 'time', 'opt_gap',
                                              'HY:time', 'cuts', 'gap', 'Heu_sol', 'time', 'opt_gap',
                                              'ROOT_val'))
            writer1 = pd.ExcelWriter('output_final.xlsx')
            Final.to_excel(writer1, 'Sheet1')
            writer1.save()
    # geometric mean
    pi_report = [[0 for i in range(len(timeline))] for j in range(4)]
    for i in range(4): #算法
        for t in range(len(timeline)): # 时间点
            temp = []
            for n in range(len(pi_final)): # 实验
                if pi_final[n][i][t] < 1:
                    pi_final[n][i][t] = 1 ##
                temp.append(pi_final[n][i][t])
            a = np.log(temp)
            pi_report[i][t] = np.exp(a.sum()/len(a))
    df = pd.DataFrame(pi_report)
    df = df.transpose()
    writer2 = pd.ExcelWriter('PI.xlsx')
    df.to_excel(writer2, 'Sheet1')
    writer2.save()

    plt.plot(timeline,pi_report[0],
             timeline,pi_report[1],
             timeline,pi_report[2],
             timeline,pi_report[3])
    plt.show()


except:
    a = 0
    # em.send(0)
else:
    a = 0
    # em.send(1)
