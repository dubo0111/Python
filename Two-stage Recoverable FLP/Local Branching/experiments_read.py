'''
Algorithm camparison

 '''
import numpy as np
import pandas as pd
import pickle
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)
import matplotlib.pyplot as plt
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)
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
rnd_seed = 100
timeline = [i for i in range(10,510,10)]

pi_final = []

for n_N in ex_N:
    for n_k in ex_k:
        for n_e in range(ex_all):
            rnd_seed += 1
            if n_N == 40 and n_k != 500:
                pi_fname = 'data/'+str(n_N)+'_'+str(n_k)+'_'+str(n_e)+'_'+str(rnd_seed)
                # conv_name = 'data/'+str(n_N)+'_'+str(n_k)+'_'+str(n_e)+'_'+str(rnd_seed)+'conver'
                with open(pi_fname, 'rb') as f:
                   pi_all = pickle.load(f)
                pi_final.append(pi_all)
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

a1, = plt.plot(timeline,pi_report[0],'k-',label="BC")
a2, = plt.plot(timeline,pi_report[1],'k--',label="LB")
a3, = plt.plot(timeline,pi_report[2],'k-.',label="PR")
a4, = plt.plot(timeline,pi_report[3],'k:',label="LB+PR")
plt.legend(handles=[a1,a2,a3,a4])
plt.axis([0, 500, 0, 35])
#plt.title('|N| = 40')
plt.xlabel('Time(s)')
plt.ylabel('Primal integral')

plt.show()
