'''
Algorithm camparison

 '''
import numpy as np
import pandas as pd
import pickle
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)
import matplotlib.pyplot as plt
import primal_integral as pi

# experiments parameters
ex_N = [10, 20, 30, 40, 50]  # number of vertexes
#ex_k = [10,20,30,40,50,100,200,500]
ex_k = [10, 20, 50, 100, 200, 500]  # number of scenarios
ex_all = 5  # number of experiments for each combination
rnd_seed = 100
timeline = [i for i in range(10,510,10)]

conv_all1 = []
conv_all12 = []
conv_all13 = []
conv_all14 = []

pi_final = []

for n_N in ex_N:
    for n_k in ex_k:
        for n_e in range(ex_all):
            rnd_seed += 1
            if n_N == 40:
                # pi_fname = 'data/'+str(n_N)+'_'+str(n_k)+'_'+str(n_e)+'_'+str(rnd_seed)
                conv_name = 'data/'+str(n_N)+'_'+str(n_k)+'_'+str(n_e)+'_'+str(rnd_seed)+'conver'
                with open(conv_name, 'rb') as f:
                   [conv1,conv12,conv13,conv14] = pickle.load(f)
                conv_all1.append(conv1)
                conv_all12.append(conv12)
                conv_all13.append(conv13)
                conv_all14.append(conv14)
                pi_all,best=pi.primal_integral(conv1,conv12,conv13,conv14,timeline)
                pi_final.append(pi_all)

#                a1, = plt.plot(timeline,pi_all[0],'-',label="BC")
#                a2, = plt.plot(timeline,pi_all[1],'--',label="LB")
#                a3, = plt.plot(timeline,pi_all[2],'-.',label="PR")
#                a4, = plt.plot(timeline,pi_all[3],':',label="LB+PR")
#                plt.legend(handles=[a1,a2,a3,a4])
#                plt.show()
#
#                b1, = plt.plot(conv1[2],conv1[0],'-',label="BC")
#                b2, = plt.plot(conv12[2],conv12[0],'--',label="LB")
#                b3, = plt.plot(conv13[2],conv13[0],'-.',label="PR")
#                b4, = plt.plot(conv14[2],conv14[0],':',label="LB+PR")
#                plt.legend(handles=[b1,b2,b3,b4])
#                plt.show()

#                aaa = 1


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
# plt.axis([0, 500, 0, 35])
#plt.title('|N| = 40')
plt.xlabel('Time(s)')
plt.ylabel('Primal integral')

plt.show()
