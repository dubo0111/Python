"""
Created on Thu Feb 28 2019

@author: DUBO

* Experiments pqcenter

"""
import pqcenter as pc
import pandas as pd
import numpy as np

# TEST
#for n in range(100):
#    print("*************************n = ",n)
# pc.pqcenter_compare(100,25) #17 (20,12)


# Experiments 1:
result = []
result_2 = []
bb = np.empty((0,8), float)
rnd = 0
size = [20,50,100,200]
for ni in size:
    num_p = [round(ni/5),round(ni/4),round(ni/3)]
    for p in num_p:
        for n in range(10):
            print('------------------ ni = ',ni,'; p = ',p,'; rnd = ',rnd,' ------------------')
            obj0,T0,Gap0,obj1,T1,Gap1 = pc.pqcenter_compare(ni,p,rnd)
            rnd += 1
            result_i = [ni,p,obj0,T0,Gap0,obj1,T1,Gap1]
            result.append(result_i)
            result_2.append(result_i)
            # sum_cols = [sum(x) for x in zip(*result)]
            result_pd = pd.DataFrame(result, columns=('|I|', 'p', 'obj0', 'time0', 'gap0', 'obj1', 'time1', 'gap1'))
            print(result_i)
            # writer = pd.ExcelWriter('output.xlsx')
            # result_pd.to_excel(writer,'Sheet1')
        aa = np.array(result_2).sum(axis=0)/10
        bb = np.vstack((bb, aa)) # restore final output
        result_2 = []
        Final = pd.DataFrame(bb, columns=('|I|', 'p', 'obj0', 'time0', 'gap0', 'obj1', 'time1', 'gap1'))
        writer1 = pd.ExcelWriter('output_final.xlsx')
        Final.to_excel(writer1,'Sheet1')
        print(Final)



#
