a = [2,3,5,1,7]
index_min = min(range(len(a)), key=a.__getitem__)

b = [2,3,4,5,6]

c= [x for x in b if x <=4 ]


#plt.plot(timeline,pi_all[0],timeline,pi_all[1],timeline,pi_all[2],timeline,pi_all[3])

table = [[[1 , 2], [3, 4]],[[1 , 2], [3, 4]]]
# df = DataFrame(table)
# df = df.transpose()
# df.columns = ['Heading1', 'Heading2']

# import pickle
# with open('data/file', 'wb') as f:
#    pickle.dump(table, f)
# with open('data/file', 'rb') as f:
#     table1 = pickle.load(f)
import numpy as np
temp = [6,12,36,24]
#temp= [0.14,0.17]
a = np.log(temp)
b = np.exp(a.sum()/len(a))
