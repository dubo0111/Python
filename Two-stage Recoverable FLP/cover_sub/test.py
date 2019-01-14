#a= 'love'
#b ='i love you'
#a in b
#
#if any(name in b for name in a):
#    print(1)
#a = [5,3,4,2]
#b = 4
#num = [[0 for i in range(a[j])] for j in range(b)]
#print(num)

from gurobipy import *
import numpy as np
#
#m = Model("p-center-cover")
#ne = 5
#nk = 4
## Create variables
## z:ordered cost, y:location
#z = m.addVars(ne,vtype=GRB.BINARY, name="z")
#
##z=[[[] for e in range(ne)] for k in range(nk)]
##for k in range(nk):
###    for e in range(ne):
##    z[k] = m.addVar(vtype=GRB.BINARY, name="z")
#
#m.update()
#print(m.getVars())

test = np.array([[[1,2], 2], [3, 4], [5, 6]])
#test[:,0] = 0
print(test[0,0])