# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 20:24:05 2019

@author: User
"""

import math
import numpy as np

a= [[72,41],[100,87],[74,25],[53,72],[82,26],[23,44]]
#  72.0         41.0
# 100.0         87.0
#  74.0         25.0
#  53.0         72.0
#  90.0         35.0(82,26)
#  23.0         44.0
#result = math.sqrt((abs(53-74))**2 + (abs(72-25))**2)
#print(result)

A = [[0 for i in range(len(a))] for j in range(len(a))]
for i in range(len(a)):
    for j in range(len(a)):
        A[i][j] = round(math.sqrt((abs(a[i][0]-a[j][0]))**2 + (abs(a[i][1]-a[j][1]))**2),2)
A = np.array(A)