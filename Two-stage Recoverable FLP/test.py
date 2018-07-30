# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 09:00:30 2018

@author: DUBO
"""
import random
a= random.sample(range(10),5)
a.sort()
b = [[0 for i in range(10)] for j in range(3)]
for i in range(3):
    for j in a:
        b[j][i] = 1