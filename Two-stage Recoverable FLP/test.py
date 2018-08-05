# -*- coding= utf-8 -*-
"""
Created on Fri Jul 27 09=00=30 2018

@author= DUBO
"""
ni = 3
beta=[[0 for i in range(ni)]]
gamma=[[[0 for i in range(ni)] for j in range(ni)]]
delta=[[[0 for i in range(ni)] for j in range(ni)]]
epsilon=[[0 for i in range(ni)]]
lamda = [[0 for j in range(ni)]]
mu = [[0 for i in range(ni)]]
nu = [0]

beta[0][0] = -0.3044989171426794
beta[0][2] = -0.6955010828573206
gamma[0][0][0] = -842.3231484589844
gamma[0][2][1] = -842.3231484589844
delta[0][0][2] = -283.03673777616837
delta[0][2][2] = -1304.5734362289518
epsilon[0][0] = 842.3231484589844
epsilon[0][2] = 1304.5734362289518
nu[0] = -842.3231484589844
####

q = - gamma[0][0][2] - gamma[0][1][2] - gamma[0][2][2] - delta[0][0][0]\
   - delta[0][0][1] - delta[0][1][0] - delta[0][1][1] - delta[0][2][0]\
   - delta[0][2][1] - epsilon[0][0] - epsilon[0][1] - epsilon[0][2]\
   - lamda[0][0] - lamda[0][1] - mu[0][0] - mu[0][1] - nu[0]

u=[[0 for i in range(ni)]]
v=[[[0 for i in range(ni)] for j in range(ni)]]
L3=[0]

u[0][0]= - gamma[0][0][0] - gamma[0][1][0] - gamma[0][2][0] + lamda[0][0]+ mu[0][0] + nu[0] <= 0
u[0][1]= - gamma[0][0][1] - gamma[0][1][1] - gamma[0][2][1] + lamda[0][1]+ mu[0][1] + nu[0] <= 0
u[0][2]= - gamma[0][0][2] - gamma[0][1][2] - gamma[0][2][2] + lamda[0][2]+ mu[0][2] + nu[0] <= 0
v[0][0][0]= gamma[0][0][0] + delta[0][0][0] + epsilon[0][0] <= 0
v[0][0][1]= 2766.259914363821* beta[0][0] + gamma[0][0][1] + delta[0][0][1]+ epsilon[0][0] <= 0
v[0][0][2]= 1836.743512689572* beta[0][0] + gamma[0][0][2] + delta[0][0][2] + epsilon[0][0] <= 0
v[0][1][0]= 1789.673772413146* beta[0][1] + gamma[0][1][0] + delta[0][1][0] + epsilon[0][1] <= 0
v[0][1][1]= gamma[0][1][1] + delta[0][1][1] + epsilon[0][1] <= 0
v[0][1][2]= 987.4596306348664* beta[0][1] + gamma[0][1][2] + delta[0][1][2] + epsilon[0][1] <= 0
v[0][2][0]= 1875.731711113065* beta[0][2] + gamma[0][2][0] + delta[0][2][0] + epsilon[0][2] <= 0
v[0][2][1]= 664.6291417274418* beta[0][2] + gamma[0][2][1] + delta[0][2][1]+ epsilon[0][2] <= 0
v[0][2][2]= gamma[0][2][2] + delta[0][2][2] + epsilon[0][2] <= 0
L3[0]= - beta[0][0] - beta[0][1] - beta[0][2] <= 1
