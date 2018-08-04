# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 15:16:34 2018

@author: DUBO

Instance Generator + US-250city dataset
"""
import scipy.io
import numpy as np
import math
import random
#import pandas as pd
#from pandas import ExcelWriter
#from pandas import ExcelFile
#read .mat file
mat = scipy.io.loadmat('city_250.mat')
a = mat['city_250']
# coordinate of each city
coordinate = a[:,1:3]
cxy = np.vstack(coordinate).astype(np.float)
#demand of each city
demand = a[:,3]
demand = np.vstack(demand).astype(np.float)

# small instance
def ins_small(nk=0):
    p = 1
    C = np.matrix([[0,20,50],[30,0,40],[30,50,0]])
    D = np.matrix([100,200,150])
    cd = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            cd[i,j] = C[i,j]*D[0,i]
    if nk != 0:
        print('nk = ',nk)
    cd = cd.tolist()
    return p,cd
#p,cd,c00 = ins_small()

#function for calculating distance and total cost(CD)
def cal_cost(coor,d):
    ni = coor.shape[0]
    cost = np.zeros((ni,ni))
    total_cost = cost.copy()
    for i in range(ni):
        for j in range(ni):
           cost[i,j] = math.sqrt((cxy[i,0]-cxy[j,0])**2 + (cxy[i,1]-cxy[j,1])**2)
           total_cost[i,j] = cost[i,j]*d[i,0]
    return cost,total_cost
#calculate 250city cost matrix
#total_cost = cal_cost(cxy,demand)

# instance random generation from 250 city
#ni, nk, sumk: number of nodes, scenarios and maximum disrupted facilities
def ins_generator(ni,nk=0,sumk=0):
    list250 = list(range(0,250))
    random.seed(3)
    city = random.sample(list250,ni)
    city.sort()
    newco = np.zeros((ni,2))
    newdemand = np.zeros((ni,1))
    for x in range(ni):
        x1 = city[x]
        newco[x,:]=cxy[x1,:]
        newdemand[x,0]=cxy[x1,0]
    if nk != 0:
        print('Generating',nk,'senarios......')
    return newco,newdemand

# big instances
def ins_big(ni):
    aa,bb = ins_generator(ni)
    cd = cal_cost(aa,bb)[1]
    cd = cd.tolist()
    p = round(ni*1/3)
    return p,cd

# scenarios: disrupted facilities (assumption: sumk = p-1)
def ins_kdisrupt(ni,nk,sumk):
    random.seed(2)
    b = [[0 for i in range(ni)] for k in range(nk)]
    for k in range(nk):
        a = random.sample(range(ni),sumk)
        a.sort()
        for i in a:
            b[k][i] = 1
    return b
#sk = ins_kdisrupt(5,3,2)
# instance with |k| scenarios
def ins_k(ni,nk):
    co,d = ins_generator(ni)
    c,cd = cal_cost(co,d)
    cdk = [[[0 for j in range(ni)] for i in range(ni)] for k in range(nk)]
    np.random.seed(2)
    for k in range(nk):
        randc2 = np.random.rand(ni,ni)+0.5
        randd2 = np.random.rand(ni,1)+0.5
        c2 = c*randc2
        d2 = d*randd2
        cd2 = np.zeros((ni,ni))
        for i in range(ni):
            for j in range(ni):
                cd2[i,j] = c2[i,j]*d2[i,0]
        cd2.tolist()
        cdk[k] = cd2
        cdk[k] = cdk[k].tolist()
        p = round(ni*1/3)
        if p == 1:
            sumk = 1
        else:
            sumk = p - 1
        sk = ins_kdisrupt(ni,nk,sumk)
    return p,cd,cdk,sk

#ins_k(5,3)
