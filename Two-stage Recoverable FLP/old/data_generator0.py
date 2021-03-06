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

class dataGen:
    #read .mat file
    mat = scipy.io.loadmat('city_250.mat')
    a = mat['city_250']
    # coordinate of each city
    coordinate = a[:,1:3]
    cxy = np.vstack(coordinate).astype(np.float)
    #demand of each city
    demand = a[:,3]
    demand = np.vstack(demand).astype(np.float)
    def __init__(self,ni,p=1):
        self.ni = ni
        self.p = p
    #function for calculating distance and total cost
    def cal_cost(coordinate):
        ni = coordinate.shape[0]
        cost = np.zeros((ni,ni))
        total_cost = cost
        for i in range(ni):
            for j in range(ni):
               cost[i,j] = math.sqrt((self.cxy[i,0]-self.cxy[j,0])**2 + (self.cxy[i,1]-self.cxy[j,1])**2)
               total_cost[i,j] = cost[i,j]*self.demand[i,0]
        return total_cost
    #calculate 250city cost matrix
    #total_cost = cal_cost(cxy,demand)

    #benchmark instance generation
    #ni, nk, sumk: number of nodes, scenarios and maximum disrupted facilities
    def ins_generator(ni=3,nk=0,sumk=0):
        list250 = list(range(0,250))
        random.seed()
        city = random.sample(list250,ni)
        city.sort()
        newcost = np.zeros((ni,2))
        newdemand = np.zeros((ni,1))
        for x in range(ni):
            x1 = city[x],
            newcost[x,:]=self.cxy[x1,:]
            newdemand[x,0]=self.cxy[x1,0]
        # if nk != 0:
        #     print('Generating',nk,'senarios......')
        # return newcost,newdemand
    #aa,bb = ins_generator(cxy,demand,5)

    # big instances
    def ins_big(ni):
        return

    # small instance
    def ins_small(self.p,nk=0):
        C = np.matrix([[0,20,50],[30,0,40],[30,50,0]])
        D = np.matrix([100,200,150])
        cd = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                cd[i,j] = C[i,j]*D[0,i]
        # if nk != 0:
        #     print('nk = ',nk)
        cd = cd.tolist()
        return p,cd
    #p,cd = ins_small()

    #experiment instances ...to be continue
