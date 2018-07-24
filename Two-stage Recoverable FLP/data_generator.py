# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 15:16:34 2018

@author: DUBO

Instance Generator + US-250city dataset
"""
import scipy.io
import numpy as np
import math
# read .mat file
mat = scipy.io.loadmat('city_250.mat')
a = mat['city_250']
# coordinate of each city
coordinate = a[:,1:3]
cxy = np.vstack(coordinate).astype(np.float)
#demand of each city
demand = a[:,3]
demand = np.vstack(demand).astype(np.float)
#function for calculating distance and total cost
def cal_cost(coor,d):
    ni = coor.shape[0]
    cost = np.zeros((ni,ni))
    total_cost = cost
    for i in range(ni):
        for j in range(ni):
           cost[i,j] = math.sqrt((cxy[i,0]-cxy[j,0])**2 + (cxy[i,1]-cxy[j,1])**2)
           total_cost[i,j] = cost[i,j]*d[i,0]
    return total_cost
#calculate 250city cost matrix
#total_cost = cal_cost(cxy,demand)

#benchmark instance generation
#number of nodes, scenarios and maximum disrupted facilities
def ins_generator(ni,nk,sumk): 
    np.arrange
    
#small instance

#experiment instances
    
    
    
    
    
    
    