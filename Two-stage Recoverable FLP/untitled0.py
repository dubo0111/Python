# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 02:25:01 2018

@author: DUBO
"""
import data_generator1 as dg
import random
#import scipy.io
#import numpy as np
#import math
#p, cd, cdk, sk = dg.ins_k_alldiff(3,3)
random.seed(1)
co1,d1 = dg.ins_generator(3,1)
c1,cd1 = dg.cal_cost(co1,d1)
random.seed(2)
co2,d2 = dg.ins_generator(3,2)
c2,cd2 = dg.cal_cost(co2,d2)