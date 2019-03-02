import numpy as np
import math
import copy

class data_gen:
    # INPUT:
    ni = 0
    nk = 0
    rnd = None
    # Attributes
    cd = []
    cdk = []
    sk = []
    p = 0
    sum_k = 0

    def __init__(self, ni, nk, rnd=None):
        self.ni = ni
        self.nk = nk
        self.rnd = rnd
        self.cd = [[0 for j in range(ni)] for j in range(ni)]
        self.cdk = [[[0 for j in range(ni)]
                     for i in range(ni)] for k in range(nk)]
        self.sk = [[0 for i in range(ni)] for k in range(nk)]
        self.p = round(ni / 3)
        if self.p == 1:
            self.sum_k = 1
        else:
            self.sum_k = int(round(self.p/2))

    def ins(self):
        np.random.seed(self.rnd)
        for i in range(self.ni):
            for j in range(self.ni):
                if i != j:
                    self.cd[i][j] = round(np.random.rand() * 100)#+100
                    for k in range(self.nk):
                        self.cdk[k][i][j] = round(np.random.rand() * 100)#+100
                        # self.cdk[k][i][j] = self.cd[i][j]

    def ins_Kdisrupt(self):
        # np.random.seed(self.rnd)
        for k in range(self.nk):
            a = np.random.choice(range(self.ni), self.sum_k, replace=False)
            for i in a:
                self.sk[k][i] = 1

    def data(self):
        self.ins()
        self.ins_Kdisrupt()
        return self.p, self.cd, self.cdk, self.sk
    #
    def illustrative(self):
        np.random.seed(self.rnd)
        coordinate = [[0,0] for n in range(self.ni)]
        for x in coordinate:
            x[0] = round(np.random.rand() * 100)
            x[1] = round(np.random.rand() * 100)
        for i in range(self.ni):
            for j in range(self.ni):
                self.cd[i][j]=math.sqrt((abs(coordinate[i][0]-coordinate[j][0]))**2 + (abs(coordinate[i][1]-coordinate[j][1]))**2)

        for k in range(self.ni):
            self.cdk[k] = copy.deepcopy(self.cd) # 不考虑需求不确定
        self.sk = [[0 for j in range(self.ni)] for i in range(self.ni)] #
        for i in range(self.ni):
            self.sk [i][i] = 1
        return  self.cd, self.cdk, self.sk, coordinate









        #
