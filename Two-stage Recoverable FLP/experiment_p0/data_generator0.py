import numpy as np
import math

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
            self.sum_k = int(round(self.p*0.2))

    def ins(self):
        np.random.seed(self.rnd)
        for i in range(self.ni):
            for j in range(self.ni):
                if i != j:
                    self.cd[i][j] = round(np.random.rand() * 1000)+100
                    for k in range(self.nk):
                        self.cdk[k][i][j] = round(np.random.rand() * 1000)+100
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
