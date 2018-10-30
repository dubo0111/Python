"""
Created on Tue Sep 27 2018

@author: DUBO

Pickup and delivery
"""
from gurobipy import *

# ================== Input parameter ======================
distance = [[0, 20, 10, 3, 15, 8, 17, 17, 8, 10],
            [20, 0, 24, 19, 12, 12, 7, 13, 18, 21],
            [10, 24, 0, 13, 14, 14, 24, 15, 6, 3],
            [3, 19, 13, 0, 16, 7, 15, 18, 10, 13],
            [15, 12, 14, 16, 0, 9, 15, 3, 9, 10],
            [8, 12, 14, 7, 9, 0, 10, 12, 8, 11],
            [17, 7, 24, 15, 15, 10, 0, 18, 18, 21],
            [17, 13, 15, 18, 3, 12, 18, 0, 10, 11],
            [8, 18, 6, 10, 9, 8, 18, 10, 0, 3],
            [10, 21, 3, 13, 10, 11, 21, 11, 3, 0]]
order = [[9, 2, 37], [4, 10, 36], [2, 5, 26], [5, 3, 48], [7, 10, 39], [9, 7, 38],
         [7, 10, 14], [7, 8, 17], [5, 10, 17], [4, 9, 12], [1, 6, 38], [4, 3, 46]]
vehicle = [[4, 4, 0], [10, 10, 0], [2, 2, 0]] # virtual order
end = [11,11,1e5] # virtual end/order
bigM = 1e5
# make new order matrix
n = len(order)
num_v = len(vehicle)
for x in vehicle:
    order.append(x)
order.append(end)
num_order = n+num_v+1 #length of virtual order⁠
# make new distance matrix
distance.append([0 for i in range(n)])
for x in distance:
    x.append(0)
# c_j: distance in a order j
c=[]
c1=[]
for j in range(num_order):
    c1.append(distance[order[j][0]-1][order[j][1]-1])
    for v in range(num_v):
        c.append(distance[order[j][0]-1][order[j][1]-1])

# =================== Model =========================
try:

    m = Model("Pickup-and-delivery")

    # Create variables: x_ijm; T_im
    x = m.addVars(num_order, num_order,num_v, vtype=GRB.BINARY, name="x")
    T = m.addVars(num_order,num_v, vtype=GRB.CONTINUOUS, name="T")
    Q = m.addVars(num_order, num_order,num_v, vtype=GRB.CONTINUOUS, name="Q")


    # (1) Objecetive funtion
    m.modelSense = GRB.MAXIMIZE
    m.setObjective(LinExpr(c*num_order,x.select())-LinExpr([1e-5]*num_order*num_v,T.select()))

    # Constraints
    # (2) \sum_{m\in M}\sum_{j\in K} x_{ijm} \leqslant 1 & \forall i \in K\cup S
    # Each order is only served once
    for j in range(n+num_v):
        sum_x = 0
        for i in range(n+num_v):
            for v in range(num_v):
                sum_x += x[i,j,v]
        m.addConstr(
            (sum_x <= 1),
            "c1_"+str(i))
    # (2.1)
    for i in range(num_v):
        sum_x = 0
        for j in range(num_v+1):
            for v in range(num_v):
                sum_x += x[n+i,n+j,v]
        m.addConstr(
            (sum_x == 0),
            "c11_"+str(i))
    # (3) \sum_{j. \in K} x_{v_{m}jm} = 1 &\forall m\in M
    # All routes start from initial vehicle node
    for v in range(num_v):
        sum_x = 0
        for j in range(n):
            sum_x += x[n+v,j,v]
        m.addConstr(
            (sum_x == 1),
            "c2_"+str(v))
    # (4) \sum_{j \in U} x_{jim} - \sum_{j \in U} x_{ijm}=0 &\forall i \in K, m\i
    # For each order: "in = out"
    for i in range(n):
        for v in range(num_v):
            m.addConstr(
                (quicksum(x.select('*',i,v)) - quicksum(x.select(i,'*',v)) == 0),
                "c3_"+str(i))
    # (5) \sum_{i\in K} x_{i0m} = 1 & \forall m\in M
    # every vehicle end in 0
    for v in range(num_v):
        sum_x = 0
        for i in range(n):
            sum_x += x[i,n+num_v,v]
        m.addConstr(
            (sum_x == 1),
            "c4_"+str(v))
    # (6) T_{jm} \geqslant (T_{im}+c_i+t_{ij})x_{ijm} \forall i\in U, j\in U,m\in M
    # consistency of time
    m.addConstrs(
        (T[j,m] >= (c1[i]+distance[order[i][1]-1][order[j][0]-1])*x[i,j,m]+Q[i,j,m]
        for i in range(num_order) for j in range(num_order) for m in range(num_v)),
        "c11")
    # (12) Q_{ijm} <= \mathbb{M} x_{ijm} &\forall i\in U, j\in U,m\in M
    m.addConstrs(
        (Q[i,j,m] <= bigM*x[i,j,m]
        for i in range(num_order) for j in range(num_order) for m in range(num_v)),
        'c12')
    # (13) Q_{ijm} <= T_{im} &\forall i\in U,m\in M
    m.addConstrs(
        (Q[i,j,m] <= T[i,m]
        for i in range(num_order) for j in range(num_order) for m in range(num_v)),
        'c13')
    # (14) Q_{ijm} \geqslant T_{im} - (1-x_{ijm})\mathbb{M} &\forall i\in U, j\in U,m\in M
    m.addConstrs(
        (Q[i,j,m] >= T[i,m]-(1-x[i,j,m])*bigM
        for i in range(num_order) for j in range(num_order) for m in range(num_v)),
        'c14')
    # (7) T_{im} \leqslant u_i^p  \forall i \in U,m\in M
    # Time windows
    m.addConstrs(
        (T[i,m] <= order[i][2] for i in range(num_order) for m in range(num_v)),
        "c6")

    # save model and optimize
    # m.write('.\pd.lp')
    m.params.OutputFlag = 0
    m.optimize()

    #========================= Output =========================
    Result = [[] for i in range(num_v)]
    Route = [[] for i in range(num_v)]
    Time = [[] for i in range(num_v)]
    for i in range(num_order):
        for j in range(num_order):
            for v in range(num_v):
                if x[i,j,v].x==1:
                    for u in range(num_v):
                        if u == int(x[i,j,v].varName[-2]):
                            Result[u].append([i,j])
                            Time[u].append(T[i,v].x)
    for u in range(num_v):
        del Time[u][-1]
        Time[u] = sorted((value,counter) for counter, value in enumerate(Time[u]))
        n = len(Time[u])
        Route[u].append(order[Result[u][n][0]][0])
        Route[u].append(0.0)
        for v in Time[u]:
            if v[0] != Route[u][-1]:
                Route[u].append(order[Result[u][v[1]][0]][0])
                Route[u].append(v[0])
            Route[u].append(order[Result[u][v[1]][0]][1])
            Route[u].append(v[0]+c1[Result[u][v[1]][0]])
    obj = LinExpr(c*num_order,x.select()).getValue()
    print('======== Objective Value Found ',str(obj),' =============')
    for v in range(num_v):
        display = ''
        for i in Route[v]:
            if isinstance(i,int):
                display += str(i)
            else:
                display +='('+str(i)+')--'
        print(display)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))
except AttributeError:
    print('Encountered an attribute error')
