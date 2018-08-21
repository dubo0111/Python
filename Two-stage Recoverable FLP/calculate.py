u = [0 for i in range(ni)]
v = [[0 for i in range(ni)] for i in range(ni)]
Q = gamma[0][0] + gamma[1][0] + gamma[2][0] + delta[0][0] + delta[0][1] + delta[1][0] + delta[1][1] + delta[2][0] + delta[2][1] + epsilon[0] + epsilon[1] + epsilon[2] + lamda[1] + lamda[2] + mu[0] + mu[1]
u[0]= -gamma[0][0] - gamma[1][0] - gamma[2][0] + lamda[0] + mu[0] + nu
u[1]= -gamma[0][1] - gamma[1][1] - gamma[2][1] + lamda[1] + mu[1] + nu
u[2]= -gamma[0][2] - gamma[1][2] - gamma[2][2] + lamda[2] + mu[2] + nu
v[0][0] = gamma[0][0] + delta[0][0] + epsilon[0]
v[0][1]= 517 * beta[0] + gamma[0][1] + delta[0][1] + epsilon[0]
v[0][2]= 820 * beta[0] + gamma[0][2] + delta[0][2] + epsilon[0]
v[1][0]= 100 * beta[1] + gamma[1][0] + delta[1][0] + epsilon[1]
v[1][1] = gamma[1][1] + delta[1][1] + epsilon[1]
v[1][2]= 402 * beta[1] + gamma[1][2] + delta[1][2] + epsilon[1]
v[2][0]= 247 * beta[2] + gamma[2][0] + delta[2][0] + epsilon[2]
v[2][1]= 192 * beta[2] + gamma[2][1] + delta[2][1] + epsilon[2]
v[2][2] = gamma[2][2] + delta[2][2] + epsilon[2]
L3 = -beta[0] - beta[1] - beta[2]

Q = gamma[0][1] + gamma[1][1] + gamma[2][1] + delta[0][0]+\
delta[0][1] + delta[1][0] + delta[1][1] + delta[2][0]+\
delta[2][1] + epsilon[0][0] + epsilon[0][1] + epsilon[0][2]+\
lamda[0][0] + lamda[0][2] + mu[0][0] + mu[0][1]