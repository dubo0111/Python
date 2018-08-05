\ Model sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0]
Subject To
 Q(k)[0]: - gamma[0][0][2] - gamma[0][1][2] - gamma[0][2][2] - delta[0][0][0]\
   - delta[0][0][1] - delta[0][1][0] - delta[0][1][1] - delta[0][2][0]\
   - delta[0][2][1] - epsilon[0][0] - epsilon[0][1] - epsilon[0][2]\
   - lamda[0][0] - lamda[0][1] - mu[0][0] - mu[0][1] - nu[0] + Qk[0] = 0
 u[0][0]: - gamma[0][0][0] - gamma[0][1][0] - gamma[0][2][0] + lamda[0][0]
   + mu[0][0] + nu[0] <= 0
 u[0][1]: - gamma[0][0][1] - gamma[0][1][1] - gamma[0][2][1] + lamda[0][1]
   + mu[0][1] + nu[0] <= 0
 u[0][2]: - gamma[0][0][2] - gamma[0][1][2] - gamma[0][2][2] + lamda[0][2]
   + mu[0][2] + nu[0] <= 0
 v[0][0][0]: gamma[0][0][0] + delta[0][0][0] + epsilon[0][0] <= 0
 v[0][0][1]: 2272.940342393084 beta[0][0] + gamma[0][0][1] + delta[0][0][1]
   + epsilon[0][0] <= 0
 v[0][0][2]: 864.7174334865325 beta[0][0] + gamma[0][0][2] + delta[0][0][2]
   + epsilon[0][0] <= 0
 v[0][1][0]: 3153.719990881915 beta[0][1] + gamma[0][1][0] + delta[0][1][0]
   + epsilon[0][1] <= 0
 v[0][1][1]: gamma[0][1][1] + delta[0][1][1] + epsilon[0][1] <= 0
 v[0][1][2]: 1121.977024482929 beta[0][1] + gamma[0][1][2] + delta[0][1][2]
   + epsilon[0][1] <= 0
 v[0][2][0]: 3101.800431999184 beta[0][2] + gamma[0][2][0] + delta[0][2][0]
   + epsilon[0][2] <= 0
 v[0][2][1]: 1535.200989421115 beta[0][2] + gamma[0][2][1] + delta[0][2][1]
   + epsilon[0][2] <= 0
 v[0][2][2]: gamma[0][2][2] + delta[0][2][2] + epsilon[0][2] <= 0
 L3[0]: - beta[0][0] - beta[0][1] - beta[0][2] <= 1
 R14: beta[0][0] = -0.4080175047827614
 R15: beta[0][2] = -0.5919824952172386
 R16: gamma[0][0][0] = -927.3994470233015
 R17: gamma[0][2][1] = -927.3994470233015
 R18: delta[0][0][2] = -574.579597469973
 R19: delta[0][2][2] = -1836.211559400786
 R20: epsilon[0][0] = 927.3994470233015
 R21: epsilon[0][2] = 1836.211559400786
 R22: nu[0] = -927.3994470233015
Bounds
 -infinity <= beta[0][0] <= 0
 -infinity <= beta[0][1] <= 0
 -infinity <= beta[0][2] <= 0
 -infinity <= gamma[0][0][0] <= 0
 -infinity <= gamma[0][0][1] <= 0
 -infinity <= gamma[0][0][2] <= 0
 -infinity <= gamma[0][1][0] <= 0
 -infinity <= gamma[0][1][1] <= 0
 -infinity <= gamma[0][1][2] <= 0
 -infinity <= gamma[0][2][0] <= 0
 -infinity <= gamma[0][2][1] <= 0
 -infinity <= gamma[0][2][2] <= 0
 -infinity <= delta[0][0][0] <= 0
 -infinity <= delta[0][0][1] <= 0
 -infinity <= delta[0][0][2] <= 0
 -infinity <= delta[0][1][0] <= 0
 -infinity <= delta[0][1][1] <= 0
 -infinity <= delta[0][1][2] <= 0
 -infinity <= delta[0][2][0] <= 0
 -infinity <= delta[0][2][1] <= 0
 -infinity <= delta[0][2][2] <= 0
 -infinity <= lamda[0][0] <= 0
 -infinity <= lamda[0][1] <= 0
 -infinity <= lamda[0][2] <= 0
 -infinity <= mu[0][0] <= 0
 -infinity <= mu[0][1] <= 0
 -infinity <= mu[0][2] <= 0
End
