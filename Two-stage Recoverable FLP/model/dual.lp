\ Model sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0]
Subject To
 Q(k)[0]: - gamma[0,0,2] - gamma[0,1,2] - gamma[0,2,2] - delta[0,0,0]
   - delta[0,0,2] - delta[0,1,0] - delta[0,1,2] - delta[0,2,0]
   - delta[0,2,2] - epsilon[0,0] - epsilon[0,1] - epsilon[0,2]
   - lambda[0,0] - lambda[0,1] - mu[0,0] - mu[0,2] + Qk[0] = 0
 u[0,0]: - gamma[0,0,0] - gamma[0,1,0] - gamma[0,2,0] + lambda[0,0]
   + mu[0,0] + nu[0] <= 0
 u[0,1]: - gamma[0,0,1] - gamma[0,1,1] - gamma[0,2,1] + lambda[0,1]
   + mu[0,1] + nu[0] <= 0
 u[0,2]: - gamma[0,0,2] - gamma[0,1,2] - gamma[0,2,2] + lambda[0,2]
   + mu[0,2] + nu[0] <= 0
 v[0,0,0]: gamma[0,0,0] + delta[0,0,0] + epsilon[0,0] <= 0
 v[0,0,1]: 820 beta[0,0] + gamma[0,0,1] + delta[0,0,1] + epsilon[0,0] <= 0
 v[0,0,2]: 402 beta[0,0] + gamma[0,0,2] + delta[0,0,2] + epsilon[0,0] <= 0
 v[0,1,0]: 192 beta[0,1] + gamma[0,1,0] + delta[0,1,0] + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 446 beta[0,1] + gamma[0,1,2] + delta[0,1,2] + epsilon[0,1] <= 0
 v[0,2,0]: 639 beta[0,2] + gamma[0,2,0] + delta[0,2,0] + epsilon[0,2] <= 0
 v[0,2,1]: 785 beta[0,2] + gamma[0,2,1] + delta[0,2,1] + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 L3[0]: - beta[0,0] - beta[0,1] - beta[0,2] <= 1
Bounds
 -infinity <= beta[0,0] <= 0
 -infinity <= beta[0,1] <= 0
 -infinity <= beta[0,2] <= 0
 -infinity <= gamma[0,0,0] <= 0
 -infinity <= gamma[0,0,1] <= 0
 -infinity <= gamma[0,0,2] <= 0
 -infinity <= gamma[0,1,0] <= 0
 -infinity <= gamma[0,1,1] <= 0
 -infinity <= gamma[0,1,2] <= 0
 -infinity <= gamma[0,2,0] <= 0
 -infinity <= gamma[0,2,1] <= 0
 -infinity <= gamma[0,2,2] <= 0
 -infinity <= delta[0,0,0] <= 0
 -infinity <= delta[0,0,1] <= 0
 -infinity <= delta[0,0,2] <= 0
 -infinity <= delta[0,1,0] <= 0
 -infinity <= delta[0,1,1] <= 0
 -infinity <= delta[0,1,2] <= 0
 -infinity <= delta[0,2,0] <= 0
 -infinity <= delta[0,2,1] <= 0
 -infinity <= delta[0,2,2] <= 0
 -infinity <= lambda[0,0] <= 0
 -infinity <= lambda[0,1] <= 0
 -infinity <= lambda[0,2] <= 0
 -infinity <= mu[0,0] <= 0
 -infinity <= mu[0,1] <= 0
 -infinity <= mu[0,2] <= 0
End
