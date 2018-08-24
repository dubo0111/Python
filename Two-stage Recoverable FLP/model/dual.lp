\ Model sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0] + Qk[1]
Subject To
 Q(k)[0]: - gamma[0,0,0] - gamma[0,1,0] - gamma[0,2,0] - delta[0,0,0]
   - delta[0,0,1] - delta[0,1,0] - delta[0,1,1] - delta[0,2,0]
   - delta[0,2,1] - epsilon[0,0] - epsilon[0,1] - epsilon[0,2]
   - lambda[0,1] - lambda[0,2] - mu[0,0] - mu[0,1] + Qk[0] = 0
 Q(k)[1]: - gamma[1,0,0] - gamma[1,1,0] - gamma[1,2,0] - delta[1,0,0]
   - delta[1,0,1] - delta[1,1,0] - delta[1,1,1] - delta[1,2,0]
   - delta[1,2,1] - epsilon[1,0] - epsilon[1,1] - epsilon[1,2]
   - lambda[1,1] - lambda[1,2] - mu[1,0] - mu[1,1] + Qk[1] = 0
 u[0,0]: - gamma[0,0,0] - gamma[0,1,0] - gamma[0,2,0] + lambda[0,0]
   + mu[0,0] + nu[0] <= 0
 u[0,1]: - gamma[0,0,1] - gamma[0,1,1] - gamma[0,2,1] + lambda[0,1]
   + mu[0,1] + nu[0] <= 0
 u[0,2]: - gamma[0,0,2] - gamma[0,1,2] - gamma[0,2,2] + lambda[0,2]
   + mu[0,2] + nu[0] <= 0
 u[1,0]: - gamma[1,0,0] - gamma[1,1,0] - gamma[1,2,0] + lambda[1,0]
   + mu[1,0] + nu[1] <= 0
 u[1,1]: - gamma[1,0,1] - gamma[1,1,1] - gamma[1,2,1] + lambda[1,1]
   + mu[1,1] + nu[1] <= 0
 u[1,2]: - gamma[1,0,2] - gamma[1,1,2] - gamma[1,2,2] + lambda[1,2]
   + mu[1,2] + nu[1] <= 0
 v[0,0,0]: gamma[0,0,0] + delta[0,0,0] + epsilon[0,0] <= 0
 v[0,0,1]: 880 beta[0,0] + gamma[0,0,1] + delta[0,0,1] + epsilon[0,0] <= 0
 v[0,0,2]: 1078 beta[0,0] + gamma[0,0,2] + delta[0,0,2] + epsilon[0,0]
   <= 0
 v[0,1,0]: 172 beta[0,1] + gamma[0,1,0] + delta[0,1,0] + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 779 beta[0,1] + gamma[0,1,2] + delta[0,1,2] + epsilon[0,1] <= 0
 v[0,2,0]: 166 beta[0,2] + gamma[0,2,0] + delta[0,2,0] + epsilon[0,2] <= 0
 v[0,2,1]: 313 beta[0,2] + gamma[0,2,1] + delta[0,2,1] + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 v[1,0,0]: gamma[1,0,0] + delta[1,0,0] + epsilon[1,0] <= 0
 v[1,0,1]: 538 beta[1,0] + gamma[1,0,1] + delta[1,0,1] + epsilon[1,0] <= 0
 v[1,0,2]: 638 beta[1,0] + gamma[1,0,2] + delta[1,0,2] + epsilon[1,0] <= 0
 v[1,1,0]: 368 beta[1,1] + gamma[1,1,0] + delta[1,1,0] + epsilon[1,1] <= 0
 v[1,1,1]: gamma[1,1,1] + delta[1,1,1] + epsilon[1,1] <= 0
 v[1,1,2]: 904 beta[1,1] + gamma[1,1,2] + delta[1,1,2] + epsilon[1,1] <= 0
 v[1,2,0]: 388 beta[1,2] + gamma[1,2,0] + delta[1,2,0] + epsilon[1,2] <= 0
 v[1,2,1]: 552 beta[1,2] + gamma[1,2,1] + delta[1,2,1] + epsilon[1,2] <= 0
 v[1,2,2]: gamma[1,2,2] + delta[1,2,2] + epsilon[1,2] <= 0
 L3[0]: - beta[0,0] - beta[0,1] - beta[0,2] <= 1
 L3[1]: - beta[1,0] - beta[1,1] - beta[1,2] <= 1
Bounds
 -infinity <= beta[0,0] <= 0
 -infinity <= beta[0,1] <= 0
 -infinity <= beta[0,2] <= 0
 -infinity <= beta[1,0] <= 0
 -infinity <= beta[1,1] <= 0
 -infinity <= beta[1,2] <= 0
 -infinity <= gamma[0,0,0] <= 0
 -infinity <= gamma[0,0,1] <= 0
 -infinity <= gamma[0,0,2] <= 0
 -infinity <= gamma[0,1,0] <= 0
 -infinity <= gamma[0,1,1] <= 0
 -infinity <= gamma[0,1,2] <= 0
 -infinity <= gamma[0,2,0] <= 0
 -infinity <= gamma[0,2,1] <= 0
 -infinity <= gamma[0,2,2] <= 0
 -infinity <= gamma[1,0,0] <= 0
 -infinity <= gamma[1,0,1] <= 0
 -infinity <= gamma[1,0,2] <= 0
 -infinity <= gamma[1,1,0] <= 0
 -infinity <= gamma[1,1,1] <= 0
 -infinity <= gamma[1,1,2] <= 0
 -infinity <= gamma[1,2,0] <= 0
 -infinity <= gamma[1,2,1] <= 0
 -infinity <= gamma[1,2,2] <= 0
 -infinity <= delta[0,0,0] <= 0
 -infinity <= delta[0,0,1] <= 0
 -infinity <= delta[0,0,2] <= 0
 -infinity <= delta[0,1,0] <= 0
 -infinity <= delta[0,1,1] <= 0
 -infinity <= delta[0,1,2] <= 0
 -infinity <= delta[0,2,0] <= 0
 -infinity <= delta[0,2,1] <= 0
 -infinity <= delta[0,2,2] <= 0
 -infinity <= delta[1,0,0] <= 0
 -infinity <= delta[1,0,1] <= 0
 -infinity <= delta[1,0,2] <= 0
 -infinity <= delta[1,1,0] <= 0
 -infinity <= delta[1,1,1] <= 0
 -infinity <= delta[1,1,2] <= 0
 -infinity <= delta[1,2,0] <= 0
 -infinity <= delta[1,2,1] <= 0
 -infinity <= delta[1,2,2] <= 0
 epsilon[0,0] free
 epsilon[0,1] free
 epsilon[0,2] free
 epsilon[1,0] free
 epsilon[1,1] free
 epsilon[1,2] free
 -infinity <= lambda[0,0] <= 0
 -infinity <= lambda[0,1] <= 0
 -infinity <= lambda[0,2] <= 0
 -infinity <= lambda[1,0] <= 0
 -infinity <= lambda[1,1] <= 0
 -infinity <= lambda[1,2] <= 0
 -infinity <= mu[0,0] <= 0
 -infinity <= mu[0,1] <= 0
 -infinity <= mu[0,2] <= 0
 -infinity <= mu[1,0] <= 0
 -infinity <= mu[1,1] <= 0
 -infinity <= mu[1,2] <= 0
 nu[0] free
 nu[1] free
End
