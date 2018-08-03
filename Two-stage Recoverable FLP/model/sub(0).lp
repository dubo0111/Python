\ Model Sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0] + Qk[1]
Subject To
 Q(k)[0]: - gamma[0,0,2] - gamma[0,1,2] - gamma[0,2,2] - delta[0,0,1]
   - delta[0,0,2] - delta[0,1,1] - delta[0,1,2] - delta[0,2,1]
   - delta[0,2,2] - epsilon[0,0] - epsilon[0,1] - epsilon[0,2]
   - lambda[0,0] - lambda[0,1] - mu[0,1] - mu[0,2] + Qk[0] = 0
 Q(k)[1]: - gamma[1,0,2] - gamma[1,1,2] - gamma[1,2,2] - delta[1,0,0]
   - delta[1,0,1] - delta[1,1,0] - delta[1,1,1] - delta[1,2,0]
   - delta[1,2,1] - epsilon[1,0] - epsilon[1,1] - epsilon[1,2]
   - lambda[1,0] - lambda[1,1] - mu[1,0] - mu[1,1] - nu[1] + Qk[1] = 0
 u[0,0,0]: - gamma[0,0,0] + lambda[0,0] + mu[0,0] + nu[0] <= 0
 u[0,0,1]: - gamma[0,0,1] + lambda[0,1] + mu[0,1] + nu[0] <= 0
 u[0,0,2]: - gamma[0,0,2] + lambda[0,2] + mu[0,2] + nu[0] <= 0
 u[0,1,0]: - gamma[0,1,0] + lambda[0,0] + mu[0,0] + nu[0] <= 0
 u[0,1,1]: - gamma[0,1,1] + lambda[0,1] + mu[0,1] + nu[0] <= 0
 u[0,1,2]: - gamma[0,1,2] + lambda[0,2] + mu[0,2] + nu[0] <= 0
 u[0,2,0]: - gamma[0,2,0] + lambda[0,0] + mu[0,0] + nu[0] <= 0
 u[0,2,1]: - gamma[0,2,1] + lambda[0,1] + mu[0,1] + nu[0] <= 0
 u[0,2,2]: - gamma[0,2,2] + lambda[0,2] + mu[0,2] + nu[0] <= 0
 u[1,0,0]: - gamma[1,0,0] + lambda[1,0] + mu[1,0] + nu[1] <= 0
 u[1,0,1]: - gamma[1,0,1] + lambda[1,1] + mu[1,1] + nu[1] <= 0
 u[1,0,2]: - gamma[1,0,2] + lambda[1,2] + mu[1,2] + nu[1] <= 0
 u[1,1,0]: - gamma[1,1,0] + lambda[1,0] + mu[1,0] + nu[1] <= 0
 u[1,1,1]: - gamma[1,1,1] + lambda[1,1] + mu[1,1] + nu[1] <= 0
 u[1,1,2]: - gamma[1,1,2] + lambda[1,2] + mu[1,2] + nu[1] <= 0
 u[1,2,0]: - gamma[1,2,0] + lambda[1,0] + mu[1,0] + nu[1] <= 0
 u[1,2,1]: - gamma[1,2,1] + lambda[1,1] + mu[1,1] + nu[1] <= 0
 u[1,2,2]: - gamma[1,2,2] + lambda[1,2] + mu[1,2] + nu[1] <= 0
 v[0,0,0]: gamma[0,0,0] + delta[0,0,0] + epsilon[0,0] <= 0
 v[0,0,1]: 3246.639321204338 beta[0,0] + gamma[0,0,1] + delta[0,0,1]
   + epsilon[0,0] <= 0
 v[0,0,2]: 981.7931943398505 beta[0,0] + gamma[0,0,2] + delta[0,0,2]
   + epsilon[0,0] <= 0
 v[0,1,0]: 1940.557740501082 beta[0,1] + gamma[0,1,0] + delta[0,1,0]
   + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 650.3715921020806 beta[0,1] + gamma[0,1,2] + delta[0,1,2]
   + epsilon[0,1] <= 0
 v[0,2,0]: 1682.100427706814 beta[0,2] + gamma[0,2,0] + delta[0,2,0]
   + epsilon[0,2] <= 0
 v[0,2,1]: 1275.073818342067 beta[0,2] + gamma[0,2,1] + delta[0,2,1]
   + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 v[1,0,0]: gamma[1,0,0] + delta[1,0,0] + epsilon[1,0] <= 0
 v[1,0,1]: 5182.143562254829 beta[1,0] + gamma[1,0,1] + delta[1,0,1]
   + epsilon[1,0] <= 0
 v[1,0,2]: 1463.339463868784 beta[1,0] + gamma[1,0,2] + delta[1,0,2]
   + epsilon[1,0] <= 0
 v[1,1,0]: 2505.192841739722 beta[1,1] + gamma[1,1,0] + delta[1,1,0]
   + epsilon[1,1] <= 0
 v[1,1,1]: gamma[1,1,1] + delta[1,1,1] + epsilon[1,1] <= 0
 v[1,1,2]: 1028.655438081741 beta[1,1] + gamma[1,1,2] + delta[1,1,2]
   + epsilon[1,1] <= 0
 v[1,2,0]: 1579.067009863555 beta[1,2] + gamma[1,2,0] + delta[1,2,0]
   + epsilon[1,2] <= 0
 v[1,2,1]: 1059.019806043569 beta[1,2] + gamma[1,2,1] + delta[1,2,1]
   + epsilon[1,2] <= 0
 v[1,2,2]: gamma[1,2,2] + delta[1,2,2] + epsilon[1,2] <= 0
 L3[0]: - beta[0,0] - beta[0,1] - beta[0,2] = 1
 L3[1]: - beta[1,0] - beta[1,1] - beta[1,2] = 1
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
End
