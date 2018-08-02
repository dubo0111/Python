\ Model Sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0] + Qk[1]
Subject To
 Q(k)[0]: - gamma[0,0,0] - gamma[0,1,0] - gamma[0,2,0] - delta[0,0,0]
   - delta[0,0,1] - delta[0,1,0] - delta[0,1,1] - delta[0,2,0]
   - delta[0,2,1] - epsilon[0,0] - epsilon[0,1] - epsilon[0,2]
   - lambda[0,1] - lambda[0,2] - mu[0,0] - mu[0,1] + Qk[0] = 0
 Q(k)[1]: - gamma[1,0,0] - gamma[1,1,0] - gamma[1,2,0] - delta[1,0,0]
   - delta[1,0,2] - delta[1,1,0] - delta[1,1,2] - delta[1,2,0]
   - delta[1,2,2] - epsilon[1,0] - epsilon[1,1] - epsilon[1,2]
   - lambda[1,1] - lambda[1,2] - mu[1,0] - mu[1,2] + Qk[1] = 0
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
 v[0,0,1]: 2476.883990998696 beta[0,0] + gamma[0,0,1] + delta[0,0,1]
   + epsilon[0,0] <= 0
 v[0,0,2]: 1460.381032261265 beta[0,0] + gamma[0,0,2] + delta[0,0,2]
   + epsilon[0,0] <= 0
 v[0,1,0]: 735.1241767210855 beta[0,1] + gamma[0,1,0] + delta[0,1,0]
   + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 884.827173373649 beta[0,1] + gamma[0,1,2] + delta[0,1,2]
   + epsilon[0,1] <= 0
 v[0,2,0]: 2433.583308578084 beta[0,2] + gamma[0,2,0] + delta[0,2,0]
   + epsilon[0,2] <= 0
 v[0,2,1]: 1960.501054240609 beta[0,2] + gamma[0,2,1] + delta[0,2,1]
   + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 v[1,0,0]: gamma[1,0,0] + delta[1,0,0] + epsilon[1,0] <= 0
 v[1,0,1]: 3484.877751521346 beta[1,0] + gamma[1,0,1] + delta[1,0,1]
   + epsilon[1,0] <= 0
 v[1,0,2]: 2695.387668233197 beta[1,0] + gamma[1,0,2] + delta[1,0,2]
   + epsilon[1,0] <= 0
 v[1,1,0]: 3504.317238652672 beta[1,1] + gamma[1,1,0] + delta[1,1,0]
   + epsilon[1,1] <= 0
 v[1,1,1]: gamma[1,1,1] + delta[1,1,1] + epsilon[1,1] <= 0
 v[1,1,2]: 1995.951509208602 beta[1,1] + gamma[1,1,2] + delta[1,1,2]
   + epsilon[1,1] <= 0
 v[1,2,0]: 1402.50402833833 beta[1,2] + gamma[1,2,0] + delta[1,2,0]
   + epsilon[1,2] <= 0
 v[1,2,1]: 1001.011982981192 beta[1,2] + gamma[1,2,1] + delta[1,2,1]
   + epsilon[1,2] <= 0
 v[1,2,2]: gamma[1,2,2] + delta[1,2,2] + epsilon[1,2] <= 0
 L3[0]: - beta[0,0] - beta[0,1] - beta[0,2] <= 1
 L3[1]: - beta[1,0] - beta[1,1] - beta[1,2] <= 1
Bounds
 beta[0,0] = 0
 beta[0,1] = 0
 beta[0,2] = 0
 beta[1,0] = 0
 beta[1,1] = 0
 beta[1,2] = 0
 gamma[0,0,0] = 0
 gamma[0,0,1] = 0
 gamma[0,0,2] = 0
 gamma[0,1,0] = 0
 gamma[0,1,1] = 0
 gamma[0,1,2] = 0
 gamma[0,2,0] = 0
 gamma[0,2,1] = 0
 gamma[0,2,2] = 0
 gamma[1,0,0] = 0
 gamma[1,0,1] = 0
 gamma[1,0,2] = 0
 gamma[1,1,0] = 0
 gamma[1,1,1] = 0
 gamma[1,1,2] = 0
 gamma[1,2,0] = 0
 gamma[1,2,1] = 0
 gamma[1,2,2] = 0
 delta[0,0,0] = 0
 delta[0,0,1] = 0
 delta[0,0,2] = 0
 delta[0,1,0] = 0
 delta[0,1,1] = 0
 delta[0,1,2] = 0
 delta[0,2,0] = 0
 delta[0,2,1] = 0
 delta[0,2,2] = 0
 delta[1,0,0] = 0
 delta[1,0,1] = 0
 delta[1,0,2] = 0
 delta[1,1,0] = 0
 delta[1,1,1] = 0
 delta[1,1,2] = 0
 delta[1,2,0] = 0
 delta[1,2,1] = 0
 delta[1,2,2] = 0
 lambda[0,0] = 0
 lambda[0,1] = 0
 lambda[0,2] = 0
 lambda[1,0] = 0
 lambda[1,1] = 0
 lambda[1,2] = 0
 mu[0,0] = 0
 mu[0,1] = 0
 mu[0,2] = 0
 mu[1,0] = 0
 mu[1,1] = 0
 mu[1,2] = 0
End
