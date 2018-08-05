\ Model sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0] + Qk[1] + Qk[2] + Qk[3] + Qk[4]
Subject To
 Q(k)[0]: - gamma[0,0,0] - gamma[0,0,7] - gamma[0,0,8] - gamma[0,1,0]
   - gamma[0,1,7] - gamma[0,1,8] - gamma[0,2,0] - gamma[0,2,7]
   - gamma[0,2,8] - gamma[0,3,0] - gamma[0,3,7] - gamma[0,3,8]
   - gamma[0,4,0] - gamma[0,4,7] - gamma[0,4,8] - gamma[0,5,0]
   - gamma[0,5,7] - gamma[0,5,8] - gamma[0,6,0] - gamma[0,6,7]
   - gamma[0,6,8] - gamma[0,7,0] - gamma[0,7,7] - gamma[0,7,8]
   - gamma[0,8,0] - gamma[0,8,7] - gamma[0,8,8] - gamma[0,9,0]
   - gamma[0,9,7] - gamma[0,9,8] - delta[0,0,0] - delta[0,0,1]
   - delta[0,0,2] - delta[0,0,4] - delta[0,0,5] - delta[0,0,6]
   - delta[0,0,7] - delta[0,0,8] - delta[0,1,0] - delta[0,1,1]
   - delta[0,1,2] - delta[0,1,4] - delta[0,1,5] - delta[0,1,6]
   - delta[0,1,7] - delta[0,1,8] - delta[0,2,0] - delta[0,2,1]
   - delta[0,2,2] - delta[0,2,4] - delta[0,2,5] - delta[0,2,6]
   - delta[0,2,7] - delta[0,2,8] - delta[0,3,0] - delta[0,3,1]
   - delta[0,3,2] - delta[0,3,4] - delta[0,3,5] - delta[0,3,6]
   - delta[0,3,7] - delta[0,3,8] - delta[0,4,0] - delta[0,4,1]
   - delta[0,4,2] - delta[0,4,4] - delta[0,4,5] - delta[0,4,6]
   - delta[0,4,7] - delta[0,4,8] - delta[0,5,0] - delta[0,5,1]
   - delta[0,5,2] - delta[0,5,4] - delta[0,5,5] - delta[0,5,6]
   - delta[0,5,7] - delta[0,5,8] - delta[0,6,0] - delta[0,6,1]
   - delta[0,6,2] - delta[0,6,4] - delta[0,6,5] - delta[0,6,6]
   - delta[0,6,7] - delta[0,6,8] - delta[0,7,0] - delta[0,7,1]
   - delta[0,7,2] - delta[0,7,4] - delta[0,7,5] - delta[0,7,6]
   - delta[0,7,7] - delta[0,7,8] - delta[0,8,0] - delta[0,8,1]
   - delta[0,8,2] - delta[0,8,4] - delta[0,8,5] - delta[0,8,6]
   - delta[0,8,7] - delta[0,8,8] - delta[0,9,0] - delta[0,9,1]
   - delta[0,9,2] - delta[0,9,4] - delta[0,9,5] - delta[0,9,6]
   - delta[0,9,7] - delta[0,9,8] - epsilon[0,0] - epsilon[0,1]
   - epsilon[0,2] - epsilon[0,3] - epsilon[0,4] - epsilon[0,5]
   - epsilon[0,6] - epsilon[0,7] - epsilon[0,8] - epsilon[0,9]
   - lambda[0,1] - lambda[0,2] - lambda[0,3] - lambda[0,4] - lambda[0,5]
   - lambda[0,6] - lambda[0,9] - mu[0,0] - mu[0,1] - mu[0,2] - mu[0,4]
   - mu[0,5] - mu[0,6] - mu[0,7] - mu[0,8] + Qk[0] = 0
 Q(k)[1]: - gamma[1,0,0] - gamma[1,0,7] - gamma[1,0,8] - gamma[1,1,0]
   - gamma[1,1,7] - gamma[1,1,8] - gamma[1,2,0] - gamma[1,2,7]
   - gamma[1,2,8] - gamma[1,3,0] - gamma[1,3,7] - gamma[1,3,8]
   - gamma[1,4,0] - gamma[1,4,7] - gamma[1,4,8] - gamma[1,5,0]
   - gamma[1,5,7] - gamma[1,5,8] - gamma[1,6,0] - gamma[1,6,7]
   - gamma[1,6,8] - gamma[1,7,0] - gamma[1,7,7] - gamma[1,7,8]
   - gamma[1,8,0] - gamma[1,8,7] - gamma[1,8,8] - gamma[1,9,0]
   - gamma[1,9,7] - gamma[1,9,8] - delta[1,0,0] - delta[1,0,1]
   - delta[1,0,2] - delta[1,0,3] - delta[1,0,4] - delta[1,0,5]
   - delta[1,0,6] - delta[1,0,9] - delta[1,1,0] - delta[1,1,1]
   - delta[1,1,2] - delta[1,1,3] - delta[1,1,4] - delta[1,1,5]
   - delta[1,1,6] - delta[1,1,9] - delta[1,2,0] - delta[1,2,1]
   - delta[1,2,2] - delta[1,2,3] - delta[1,2,4] - delta[1,2,5]
   - delta[1,2,6] - delta[1,2,9] - delta[1,3,0] - delta[1,3,1]
   - delta[1,3,2] - delta[1,3,3] - delta[1,3,4] - delta[1,3,5]
   - delta[1,3,6] - delta[1,3,9] - delta[1,4,0] - delta[1,4,1]
   - delta[1,4,2] - delta[1,4,3] - delta[1,4,4] - delta[1,4,5]
   - delta[1,4,6] - delta[1,4,9] - delta[1,5,0] - delta[1,5,1]
   - delta[1,5,2] - delta[1,5,3] - delta[1,5,4] - delta[1,5,5]
   - delta[1,5,6] - delta[1,5,9] - delta[1,6,0] - delta[1,6,1]
   - delta[1,6,2] - delta[1,6,3] - delta[1,6,4] - delta[1,6,5]
   - delta[1,6,6] - delta[1,6,9] - delta[1,7,0] - delta[1,7,1]
   - delta[1,7,2] - delta[1,7,3] - delta[1,7,4] - delta[1,7,5]
   - delta[1,7,6] - delta[1,7,9] - delta[1,8,0] - delta[1,8,1]
   - delta[1,8,2] - delta[1,8,3] - delta[1,8,4] - delta[1,8,5]
   - delta[1,8,6] - delta[1,8,9] - delta[1,9,0] - delta[1,9,1]
   - delta[1,9,2] - delta[1,9,3] - delta[1,9,4] - delta[1,9,5]
   - delta[1,9,6] - delta[1,9,9] - epsilon[1,0] - epsilon[1,1]
   - epsilon[1,2] - epsilon[1,3] - epsilon[1,4] - epsilon[1,5]
   - epsilon[1,6] - epsilon[1,7] - epsilon[1,8] - epsilon[1,9]
   - lambda[1,1] - lambda[1,2] - lambda[1,3] - lambda[1,4] - lambda[1,5]
   - lambda[1,6] - lambda[1,9] - mu[1,0] - mu[1,1] - mu[1,2] - mu[1,3]
   - mu[1,4] - mu[1,5] - mu[1,6] - mu[1,9] - 2 nu[1] + Qk[1] = 0
 Q(k)[2]: - gamma[2,0,0] - gamma[2,0,7] - gamma[2,0,8] - gamma[2,1,0]
   - gamma[2,1,7] - gamma[2,1,8] - gamma[2,2,0] - gamma[2,2,7]
   - gamma[2,2,8] - gamma[2,3,0] - gamma[2,3,7] - gamma[2,3,8]
   - gamma[2,4,0] - gamma[2,4,7] - gamma[2,4,8] - gamma[2,5,0]
   - gamma[2,5,7] - gamma[2,5,8] - gamma[2,6,0] - gamma[2,6,7]
   - gamma[2,6,8] - gamma[2,7,0] - gamma[2,7,7] - gamma[2,7,8]
   - gamma[2,8,0] - gamma[2,8,7] - gamma[2,8,8] - gamma[2,9,0]
   - gamma[2,9,7] - gamma[2,9,8] - delta[2,0,0] - delta[2,0,1]
   - delta[2,0,2] - delta[2,0,3] - delta[2,0,4] - delta[2,0,5]
   - delta[2,0,6] - delta[2,0,9] - delta[2,1,0] - delta[2,1,1]
   - delta[2,1,2] - delta[2,1,3] - delta[2,1,4] - delta[2,1,5]
   - delta[2,1,6] - delta[2,1,9] - delta[2,2,0] - delta[2,2,1]
   - delta[2,2,2] - delta[2,2,3] - delta[2,2,4] - delta[2,2,5]
   - delta[2,2,6] - delta[2,2,9] - delta[2,3,0] - delta[2,3,1]
   - delta[2,3,2] - delta[2,3,3] - delta[2,3,4] - delta[2,3,5]
   - delta[2,3,6] - delta[2,3,9] - delta[2,4,0] - delta[2,4,1]
   - delta[2,4,2] - delta[2,4,3] - delta[2,4,4] - delta[2,4,5]
   - delta[2,4,6] - delta[2,4,9] - delta[2,5,0] - delta[2,5,1]
   - delta[2,5,2] - delta[2,5,3] - delta[2,5,4] - delta[2,5,5]
   - delta[2,5,6] - delta[2,5,9] - delta[2,6,0] - delta[2,6,1]
   - delta[2,6,2] - delta[2,6,3] - delta[2,6,4] - delta[2,6,5]
   - delta[2,6,6] - delta[2,6,9] - delta[2,7,0] - delta[2,7,1]
   - delta[2,7,2] - delta[2,7,3] - delta[2,7,4] - delta[2,7,5]
   - delta[2,7,6] - delta[2,7,9] - delta[2,8,0] - delta[2,8,1]
   - delta[2,8,2] - delta[2,8,3] - delta[2,8,4] - delta[2,8,5]
   - delta[2,8,6] - delta[2,8,9] - delta[2,9,0] - delta[2,9,1]
   - delta[2,9,2] - delta[2,9,3] - delta[2,9,4] - delta[2,9,5]
   - delta[2,9,6] - delta[2,9,9] - epsilon[2,0] - epsilon[2,1]
   - epsilon[2,2] - epsilon[2,3] - epsilon[2,4] - epsilon[2,5]
   - epsilon[2,6] - epsilon[2,7] - epsilon[2,8] - epsilon[2,9]
   - lambda[2,1] - lambda[2,2] - lambda[2,3] - lambda[2,4] - lambda[2,5]
   - lambda[2,6] - lambda[2,9] - mu[2,0] - mu[2,1] - mu[2,2] - mu[2,3]
   - mu[2,4] - mu[2,5] - mu[2,6] - mu[2,9] - 2 nu[2] + Qk[2] = 0
 Q(k)[3]: - gamma[3,0,0] - gamma[3,0,7] - gamma[3,0,8] - gamma[3,1,0]
   - gamma[3,1,7] - gamma[3,1,8] - gamma[3,2,0] - gamma[3,2,7]
   - gamma[3,2,8] - gamma[3,3,0] - gamma[3,3,7] - gamma[3,3,8]
   - gamma[3,4,0] - gamma[3,4,7] - gamma[3,4,8] - gamma[3,5,0]
   - gamma[3,5,7] - gamma[3,5,8] - gamma[3,6,0] - gamma[3,6,7]
   - gamma[3,6,8] - gamma[3,7,0] - gamma[3,7,7] - gamma[3,7,8]
   - gamma[3,8,0] - gamma[3,8,7] - gamma[3,8,8] - gamma[3,9,0]
   - gamma[3,9,7] - gamma[3,9,8] - delta[3,0,0] - delta[3,0,1]
   - delta[3,0,2] - delta[3,0,3] - delta[3,0,4] - delta[3,0,6]
   - delta[3,0,8] - delta[3,0,9] - delta[3,1,0] - delta[3,1,1]
   - delta[3,1,2] - delta[3,1,3] - delta[3,1,4] - delta[3,1,6]
   - delta[3,1,8] - delta[3,1,9] - delta[3,2,0] - delta[3,2,1]
   - delta[3,2,2] - delta[3,2,3] - delta[3,2,4] - delta[3,2,6]
   - delta[3,2,8] - delta[3,2,9] - delta[3,3,0] - delta[3,3,1]
   - delta[3,3,2] - delta[3,3,3] - delta[3,3,4] - delta[3,3,6]
   - delta[3,3,8] - delta[3,3,9] - delta[3,4,0] - delta[3,4,1]
   - delta[3,4,2] - delta[3,4,3] - delta[3,4,4] - delta[3,4,6]
   - delta[3,4,8] - delta[3,4,9] - delta[3,5,0] - delta[3,5,1]
   - delta[3,5,2] - delta[3,5,3] - delta[3,5,4] - delta[3,5,6]
   - delta[3,5,8] - delta[3,5,9] - delta[3,6,0] - delta[3,6,1]
   - delta[3,6,2] - delta[3,6,3] - delta[3,6,4] - delta[3,6,6]
   - delta[3,6,8] - delta[3,6,9] - delta[3,7,0] - delta[3,7,1]
   - delta[3,7,2] - delta[3,7,3] - delta[3,7,4] - delta[3,7,6]
   - delta[3,7,8] - delta[3,7,9] - delta[3,8,0] - delta[3,8,1]
   - delta[3,8,2] - delta[3,8,3] - delta[3,8,4] - delta[3,8,6]
   - delta[3,8,8] - delta[3,8,9] - delta[3,9,0] - delta[3,9,1]
   - delta[3,9,2] - delta[3,9,3] - delta[3,9,4] - delta[3,9,6]
   - delta[3,9,8] - delta[3,9,9] - epsilon[3,0] - epsilon[3,1]
   - epsilon[3,2] - epsilon[3,3] - epsilon[3,4] - epsilon[3,5]
   - epsilon[3,6] - epsilon[3,7] - epsilon[3,8] - epsilon[3,9]
   - lambda[3,1] - lambda[3,2] - lambda[3,3] - lambda[3,4] - lambda[3,5]
   - lambda[3,6] - lambda[3,9] - mu[3,0] - mu[3,1] - mu[3,2] - mu[3,3]
   - mu[3,4] - mu[3,6] - mu[3,8] - mu[3,9] - nu[3] + Qk[3] = 0
 Q(k)[4]: - gamma[4,0,0] - gamma[4,0,7] - gamma[4,0,8] - gamma[4,1,0]
   - gamma[4,1,7] - gamma[4,1,8] - gamma[4,2,0] - gamma[4,2,7]
   - gamma[4,2,8] - gamma[4,3,0] - gamma[4,3,7] - gamma[4,3,8]
   - gamma[4,4,0] - gamma[4,4,7] - gamma[4,4,8] - gamma[4,5,0]
   - gamma[4,5,7] - gamma[4,5,8] - gamma[4,6,0] - gamma[4,6,7]
   - gamma[4,6,8] - gamma[4,7,0] - gamma[4,7,7] - gamma[4,7,8]
   - gamma[4,8,0] - gamma[4,8,7] - gamma[4,8,8] - gamma[4,9,0]
   - gamma[4,9,7] - gamma[4,9,8] - delta[4,0,1] - delta[4,0,2]
   - delta[4,0,3] - delta[4,0,4] - delta[4,0,5] - delta[4,0,6]
   - delta[4,0,7] - delta[4,0,8] - delta[4,1,1] - delta[4,1,2]
   - delta[4,1,3] - delta[4,1,4] - delta[4,1,5] - delta[4,1,6]
   - delta[4,1,7] - delta[4,1,8] - delta[4,2,1] - delta[4,2,2]
   - delta[4,2,3] - delta[4,2,4] - delta[4,2,5] - delta[4,2,6]
   - delta[4,2,7] - delta[4,2,8] - delta[4,3,1] - delta[4,3,2]
   - delta[4,3,3] - delta[4,3,4] - delta[4,3,5] - delta[4,3,6]
   - delta[4,3,7] - delta[4,3,8] - delta[4,4,1] - delta[4,4,2]
   - delta[4,4,3] - delta[4,4,4] - delta[4,4,5] - delta[4,4,6]
   - delta[4,4,7] - delta[4,4,8] - delta[4,5,1] - delta[4,5,2]
   - delta[4,5,3] - delta[4,5,4] - delta[4,5,5] - delta[4,5,6]
   - delta[4,5,7] - delta[4,5,8] - delta[4,6,1] - delta[4,6,2]
   - delta[4,6,3] - delta[4,6,4] - delta[4,6,5] - delta[4,6,6]
   - delta[4,6,7] - delta[4,6,8] - delta[4,7,1] - delta[4,7,2]
   - delta[4,7,3] - delta[4,7,4] - delta[4,7,5] - delta[4,7,6]
   - delta[4,7,7] - delta[4,7,8] - delta[4,8,1] - delta[4,8,2]
   - delta[4,8,3] - delta[4,8,4] - delta[4,8,5] - delta[4,8,6]
   - delta[4,8,7] - delta[4,8,8] - delta[4,9,1] - delta[4,9,2]
   - delta[4,9,3] - delta[4,9,4] - delta[4,9,5] - delta[4,9,6]
   - delta[4,9,7] - delta[4,9,8] - epsilon[4,0] - epsilon[4,1]
   - epsilon[4,2] - epsilon[4,3] - epsilon[4,4] - epsilon[4,5]
   - epsilon[4,6] - epsilon[4,7] - epsilon[4,8] - epsilon[4,9]
   - lambda[4,1] - lambda[4,2] - lambda[4,3] - lambda[4,4] - lambda[4,5]
   - lambda[4,6] - lambda[4,9] - mu[4,1] - mu[4,2] - mu[4,3] - mu[4,4]
   - mu[4,5] - mu[4,6] - mu[4,7] - mu[4,8] - nu[4] + Qk[4] = 0
 u[0,0]: - gamma[0,0,0] - gamma[0,1,0] - gamma[0,2,0] - gamma[0,3,0]
   - gamma[0,4,0] - gamma[0,5,0] - gamma[0,6,0] - gamma[0,7,0]
   - gamma[0,8,0] - gamma[0,9,0] + lambda[0,0] + mu[0,0] + nu[0] <= 0
 u[0,1]: - gamma[0,0,1] - gamma[0,1,1] - gamma[0,2,1] - gamma[0,3,1]
   - gamma[0,4,1] - gamma[0,5,1] - gamma[0,6,1] - gamma[0,7,1]
   - gamma[0,8,1] - gamma[0,9,1] + lambda[0,1] + mu[0,1] + nu[0] <= 0
 u[0,2]: - gamma[0,0,2] - gamma[0,1,2] - gamma[0,2,2] - gamma[0,3,2]
   - gamma[0,4,2] - gamma[0,5,2] - gamma[0,6,2] - gamma[0,7,2]
   - gamma[0,8,2] - gamma[0,9,2] + lambda[0,2] + mu[0,2] + nu[0] <= 0
 u[0,3]: - gamma[0,0,3] - gamma[0,1,3] - gamma[0,2,3] - gamma[0,3,3]
   - gamma[0,4,3] - gamma[0,5,3] - gamma[0,6,3] - gamma[0,7,3]
   - gamma[0,8,3] - gamma[0,9,3] + lambda[0,3] + mu[0,3] + nu[0] <= 0
 u[0,4]: - gamma[0,0,4] - gamma[0,1,4] - gamma[0,2,4] - gamma[0,3,4]
   - gamma[0,4,4] - gamma[0,5,4] - gamma[0,6,4] - gamma[0,7,4]
   - gamma[0,8,4] - gamma[0,9,4] + lambda[0,4] + mu[0,4] + nu[0] <= 0
 u[0,5]: - gamma[0,0,5] - gamma[0,1,5] - gamma[0,2,5] - gamma[0,3,5]
   - gamma[0,4,5] - gamma[0,5,5] - gamma[0,6,5] - gamma[0,7,5]
   - gamma[0,8,5] - gamma[0,9,5] + lambda[0,5] + mu[0,5] + nu[0] <= 0
 u[0,6]: - gamma[0,0,6] - gamma[0,1,6] - gamma[0,2,6] - gamma[0,3,6]
   - gamma[0,4,6] - gamma[0,5,6] - gamma[0,6,6] - gamma[0,7,6]
   - gamma[0,8,6] - gamma[0,9,6] + lambda[0,6] + mu[0,6] + nu[0] <= 0
 u[0,7]: - gamma[0,0,7] - gamma[0,1,7] - gamma[0,2,7] - gamma[0,3,7]
   - gamma[0,4,7] - gamma[0,5,7] - gamma[0,6,7] - gamma[0,7,7]
   - gamma[0,8,7] - gamma[0,9,7] + lambda[0,7] + mu[0,7] + nu[0] <= 0
 u[0,8]: - gamma[0,0,8] - gamma[0,1,8] - gamma[0,2,8] - gamma[0,3,8]
   - gamma[0,4,8] - gamma[0,5,8] - gamma[0,6,8] - gamma[0,7,8]
   - gamma[0,8,8] - gamma[0,9,8] + lambda[0,8] + mu[0,8] + nu[0] <= 0
 u[0,9]: - gamma[0,0,9] - gamma[0,1,9] - gamma[0,2,9] - gamma[0,3,9]
   - gamma[0,4,9] - gamma[0,5,9] - gamma[0,6,9] - gamma[0,7,9]
   - gamma[0,8,9] - gamma[0,9,9] + lambda[0,9] + mu[0,9] + nu[0] <= 0
 u[1,0]: - gamma[1,0,0] - gamma[1,1,0] - gamma[1,2,0] - gamma[1,3,0]
   - gamma[1,4,0] - gamma[1,5,0] - gamma[1,6,0] - gamma[1,7,0]
   - gamma[1,8,0] - gamma[1,9,0] + lambda[1,0] + mu[1,0] + nu[1] <= 0
 u[1,1]: - gamma[1,0,1] - gamma[1,1,1] - gamma[1,2,1] - gamma[1,3,1]
   - gamma[1,4,1] - gamma[1,5,1] - gamma[1,6,1] - gamma[1,7,1]
   - gamma[1,8,1] - gamma[1,9,1] + lambda[1,1] + mu[1,1] + nu[1] <= 0
 u[1,2]: - gamma[1,0,2] - gamma[1,1,2] - gamma[1,2,2] - gamma[1,3,2]
   - gamma[1,4,2] - gamma[1,5,2] - gamma[1,6,2] - gamma[1,7,2]
   - gamma[1,8,2] - gamma[1,9,2] + lambda[1,2] + mu[1,2] + nu[1] <= 0
 u[1,3]: - gamma[1,0,3] - gamma[1,1,3] - gamma[1,2,3] - gamma[1,3,3]
   - gamma[1,4,3] - gamma[1,5,3] - gamma[1,6,3] - gamma[1,7,3]
   - gamma[1,8,3] - gamma[1,9,3] + lambda[1,3] + mu[1,3] + nu[1] <= 0
 u[1,4]: - gamma[1,0,4] - gamma[1,1,4] - gamma[1,2,4] - gamma[1,3,4]
   - gamma[1,4,4] - gamma[1,5,4] - gamma[1,6,4] - gamma[1,7,4]
   - gamma[1,8,4] - gamma[1,9,4] + lambda[1,4] + mu[1,4] + nu[1] <= 0
 u[1,5]: - gamma[1,0,5] - gamma[1,1,5] - gamma[1,2,5] - gamma[1,3,5]
   - gamma[1,4,5] - gamma[1,5,5] - gamma[1,6,5] - gamma[1,7,5]
   - gamma[1,8,5] - gamma[1,9,5] + lambda[1,5] + mu[1,5] + nu[1] <= 0
 u[1,6]: - gamma[1,0,6] - gamma[1,1,6] - gamma[1,2,6] - gamma[1,3,6]
   - gamma[1,4,6] - gamma[1,5,6] - gamma[1,6,6] - gamma[1,7,6]
   - gamma[1,8,6] - gamma[1,9,6] + lambda[1,6] + mu[1,6] + nu[1] <= 0
 u[1,7]: - gamma[1,0,7] - gamma[1,1,7] - gamma[1,2,7] - gamma[1,3,7]
   - gamma[1,4,7] - gamma[1,5,7] - gamma[1,6,7] - gamma[1,7,7]
   - gamma[1,8,7] - gamma[1,9,7] + lambda[1,7] + mu[1,7] + nu[1] <= 0
 u[1,8]: - gamma[1,0,8] - gamma[1,1,8] - gamma[1,2,8] - gamma[1,3,8]
   - gamma[1,4,8] - gamma[1,5,8] - gamma[1,6,8] - gamma[1,7,8]
   - gamma[1,8,8] - gamma[1,9,8] + lambda[1,8] + mu[1,8] + nu[1] <= 0
 u[1,9]: - gamma[1,0,9] - gamma[1,1,9] - gamma[1,2,9] - gamma[1,3,9]
   - gamma[1,4,9] - gamma[1,5,9] - gamma[1,6,9] - gamma[1,7,9]
   - gamma[1,8,9] - gamma[1,9,9] + lambda[1,9] + mu[1,9] + nu[1] <= 0
 u[2,0]: - gamma[2,0,0] - gamma[2,1,0] - gamma[2,2,0] - gamma[2,3,0]
   - gamma[2,4,0] - gamma[2,5,0] - gamma[2,6,0] - gamma[2,7,0]
   - gamma[2,8,0] - gamma[2,9,0] + lambda[2,0] + mu[2,0] + nu[2] <= 0
 u[2,1]: - gamma[2,0,1] - gamma[2,1,1] - gamma[2,2,1] - gamma[2,3,1]
   - gamma[2,4,1] - gamma[2,5,1] - gamma[2,6,1] - gamma[2,7,1]
   - gamma[2,8,1] - gamma[2,9,1] + lambda[2,1] + mu[2,1] + nu[2] <= 0
 u[2,2]: - gamma[2,0,2] - gamma[2,1,2] - gamma[2,2,2] - gamma[2,3,2]
   - gamma[2,4,2] - gamma[2,5,2] - gamma[2,6,2] - gamma[2,7,2]
   - gamma[2,8,2] - gamma[2,9,2] + lambda[2,2] + mu[2,2] + nu[2] <= 0
 u[2,3]: - gamma[2,0,3] - gamma[2,1,3] - gamma[2,2,3] - gamma[2,3,3]
   - gamma[2,4,3] - gamma[2,5,3] - gamma[2,6,3] - gamma[2,7,3]
   - gamma[2,8,3] - gamma[2,9,3] + lambda[2,3] + mu[2,3] + nu[2] <= 0
 u[2,4]: - gamma[2,0,4] - gamma[2,1,4] - gamma[2,2,4] - gamma[2,3,4]
   - gamma[2,4,4] - gamma[2,5,4] - gamma[2,6,4] - gamma[2,7,4]
   - gamma[2,8,4] - gamma[2,9,4] + lambda[2,4] + mu[2,4] + nu[2] <= 0
 u[2,5]: - gamma[2,0,5] - gamma[2,1,5] - gamma[2,2,5] - gamma[2,3,5]
   - gamma[2,4,5] - gamma[2,5,5] - gamma[2,6,5] - gamma[2,7,5]
   - gamma[2,8,5] - gamma[2,9,5] + lambda[2,5] + mu[2,5] + nu[2] <= 0
 u[2,6]: - gamma[2,0,6] - gamma[2,1,6] - gamma[2,2,6] - gamma[2,3,6]
   - gamma[2,4,6] - gamma[2,5,6] - gamma[2,6,6] - gamma[2,7,6]
   - gamma[2,8,6] - gamma[2,9,6] + lambda[2,6] + mu[2,6] + nu[2] <= 0
 u[2,7]: - gamma[2,0,7] - gamma[2,1,7] - gamma[2,2,7] - gamma[2,3,7]
   - gamma[2,4,7] - gamma[2,5,7] - gamma[2,6,7] - gamma[2,7,7]
   - gamma[2,8,7] - gamma[2,9,7] + lambda[2,7] + mu[2,7] + nu[2] <= 0
 u[2,8]: - gamma[2,0,8] - gamma[2,1,8] - gamma[2,2,8] - gamma[2,3,8]
   - gamma[2,4,8] - gamma[2,5,8] - gamma[2,6,8] - gamma[2,7,8]
   - gamma[2,8,8] - gamma[2,9,8] + lambda[2,8] + mu[2,8] + nu[2] <= 0
 u[2,9]: - gamma[2,0,9] - gamma[2,1,9] - gamma[2,2,9] - gamma[2,3,9]
   - gamma[2,4,9] - gamma[2,5,9] - gamma[2,6,9] - gamma[2,7,9]
   - gamma[2,8,9] - gamma[2,9,9] + lambda[2,9] + mu[2,9] + nu[2] <= 0
 u[3,0]: - gamma[3,0,0] - gamma[3,1,0] - gamma[3,2,0] - gamma[3,3,0]
   - gamma[3,4,0] - gamma[3,5,0] - gamma[3,6,0] - gamma[3,7,0]
   - gamma[3,8,0] - gamma[3,9,0] + lambda[3,0] + mu[3,0] + nu[3] <= 0
 u[3,1]: - gamma[3,0,1] - gamma[3,1,1] - gamma[3,2,1] - gamma[3,3,1]
   - gamma[3,4,1] - gamma[3,5,1] - gamma[3,6,1] - gamma[3,7,1]
   - gamma[3,8,1] - gamma[3,9,1] + lambda[3,1] + mu[3,1] + nu[3] <= 0
 u[3,2]: - gamma[3,0,2] - gamma[3,1,2] - gamma[3,2,2] - gamma[3,3,2]
   - gamma[3,4,2] - gamma[3,5,2] - gamma[3,6,2] - gamma[3,7,2]
   - gamma[3,8,2] - gamma[3,9,2] + lambda[3,2] + mu[3,2] + nu[3] <= 0
 u[3,3]: - gamma[3,0,3] - gamma[3,1,3] - gamma[3,2,3] - gamma[3,3,3]
   - gamma[3,4,3] - gamma[3,5,3] - gamma[3,6,3] - gamma[3,7,3]
   - gamma[3,8,3] - gamma[3,9,3] + lambda[3,3] + mu[3,3] + nu[3] <= 0
 u[3,4]: - gamma[3,0,4] - gamma[3,1,4] - gamma[3,2,4] - gamma[3,3,4]
   - gamma[3,4,4] - gamma[3,5,4] - gamma[3,6,4] - gamma[3,7,4]
   - gamma[3,8,4] - gamma[3,9,4] + lambda[3,4] + mu[3,4] + nu[3] <= 0
 u[3,5]: - gamma[3,0,5] - gamma[3,1,5] - gamma[3,2,5] - gamma[3,3,5]
   - gamma[3,4,5] - gamma[3,5,5] - gamma[3,6,5] - gamma[3,7,5]
   - gamma[3,8,5] - gamma[3,9,5] + lambda[3,5] + mu[3,5] + nu[3] <= 0
 u[3,6]: - gamma[3,0,6] - gamma[3,1,6] - gamma[3,2,6] - gamma[3,3,6]
   - gamma[3,4,6] - gamma[3,5,6] - gamma[3,6,6] - gamma[3,7,6]
   - gamma[3,8,6] - gamma[3,9,6] + lambda[3,6] + mu[3,6] + nu[3] <= 0
 u[3,7]: - gamma[3,0,7] - gamma[3,1,7] - gamma[3,2,7] - gamma[3,3,7]
   - gamma[3,4,7] - gamma[3,5,7] - gamma[3,6,7] - gamma[3,7,7]
   - gamma[3,8,7] - gamma[3,9,7] + lambda[3,7] + mu[3,7] + nu[3] <= 0
 u[3,8]: - gamma[3,0,8] - gamma[3,1,8] - gamma[3,2,8] - gamma[3,3,8]
   - gamma[3,4,8] - gamma[3,5,8] - gamma[3,6,8] - gamma[3,7,8]
   - gamma[3,8,8] - gamma[3,9,8] + lambda[3,8] + mu[3,8] + nu[3] <= 0
 u[3,9]: - gamma[3,0,9] - gamma[3,1,9] - gamma[3,2,9] - gamma[3,3,9]
   - gamma[3,4,9] - gamma[3,5,9] - gamma[3,6,9] - gamma[3,7,9]
   - gamma[3,8,9] - gamma[3,9,9] + lambda[3,9] + mu[3,9] + nu[3] <= 0
 u[4,0]: - gamma[4,0,0] - gamma[4,1,0] - gamma[4,2,0] - gamma[4,3,0]
   - gamma[4,4,0] - gamma[4,5,0] - gamma[4,6,0] - gamma[4,7,0]
   - gamma[4,8,0] - gamma[4,9,0] + lambda[4,0] + mu[4,0] + nu[4] <= 0
 u[4,1]: - gamma[4,0,1] - gamma[4,1,1] - gamma[4,2,1] - gamma[4,3,1]
   - gamma[4,4,1] - gamma[4,5,1] - gamma[4,6,1] - gamma[4,7,1]
   - gamma[4,8,1] - gamma[4,9,1] + lambda[4,1] + mu[4,1] + nu[4] <= 0
 u[4,2]: - gamma[4,0,2] - gamma[4,1,2] - gamma[4,2,2] - gamma[4,3,2]
   - gamma[4,4,2] - gamma[4,5,2] - gamma[4,6,2] - gamma[4,7,2]
   - gamma[4,8,2] - gamma[4,9,2] + lambda[4,2] + mu[4,2] + nu[4] <= 0
 u[4,3]: - gamma[4,0,3] - gamma[4,1,3] - gamma[4,2,3] - gamma[4,3,3]
   - gamma[4,4,3] - gamma[4,5,3] - gamma[4,6,3] - gamma[4,7,3]
   - gamma[4,8,3] - gamma[4,9,3] + lambda[4,3] + mu[4,3] + nu[4] <= 0
 u[4,4]: - gamma[4,0,4] - gamma[4,1,4] - gamma[4,2,4] - gamma[4,3,4]
   - gamma[4,4,4] - gamma[4,5,4] - gamma[4,6,4] - gamma[4,7,4]
   - gamma[4,8,4] - gamma[4,9,4] + lambda[4,4] + mu[4,4] + nu[4] <= 0
 u[4,5]: - gamma[4,0,5] - gamma[4,1,5] - gamma[4,2,5] - gamma[4,3,5]
   - gamma[4,4,5] - gamma[4,5,5] - gamma[4,6,5] - gamma[4,7,5]
   - gamma[4,8,5] - gamma[4,9,5] + lambda[4,5] + mu[4,5] + nu[4] <= 0
 u[4,6]: - gamma[4,0,6] - gamma[4,1,6] - gamma[4,2,6] - gamma[4,3,6]
   - gamma[4,4,6] - gamma[4,5,6] - gamma[4,6,6] - gamma[4,7,6]
   - gamma[4,8,6] - gamma[4,9,6] + lambda[4,6] + mu[4,6] + nu[4] <= 0
 u[4,7]: - gamma[4,0,7] - gamma[4,1,7] - gamma[4,2,7] - gamma[4,3,7]
   - gamma[4,4,7] - gamma[4,5,7] - gamma[4,6,7] - gamma[4,7,7]
   - gamma[4,8,7] - gamma[4,9,7] + lambda[4,7] + mu[4,7] + nu[4] <= 0
 u[4,8]: - gamma[4,0,8] - gamma[4,1,8] - gamma[4,2,8] - gamma[4,3,8]
   - gamma[4,4,8] - gamma[4,5,8] - gamma[4,6,8] - gamma[4,7,8]
   - gamma[4,8,8] - gamma[4,9,8] + lambda[4,8] + mu[4,8] + nu[4] <= 0
 u[4,9]: - gamma[4,0,9] - gamma[4,1,9] - gamma[4,2,9] - gamma[4,3,9]
   - gamma[4,4,9] - gamma[4,5,9] - gamma[4,6,9] - gamma[4,7,9]
   - gamma[4,8,9] - gamma[4,9,9] + lambda[4,9] + mu[4,9] + nu[4] <= 0
 v[0,0,0]: gamma[0,0,0] + delta[0,0,0] + epsilon[0,0] <= 0
 v[0,0,1]: 1862.815208437965 beta[0,0] + gamma[0,0,1] + delta[0,0,1]
   + epsilon[0,0] <= 0
 v[0,0,2]: 1236.873560460344 beta[0,0] + gamma[0,0,2] + delta[0,0,2]
   + epsilon[0,0] <= 0
 v[0,0,3]: 296.1156143692384 beta[0,0] + gamma[0,0,3] + delta[0,0,3]
   + epsilon[0,0] <= 0
 v[0,0,4]: 99.09995719663213 beta[0,0] + gamma[0,0,4] + delta[0,0,4]
   + epsilon[0,0] <= 0
 v[0,0,5]: 29.18457831202829 beta[0,0] + gamma[0,0,5] + delta[0,0,5]
   + epsilon[0,0] <= 0
 v[0,0,6]: 3315.709605415814 beta[0,0] + gamma[0,0,6] + delta[0,0,6]
   + epsilon[0,0] <= 0
 v[0,0,7]: 1656.509572875603 beta[0,0] + gamma[0,0,7] + delta[0,0,7]
   + epsilon[0,0] <= 0
 v[0,0,8]: 1173.025440341693 beta[0,0] + gamma[0,0,8] + delta[0,0,8]
   + epsilon[0,0] <= 0
 v[0,0,9]: 2123.761305412477 beta[0,0] + gamma[0,0,9] + delta[0,0,9]
   + epsilon[0,0] <= 0
 v[0,1,0]: 1838.516050576543 beta[0,1] + gamma[0,1,0] + delta[0,1,0]
   + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 720.378866837312 beta[0,1] + gamma[0,1,2] + delta[0,1,2]
   + epsilon[0,1] <= 0
 v[0,1,3]: 2115.229233878912 beta[0,1] + gamma[0,1,3] + delta[0,1,3]
   + epsilon[0,1] <= 0
 v[0,1,4]: 2463.880593517438 beta[0,1] + gamma[0,1,4] + delta[0,1,4]
   + epsilon[0,1] <= 0
 v[0,1,5]: 3034.385966082491 beta[0,1] + gamma[0,1,5] + delta[0,1,5]
   + epsilon[0,1] <= 0
 v[0,1,6]: 1783.0751887558 beta[0,1] + gamma[0,1,6] + delta[0,1,6]
   + epsilon[0,1] <= 0
 v[0,1,7]: 1089.880827746051 beta[0,1] + gamma[0,1,7] + delta[0,1,7]
   + epsilon[0,1] <= 0
 v[0,1,8]: 1068.417921363381 beta[0,1] + gamma[0,1,8] + delta[0,1,8]
   + epsilon[0,1] <= 0
 v[0,1,9]: 1482.491638485369 beta[0,1] + gamma[0,1,9] + delta[0,1,9]
   + epsilon[0,1] <= 0
 v[0,2,0]: 2284.568692180269 beta[0,2] + gamma[0,2,0] + delta[0,2,0]
   + epsilon[0,2] <= 0
 v[0,2,1]: 1403.883198013467 beta[0,2] + gamma[0,2,1] + delta[0,2,1]
   + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 v[0,2,3]: 1737.782641677564 beta[0,2] + gamma[0,2,3] + delta[0,2,3]
   + epsilon[0,2] <= 0
 v[0,2,4]: 1268.449254465487 beta[0,2] + gamma[0,2,4] + delta[0,2,4]
   + epsilon[0,2] <= 0
 v[0,2,5]: 1447.311733809632 beta[0,2] + gamma[0,2,5] + delta[0,2,5]
   + epsilon[0,2] <= 0
 v[0,2,6]: 2085.382943602161 beta[0,2] + gamma[0,2,6] + delta[0,2,6]
   + epsilon[0,2] <= 0
 v[0,2,7]: 1963.764622772775 beta[0,2] + gamma[0,2,7] + delta[0,2,7]
   + epsilon[0,2] <= 0
 v[0,2,8]: 424.7428613896919 beta[0,2] + gamma[0,2,8] + delta[0,2,8]
   + epsilon[0,2] <= 0
 v[0,2,9]: 1624.503072902518 beta[0,2] + gamma[0,2,9] + delta[0,2,9]
   + epsilon[0,2] <= 0
 v[0,3,0]: 630.0854215069969 beta[0,3] + gamma[0,3,0] + delta[0,3,0]
   + epsilon[0,3] <= 0
 v[0,3,1]: 2609.756233631956 beta[0,3] + gamma[0,3,1] + delta[0,3,1]
   + epsilon[0,3] <= 0
 v[0,3,2]: 1650.11296356393 beta[0,3] + gamma[0,3,2] + delta[0,3,2]
   + epsilon[0,3] <= 0
 v[0,3,3]: gamma[0,3,3] + delta[0,3,3] + epsilon[0,3] <= 0
 v[0,3,4]: 708.9950163095028 beta[0,3] + gamma[0,3,4] + delta[0,3,4]
   + epsilon[0,3] <= 0
 v[0,3,5]: 734.9023673668759 beta[0,3] + gamma[0,3,5] + delta[0,3,5]
   + epsilon[0,3] <= 0
 v[0,3,6]: 3262.755905786364 beta[0,3] + gamma[0,3,6] + delta[0,3,6]
   + epsilon[0,3] <= 0
 v[0,3,7]: 2219.73083115675 beta[0,3] + gamma[0,3,7] + delta[0,3,7]
   + epsilon[0,3] <= 0
 v[0,3,8]: 2303.762824302708 beta[0,3] + gamma[0,3,8] + delta[0,3,8]
   + epsilon[0,3] <= 0
 v[0,3,9]: 5198.913564326192 beta[0,3] + gamma[0,3,9] + delta[0,3,9]
   + epsilon[0,3] <= 0
 v[0,4,0]: 123.7043530181893 beta[0,4] + gamma[0,4,0] + delta[0,4,0]
   + epsilon[0,4] <= 0
 v[0,4,1]: 3769.73402678317 beta[0,4] + gamma[0,4,1] + delta[0,4,1]
   + epsilon[0,4] <= 0
 v[0,4,2]: 1705.65078120976 beta[0,4] + gamma[0,4,2] + delta[0,4,2]
   + epsilon[0,4] <= 0
 v[0,4,3]: 707.8636878743153 beta[0,4] + gamma[0,4,3] + delta[0,4,3]
   + epsilon[0,4] <= 0
 v[0,4,4]: gamma[0,4,4] + delta[0,4,4] + epsilon[0,4] <= 0
 v[0,4,5]: 82.88580471913262 beta[0,4] + gamma[0,4,5] + delta[0,4,5]
   + epsilon[0,4] <= 0
 v[0,4,6]: 2349.233975154239 beta[0,4] + gamma[0,4,6] + delta[0,4,6]
   + epsilon[0,4] <= 0
 v[0,4,7]: 3786.067999740261 beta[0,4] + gamma[0,4,7] + delta[0,4,7]
   + epsilon[0,4] <= 0
 v[0,4,8]: 2526.143765549908 beta[0,4] + gamma[0,4,8] + delta[0,4,8]
   + epsilon[0,4] <= 0
 v[0,4,9]: 3785.775375455064 beta[0,4] + gamma[0,4,9] + delta[0,4,9]
   + epsilon[0,4] <= 0
 v[0,5,0]: 26.61999085390486 beta[0,5] + gamma[0,5,0] + delta[0,5,0]
   + epsilon[0,5] <= 0
 v[0,5,1]: 2877.229525006581 beta[0,5] + gamma[0,5,1] + delta[0,5,1]
   + epsilon[0,5] <= 0
 v[0,5,2]: 2114.736513256397 beta[0,5] + gamma[0,5,2] + delta[0,5,2]
   + epsilon[0,5] <= 0
 v[0,5,3]: 224.1786295703017 beta[0,5] + gamma[0,5,3] + delta[0,5,3]
   + epsilon[0,5] <= 0
 v[0,5,4]: 107.4010333322114 beta[0,5] + gamma[0,5,4] + delta[0,5,4]
   + epsilon[0,5] <= 0
 v[0,5,5]: gamma[0,5,5] + delta[0,5,5] + epsilon[0,5] <= 0
 v[0,5,6]: 1742.448075361997 beta[0,5] + gamma[0,5,6] + delta[0,5,6]
   + epsilon[0,5] <= 0
 v[0,5,7]: 2401.972327695412 beta[0,5] + gamma[0,5,7] + delta[0,5,7]
   + epsilon[0,5] <= 0
 v[0,5,8]: 1526.42222690684 beta[0,5] + gamma[0,5,8] + delta[0,5,8]
   + epsilon[0,5] <= 0
 v[0,5,9]: 2058.93100150608 beta[0,5] + gamma[0,5,9] + delta[0,5,9]
   + epsilon[0,5] <= 0
 v[0,6,0]: 3938.539287147146 beta[0,6] + gamma[0,6,0] + delta[0,6,0]
   + epsilon[0,6] <= 0
 v[0,6,1]: 449.0077665147642 beta[0,6] + gamma[0,6,1] + delta[0,6,1]
   + epsilon[0,6] <= 0
 v[0,6,2]: 1807.552600820006 beta[0,6] + gamma[0,6,2] + delta[0,6,2]
   + epsilon[0,6] <= 0
 v[0,6,3]: 1489.153551369234 beta[0,6] + gamma[0,6,3] + delta[0,6,3]
   + epsilon[0,6] <= 0
 v[0,6,4]: 2303.62356582884 beta[0,6] + gamma[0,6,4] + delta[0,6,4]
   + epsilon[0,6] <= 0
 v[0,6,5]: 1472.481901401069 beta[0,6] + gamma[0,6,5] + delta[0,6,5]
   + epsilon[0,6] <= 0
 v[0,6,6]: gamma[0,6,6] + delta[0,6,6] + epsilon[0,6] <= 0
 v[0,6,7]: 1211.753748645113 beta[0,6] + gamma[0,6,7] + delta[0,6,7]
   + epsilon[0,6] <= 0
 v[0,6,8]: 1755.430131021978 beta[0,6] + gamma[0,6,8] + delta[0,6,8]
   + epsilon[0,6] <= 0
 v[0,6,9]: 11.451237614142 beta[0,6] + gamma[0,6,9] + delta[0,6,9]
   + epsilon[0,6] <= 0
 v[0,7,0]: 8286.680846010822 beta[0,7] + gamma[0,7,0] + delta[0,7,0]
   + epsilon[0,7] <= 0
 v[0,7,1]: 3486.461144577385 beta[0,7] + gamma[0,7,1] + delta[0,7,1]
   + epsilon[0,7] <= 0
 v[0,7,2]: 2798.9188015747 beta[0,7] + gamma[0,7,2] + delta[0,7,2]
   + epsilon[0,7] <= 0
 v[0,7,3]: 5063.143526740962 beta[0,7] + gamma[0,7,3] + delta[0,7,3]
   + epsilon[0,7] <= 0
 v[0,7,4]: 6297.204328133464 beta[0,7] + gamma[0,7,4] + delta[0,7,4]
   + epsilon[0,7] <= 0
 v[0,7,5]: 7936.004765404767 beta[0,7] + gamma[0,7,5] + delta[0,7,5]
   + epsilon[0,7] <= 0
 v[0,7,6]: 3285.668819290892 beta[0,7] + gamma[0,7,6] + delta[0,7,6]
   + epsilon[0,7] <= 0
 v[0,7,7]: gamma[0,7,7] + delta[0,7,7] + epsilon[0,7] <= 0
 v[0,7,8]: 1489.280062854308 beta[0,7] + gamma[0,7,8] + delta[0,7,8]
   + epsilon[0,7] <= 0
 v[0,7,9]: 2857.246827417122 beta[0,7] + gamma[0,7,9] + delta[0,7,9]
   + epsilon[0,7] <= 0
 v[0,8,0]: 1420.869355541892 beta[0,8] + gamma[0,8,0] + delta[0,8,0]
   + epsilon[0,8] <= 0
 v[0,8,1]: 1047.26470992235 beta[0,8] + gamma[0,8,1] + delta[0,8,1]
   + epsilon[0,8] <= 0
 v[0,8,2]: 158.7303763469253 beta[0,8] + gamma[0,8,2] + delta[0,8,2]
   + epsilon[0,8] <= 0
 v[0,8,3]: 1240.336952735973 beta[0,8] + gamma[0,8,3] + delta[0,8,3]
   + epsilon[0,8] <= 0
 v[0,8,4]: 1921.041486460027 beta[0,8] + gamma[0,8,4] + delta[0,8,4]
   + epsilon[0,8] <= 0
 v[0,8,5]: 1682.407462459523 beta[0,8] + gamma[0,8,5] + delta[0,8,5]
   + epsilon[0,8] <= 0
 v[0,8,6]: 1131.413833575083 beta[0,8] + gamma[0,8,6] + delta[0,8,6]
   + epsilon[0,8] <= 0
 v[0,8,7]: 1261.551440915978 beta[0,8] + gamma[0,8,7] + delta[0,8,7]
   + epsilon[0,8] <= 0
 v[0,8,8]: gamma[0,8,8] + delta[0,8,8] + epsilon[0,8] <= 0
 v[0,8,9]: 1768.86243166163 beta[0,8] + gamma[0,8,9] + delta[0,8,9]
   + epsilon[0,8] <= 0
 v[0,9,0]: 2111.662334109893 beta[0,9] + gamma[0,9,0] + delta[0,9,0]
   + epsilon[0,9] <= 0
 v[0,9,1]: 837.913478761828 beta[0,9] + gamma[0,9,1] + delta[0,9,1]
   + epsilon[0,9] <= 0
 v[0,9,2]: 928.3962429943935 beta[0,9] + gamma[0,9,2] + delta[0,9,2]
   + epsilon[0,9] <= 0
 v[0,9,3]: 2581.871212111922 beta[0,9] + gamma[0,9,3] + delta[0,9,3]
   + epsilon[0,9] <= 0
 v[0,9,4]: 3412.805441991485 beta[0,9] + gamma[0,9,4] + delta[0,9,4]
   + epsilon[0,9] <= 0
 v[0,9,5]: 2478.368646131803 beta[0,9] + gamma[0,9,5] + delta[0,9,5]
   + epsilon[0,9] <= 0
 v[0,9,6]: 6.850028369486729 beta[0,9] + gamma[0,9,6] + delta[0,9,6]
   + epsilon[0,9] <= 0
 v[0,9,7]: 1338.105654785963 beta[0,9] + gamma[0,9,7] + delta[0,9,7]
   + epsilon[0,9] <= 0
 v[0,9,8]: 1679.859708923254 beta[0,9] + gamma[0,9,8] + delta[0,9,8]
   + epsilon[0,9] <= 0
 v[0,9,9]: gamma[0,9,9] + delta[0,9,9] + epsilon[0,9] <= 0
 v[1,0,0]: gamma[1,0,0] + delta[1,0,0] + epsilon[1,0] <= 0
 v[1,0,1]: 1933.595323979484 beta[1,0] + gamma[1,0,1] + delta[1,0,1]
   + epsilon[1,0] <= 0
 v[1,0,2]: 1477.372351512668 beta[1,0] + gamma[1,0,2] + delta[1,0,2]
   + epsilon[1,0] <= 0
 v[1,0,3]: 404.8129233778066 beta[1,0] + gamma[1,0,3] + delta[1,0,3]
   + epsilon[1,0] <= 0
 v[1,0,4]: 60.15116639567031 beta[1,0] + gamma[1,0,4] + delta[1,0,4]
   + epsilon[1,0] <= 0
 v[1,0,5]: 40.57306784985833 beta[1,0] + gamma[1,0,5] + delta[1,0,5]
   + epsilon[1,0] <= 0
 v[1,0,6]: 3055.178650576386 beta[1,0] + gamma[1,0,6] + delta[1,0,6]
   + epsilon[1,0] <= 0
 v[1,0,7]: 1258.186784860209 beta[1,0] + gamma[1,0,7] + delta[1,0,7]
   + epsilon[1,0] <= 0
 v[1,0,8]: 1564.050576870601 beta[1,0] + gamma[1,0,8] + delta[1,0,8]
   + epsilon[1,0] <= 0
 v[1,0,9]: 2865.072523848362 beta[1,0] + gamma[1,0,9] + delta[1,0,9]
   + epsilon[1,0] <= 0
 v[1,1,0]: 1427.806019026917 beta[1,1] + gamma[1,1,0] + delta[1,1,0]
   + epsilon[1,1] <= 0
 v[1,1,1]: gamma[1,1,1] + delta[1,1,1] + epsilon[1,1] <= 0
 v[1,1,2]: 1631.174029546897 beta[1,1] + gamma[1,1,2] + delta[1,1,2]
   + epsilon[1,1] <= 0
 v[1,1,3]: 1672.744939357541 beta[1,1] + gamma[1,1,3] + delta[1,1,3]
   + epsilon[1,1] <= 0
 v[1,1,4]: 2584.670368835925 beta[1,1] + gamma[1,1,4] + delta[1,1,4]
   + epsilon[1,1] <= 0
 v[1,1,5]: 1406.317962904727 beta[1,1] + gamma[1,1,5] + delta[1,1,5]
   + epsilon[1,1] <= 0
 v[1,1,6]: 1138.357230019408 beta[1,1] + gamma[1,1,6] + delta[1,1,6]
   + epsilon[1,1] <= 0
 v[1,1,7]: 1810.276022757952 beta[1,1] + gamma[1,1,7] + delta[1,1,7]
   + epsilon[1,1] <= 0
 v[1,1,8]: 941.4614012853294 beta[1,1] + gamma[1,1,8] + delta[1,1,8]
   + epsilon[1,1] <= 0
 v[1,1,9]: 925.0948490509082 beta[1,1] + gamma[1,1,9] + delta[1,1,9]
   + epsilon[1,1] <= 0
 v[1,2,0]: 1790.379890345054 beta[1,2] + gamma[1,2,0] + delta[1,2,0]
   + epsilon[1,2] <= 0
 v[1,2,1]: 1271.724693519441 beta[1,2] + gamma[1,2,1] + delta[1,2,1]
   + epsilon[1,2] <= 0
 v[1,2,2]: gamma[1,2,2] + delta[1,2,2] + epsilon[1,2] <= 0
 v[1,2,3]: 1978.214021200922 beta[1,2] + gamma[1,2,3] + delta[1,2,3]
   + epsilon[1,2] <= 0
 v[1,2,4]: 2264.287332578305 beta[1,2] + gamma[1,2,4] + delta[1,2,4]
   + epsilon[1,2] <= 0
 v[1,2,5]: 2011.49838448103 beta[1,2] + gamma[1,2,5] + delta[1,2,5]
   + epsilon[1,2] <= 0
 v[1,2,6]: 2191.059989474819 beta[1,2] + gamma[1,2,6] + delta[1,2,6]
   + epsilon[1,2] <= 0
 v[1,2,7]: 1500.663792070507 beta[1,2] + gamma[1,2,7] + delta[1,2,7]
   + epsilon[1,2] <= 0
 v[1,2,8]: 273.4339350442244 beta[1,2] + gamma[1,2,8] + delta[1,2,8]
   + epsilon[1,2] <= 0
 v[1,2,9]: 2251.745744287417 beta[1,2] + gamma[1,2,9] + delta[1,2,9]
   + epsilon[1,2] <= 0
 v[1,3,0]: 912.9647204799142 beta[1,3] + gamma[1,3,0] + delta[1,3,0]
   + epsilon[1,3] <= 0
 v[1,3,1]: 4170.367260867379 beta[1,3] + gamma[1,3,1] + delta[1,3,1]
   + epsilon[1,3] <= 0
 v[1,3,2]: 2175.147017088202 beta[1,3] + gamma[1,3,2] + delta[1,3,2]
   + epsilon[1,3] <= 0
 v[1,3,3]: gamma[1,3,3] + delta[1,3,3] + epsilon[1,3] <= 0
 v[1,3,4]: 605.6763192987504 beta[1,3] + gamma[1,3,4] + delta[1,3,4]
   + epsilon[1,3] <= 0
 v[1,3,5]: 618.2125299583195 beta[1,3] + gamma[1,3,5] + delta[1,3,5]
   + epsilon[1,3] <= 0
 v[1,3,6]: 5601.518066678201 beta[1,3] + gamma[1,3,6] + delta[1,3,6]
   + epsilon[1,3] <= 0
 v[1,3,7]: 3236.672826307843 beta[1,3] + gamma[1,3,7] + delta[1,3,7]
   + epsilon[1,3] <= 0
 v[1,3,8]: 2611.538167745011 beta[1,3] + gamma[1,3,8] + delta[1,3,8]
   + epsilon[1,3] <= 0
 v[1,3,9]: 5069.590213725021 beta[1,3] + gamma[1,3,9] + delta[1,3,9]
   + epsilon[1,3] <= 0
 v[1,4,0]: 64.49261815926258 beta[1,4] + gamma[1,4,0] + delta[1,4,0]
   + epsilon[1,4] <= 0
 v[1,4,1]: 1470.40552757063 beta[1,4] + gamma[1,4,1] + delta[1,4,1]
   + epsilon[1,4] <= 0
 v[1,4,2]: 825.6780734211648 beta[1,4] + gamma[1,4,2] + delta[1,4,2]
   + epsilon[1,4] <= 0
 v[1,4,3]: 287.8420826536379 beta[1,4] + gamma[1,4,3] + delta[1,4,3]
   + epsilon[1,4] <= 0
 v[1,4,4]: gamma[1,4,4] + delta[1,4,4] + epsilon[1,4] <= 0
 v[1,4,5]: 41.24963426690492 beta[1,4] + gamma[1,4,5] + delta[1,4,5]
   + epsilon[1,4] <= 0
 v[1,4,6]: 2122.544364795155 beta[1,4] + gamma[1,4,6] + delta[1,4,6]
   + epsilon[1,4] <= 0
 v[1,4,7]: 2257.265932859311 beta[1,4] + gamma[1,4,7] + delta[1,4,7]
   + epsilon[1,4] <= 0
 v[1,4,8]: 605.2805339988871 beta[1,4] + gamma[1,4,8] + delta[1,4,8]
   + epsilon[1,4] <= 0
 v[1,4,9]: 1622.668518392664 beta[1,4] + gamma[1,4,9] + delta[1,4,9]
   + epsilon[1,4] <= 0
 v[1,5,0]: 31.09609727861513 beta[1,5] + gamma[1,5,0] + delta[1,5,0]
   + epsilon[1,5] <= 0
 v[1,5,1]: 2949.416085460473 beta[1,5] + gamma[1,5,1] + delta[1,5,1]
   + epsilon[1,5] <= 0
 v[1,5,2]: 1211.771878924038 beta[1,5] + gamma[1,5,2] + delta[1,5,2]
   + epsilon[1,5] <= 0
 v[1,5,3]: 565.8779053489203 beta[1,5] + gamma[1,5,3] + delta[1,5,3]
   + epsilon[1,5] <= 0
 v[1,5,4]: 58.67310970266837 beta[1,5] + gamma[1,5,4] + delta[1,5,4]
   + epsilon[1,5] <= 0
 v[1,5,5]: gamma[1,5,5] + delta[1,5,5] + epsilon[1,5] <= 0
 v[1,5,6]: 2140.274180236246 beta[1,5] + gamma[1,5,6] + delta[1,5,6]
   + epsilon[1,5] <= 0
 v[1,5,7]: 4952.450182619516 beta[1,5] + gamma[1,5,7] + delta[1,5,7]
   + epsilon[1,5] <= 0
 v[1,5,8]: 2139.803882521453 beta[1,5] + gamma[1,5,8] + delta[1,5,8]
   + epsilon[1,5] <= 0
 v[1,5,9]: 2314.174613575601 beta[1,5] + gamma[1,5,9] + delta[1,5,9]
   + epsilon[1,5] <= 0
 v[1,6,0]: 2514.073622871666 beta[1,6] + gamma[1,6,0] + delta[1,6,0]
   + epsilon[1,6] <= 0
 v[1,6,1]: 918.0722267508179 beta[1,6] + gamma[1,6,1] + delta[1,6,1]
   + epsilon[1,6] <= 0
 v[1,6,2]: 1205.744182473689 beta[1,6] + gamma[1,6,2] + delta[1,6,2]
   + epsilon[1,6] <= 0
 v[1,6,3]: 1112.299590909645 beta[1,6] + gamma[1,6,3] + delta[1,6,3]
   + epsilon[1,6] <= 0
 v[1,6,4]: 1794.903698519185 beta[1,6] + gamma[1,6,4] + delta[1,6,4]
   + epsilon[1,6] <= 0
 v[1,6,5]: 1654.457801059092 beta[1,6] + gamma[1,6,5] + delta[1,6,5]
   + epsilon[1,6] <= 0
 v[1,6,6]: gamma[1,6,6] + delta[1,6,6] + epsilon[1,6] <= 0
 v[1,6,7]: 575.3528076479838 beta[1,6] + gamma[1,6,7] + delta[1,6,7]
   + epsilon[1,6] <= 0
 v[1,6,8]: 1103.4009007951 beta[1,6] + gamma[1,6,8] + delta[1,6,8]
   + epsilon[1,6] <= 0
 v[1,6,9]: 9.532436419549775 beta[1,6] + gamma[1,6,9] + delta[1,6,9]
   + epsilon[1,6] <= 0
 v[1,7,0]: 6816.378254458494 beta[1,7] + gamma[1,7,0] + delta[1,7,0]
   + epsilon[1,7] <= 0
 v[1,7,1]: 2279.441690380021 beta[1,7] + gamma[1,7,1] + delta[1,7,1]
   + epsilon[1,7] <= 0
 v[1,7,2]: 2253.198025983241 beta[1,7] + gamma[1,7,2] + delta[1,7,2]
   + epsilon[1,7] <= 0
 v[1,7,3]: 7427.970248577854 beta[1,7] + gamma[1,7,3] + delta[1,7,3]
   + epsilon[1,7] <= 0
 v[1,7,4]: 7930.952682744772 beta[1,7] + gamma[1,7,4] + delta[1,7,4]
   + epsilon[1,7] <= 0
 v[1,7,5]: 3500.634656214528 beta[1,7] + gamma[1,7,5] + delta[1,7,5]
   + epsilon[1,7] <= 0
 v[1,7,6]: 3028.697489698611 beta[1,7] + gamma[1,7,6] + delta[1,7,6]
   + epsilon[1,7] <= 0
 v[1,7,7]: gamma[1,7,7] + delta[1,7,7] + epsilon[1,7] <= 0
 v[1,7,8]: 2340.945523371646 beta[1,7] + gamma[1,7,8] + delta[1,7,8]
   + epsilon[1,7] <= 0
 v[1,7,9]: 3686.877729898059 beta[1,7] + gamma[1,7,9] + delta[1,7,9]
   + epsilon[1,7] <= 0
 v[1,8,0]: 1214.623118795376 beta[1,8] + gamma[1,8,0] + delta[1,8,0]
   + epsilon[1,8] <= 0
 v[1,8,1]: 803.9669475147584 beta[1,8] + gamma[1,8,1] + delta[1,8,1]
   + epsilon[1,8] <= 0
 v[1,8,2]: 132.763611223576 beta[1,8] + gamma[1,8,2] + delta[1,8,2]
   + epsilon[1,8] <= 0
 v[1,8,3]: 832.2165895978626 beta[1,8] + gamma[1,8,3] + delta[1,8,3]
   + epsilon[1,8] <= 0
 v[1,8,4]: 849.6086981483306 beta[1,8] + gamma[1,8,4] + delta[1,8,4]
   + epsilon[1,8] <= 0
 v[1,8,5]: 1015.82044713956 beta[1,8] + gamma[1,8,5] + delta[1,8,5]
   + epsilon[1,8] <= 0
 v[1,8,6]: 804.4632214102749 beta[1,8] + gamma[1,8,6] + delta[1,8,6]
   + epsilon[1,8] <= 0
 v[1,8,7]: 424.4576689029857 beta[1,8] + gamma[1,8,7] + delta[1,8,7]
   + epsilon[1,8] <= 0
 v[1,8,8]: gamma[1,8,8] + delta[1,8,8] + epsilon[1,8] <= 0
 v[1,8,9]: 1342.163959535086 beta[1,8] + gamma[1,8,9] + delta[1,8,9]
   + epsilon[1,8] <= 0
 v[1,9,0]: 5794.774214049502 beta[1,9] + gamma[1,9,0] + delta[1,9,0]
   + epsilon[1,9] <= 0
 v[1,9,1]: 1876.615104781714 beta[1,9] + gamma[1,9,1] + delta[1,9,1]
   + epsilon[1,9] <= 0
 v[1,9,2]: 2519.274643249119 beta[1,9] + gamma[1,9,2] + delta[1,9,2]
   + epsilon[1,9] <= 0
 v[1,9,3]: 4572.029864880024 beta[1,9] + gamma[1,9,3] + delta[1,9,3]
   + epsilon[1,9] <= 0
 v[1,9,4]: 3925.335615407136 beta[1,9] + gamma[1,9,4] + delta[1,9,4]
   + epsilon[1,9] <= 0
 v[1,9,5]: 2219.21728113674 beta[1,9] + gamma[1,9,5] + delta[1,9,5]
   + epsilon[1,9] <= 0
 v[1,9,6]: 16.17175626418757 beta[1,9] + gamma[1,9,6] + delta[1,9,6]
   + epsilon[1,9] <= 0
 v[1,9,7]: 2283.85505003956 beta[1,9] + gamma[1,9,7] + delta[1,9,7]
   + epsilon[1,9] <= 0
 v[1,9,8]: 2942.94676430266 beta[1,9] + gamma[1,9,8] + delta[1,9,8]
   + epsilon[1,9] <= 0
 v[1,9,9]: gamma[1,9,9] + delta[1,9,9] + epsilon[1,9] <= 0
 v[2,0,0]: gamma[2,0,0] + delta[2,0,0] + epsilon[2,0] <= 0
 v[2,0,1]: 3660.638067689231 beta[2,0] + gamma[2,0,1] + delta[2,0,1]
   + epsilon[2,0] <= 0
 v[2,0,2]: 1946.964400712413 beta[2,0] + gamma[2,0,2] + delta[2,0,2]
   + epsilon[2,0] <= 0
 v[2,0,3]: 497.212636265385 beta[2,0] + gamma[2,0,3] + delta[2,0,3]
   + epsilon[2,0] <= 0
 v[2,0,4]: 93.27216564902177 beta[2,0] + gamma[2,0,4] + delta[2,0,4]
   + epsilon[2,0] <= 0
 v[2,0,5]: 55.01318426212527 beta[2,0] + gamma[2,0,5] + delta[2,0,5]
   + epsilon[2,0] <= 0
 v[2,0,6]: 6191.637471919165 beta[2,0] + gamma[2,0,6] + delta[2,0,6]
   + epsilon[2,0] <= 0
 v[2,0,7]: 4841.867801814367 beta[2,0] + gamma[2,0,7] + delta[2,0,7]
   + epsilon[2,0] <= 0
 v[2,0,8]: 2080.013090944409 beta[2,0] + gamma[2,0,8] + delta[2,0,8]
   + epsilon[2,0] <= 0
 v[2,0,9]: 4381.160761535482 beta[2,0] + gamma[2,0,9] + delta[2,0,9]
   + epsilon[2,0] <= 0
 v[2,1,0]: 2929.07938954156 beta[2,1] + gamma[2,1,0] + delta[2,1,0]
   + epsilon[2,1] <= 0
 v[2,1,1]: gamma[2,1,1] + delta[2,1,1] + epsilon[2,1] <= 0
 v[2,1,2]: 1986.227406248878 beta[2,1] + gamma[2,1,2] + delta[2,1,2]
   + epsilon[2,1] <= 0
 v[2,1,3]: 3026.04622700831 beta[2,1] + gamma[2,1,3] + delta[2,1,3]
   + epsilon[2,1] <= 0
 v[2,1,4]: 4242.861020666305 beta[2,1] + gamma[2,1,4] + delta[2,1,4]
   + epsilon[2,1] <= 0
 v[2,1,5]: 4908.004056151604 beta[2,1] + gamma[2,1,5] + delta[2,1,5]
   + epsilon[2,1] <= 0
 v[2,1,6]: 2057.076549857149 beta[2,1] + gamma[2,1,6] + delta[2,1,6]
   + epsilon[2,1] <= 0
 v[2,1,7]: 2316.239639669995 beta[2,1] + gamma[2,1,7] + delta[2,1,7]
   + epsilon[2,1] <= 0
 v[2,1,8]: 1739.152346973895 beta[2,1] + gamma[2,1,8] + delta[2,1,8]
   + epsilon[2,1] <= 0
 v[2,1,9]: 1627.634736586665 beta[2,1] + gamma[2,1,9] + delta[2,1,9]
   + epsilon[2,1] <= 0
 v[2,2,0]: 2664.990449990466 beta[2,2] + gamma[2,2,0] + delta[2,2,0]
   + epsilon[2,2] <= 0
 v[2,2,1]: 1041.948604100477 beta[2,2] + gamma[2,2,1] + delta[2,2,1]
   + epsilon[2,2] <= 0
 v[2,2,2]: gamma[2,2,2] + delta[2,2,2] + epsilon[2,2] <= 0
 v[2,2,3]: 2031.229291032819 beta[2,2] + gamma[2,2,3] + delta[2,2,3]
   + epsilon[2,2] <= 0
 v[2,2,4]: 1956.390519865439 beta[2,2] + gamma[2,2,4] + delta[2,2,4]
   + epsilon[2,2] <= 0
 v[2,2,5]: 2190.192335340239 beta[2,2] + gamma[2,2,5] + delta[2,2,5]
   + epsilon[2,2] <= 0
 v[2,2,6]: 1804.933845765483 beta[2,2] + gamma[2,2,6] + delta[2,2,6]
   + epsilon[2,2] <= 0
 v[2,2,7]: 1415.668561864667 beta[2,2] + gamma[2,2,7] + delta[2,2,7]
   + epsilon[2,2] <= 0
 v[2,2,8]: 151.4780322481716 beta[2,2] + gamma[2,2,8] + delta[2,2,8]
   + epsilon[2,2] <= 0
 v[2,2,9]: 2015.800298313411 beta[2,2] + gamma[2,2,9] + delta[2,2,9]
   + epsilon[2,2] <= 0
 v[2,3,0]: 948.8398123059794 beta[2,3] + gamma[2,3,0] + delta[2,3,0]
   + epsilon[2,3] <= 0
 v[2,3,1]: 2579.94326631866 beta[2,3] + gamma[2,3,1] + delta[2,3,1]
   + epsilon[2,3] <= 0
 v[2,3,2]: 1434.960335941569 beta[2,3] + gamma[2,3,2] + delta[2,3,2]
   + epsilon[2,3] <= 0
 v[2,3,3]: gamma[2,3,3] + delta[2,3,3] + epsilon[2,3] <= 0
 v[2,3,4]: 745.0404266983919 beta[2,3] + gamma[2,3,4] + delta[2,3,4]
   + epsilon[2,3] <= 0
 v[2,3,5]: 572.1602250907004 beta[2,3] + gamma[2,3,5] + delta[2,3,5]
   + epsilon[2,3] <= 0
 v[2,3,6]: 6481.700795021287 beta[2,3] + gamma[2,3,6] + delta[2,3,6]
   + epsilon[2,3] <= 0
 v[2,3,7]: 6885.064610459362 beta[2,3] + gamma[2,3,7] + delta[2,3,7]
   + epsilon[2,3] <= 0
 v[2,3,8]: 2209.98493279812 beta[2,3] + gamma[2,3,8] + delta[2,3,8]
   + epsilon[2,3] <= 0
 v[2,3,9]: 7201.550539450875 beta[2,3] + gamma[2,3,9] + delta[2,3,9]
   + epsilon[2,3] <= 0
 v[2,4,0]: 63.64897656888651 beta[2,4] + gamma[2,4,0] + delta[2,4,0]
   + epsilon[2,4] <= 0
 v[2,4,1]: 1895.883598414145 beta[2,4] + gamma[2,4,1] + delta[2,4,1]
   + epsilon[2,4] <= 0
 v[2,4,2]: 1186.866424870085 beta[2,4] + gamma[2,4,2] + delta[2,4,2]
   + epsilon[2,4] <= 0
 v[2,4,3]: 316.941474247536 beta[2,4] + gamma[2,4,3] + delta[2,4,3]
   + epsilon[2,4] <= 0
 v[2,4,4]: gamma[2,4,4] + delta[2,4,4] + epsilon[2,4] <= 0
 v[2,4,5]: 66.18469446078484 beta[2,4] + gamma[2,4,5] + delta[2,4,5]
   + epsilon[2,4] <= 0
 v[2,4,6]: 2342.343370037313 beta[2,4] + gamma[2,4,6] + delta[2,4,6]
   + epsilon[2,4] <= 0
 v[2,4,7]: 2464.981344927093 beta[2,4] + gamma[2,4,7] + delta[2,4,7]
   + epsilon[2,4] <= 0
 v[2,4,8]: 828.4419886220343 beta[2,4] + gamma[2,4,8] + delta[2,4,8]
   + epsilon[2,4] <= 0
 v[2,4,9]: 4396.251273616275 beta[2,4] + gamma[2,4,9] + delta[2,4,9]
   + epsilon[2,4] <= 0
 v[2,5,0]: 19.47420366994589 beta[2,5] + gamma[2,5,0] + delta[2,5,0]
   + epsilon[2,5] <= 0
 v[2,5,1]: 934.2920782445809 beta[2,5] + gamma[2,5,1] + delta[2,5,1]
   + epsilon[2,5] <= 0
 v[2,5,2]: 1414.902688219595 beta[2,5] + gamma[2,5,2] + delta[2,5,2]
   + epsilon[2,5] <= 0
 v[2,5,3]: 218.4907253936597 beta[2,5] + gamma[2,5,3] + delta[2,5,3]
   + epsilon[2,5] <= 0
 v[2,5,4]: 93.27957559751538 beta[2,5] + gamma[2,5,4] + delta[2,5,4]
   + epsilon[2,5] <= 0
 v[2,5,5]: gamma[2,5,5] + delta[2,5,5] + epsilon[2,5] <= 0
 v[2,5,6]: 2194.50137717486 beta[2,5] + gamma[2,5,6] + delta[2,5,6]
   + epsilon[2,5] <= 0
 v[2,5,7]: 1973.819308396998 beta[2,5] + gamma[2,5,7] + delta[2,5,7]
   + epsilon[2,5] <= 0
 v[2,5,8]: 1163.558972349903 beta[2,5] + gamma[2,5,8] + delta[2,5,8]
   + epsilon[2,5] <= 0
 v[2,5,9]: 1473.729524528642 beta[2,5] + gamma[2,5,9] + delta[2,5,9]
   + epsilon[2,5] <= 0
 v[2,6,0]: 4487.775305191663 beta[2,6] + gamma[2,6,0] + delta[2,6,0]
   + epsilon[2,6] <= 0
 v[2,6,1]: 1221.294796514584 beta[2,6] + gamma[2,6,1] + delta[2,6,1]
   + epsilon[2,6] <= 0
 v[2,6,2]: 1947.783415315459 beta[2,6] + gamma[2,6,2] + delta[2,6,2]
   + epsilon[2,6] <= 0
 v[2,6,3]: 3277.928619778996 beta[2,6] + gamma[2,6,3] + delta[2,6,3]
   + epsilon[2,6] <= 0
 v[2,6,4]: 2672.563712570535 beta[2,6] + gamma[2,6,4] + delta[2,6,4]
   + epsilon[2,6] <= 0
 v[2,6,5]: 2082.664248742235 beta[2,6] + gamma[2,6,5] + delta[2,6,5]
   + epsilon[2,6] <= 0
 v[2,6,6]: gamma[2,6,6] + delta[2,6,6] + epsilon[2,6] <= 0
 v[2,6,7]: 799.2469265497268 beta[2,6] + gamma[2,6,7] + delta[2,6,7]
   + epsilon[2,6] <= 0
 v[2,6,8]: 1235.428682915379 beta[2,6] + gamma[2,6,8] + delta[2,6,8]
   + epsilon[2,6] <= 0
 v[2,6,9]: 12.37791862144148 beta[2,6] + gamma[2,6,9] + delta[2,6,9]
   + epsilon[2,6] <= 0
 v[2,7,0]: 5215.994493475137 beta[2,7] + gamma[2,7,0] + delta[2,7,0]
   + epsilon[2,7] <= 0
 v[2,7,1]: 3102.042703911195 beta[2,7] + gamma[2,7,1] + delta[2,7,1]
   + epsilon[2,7] <= 0
 v[2,7,2]: 1346.518749044297 beta[2,7] + gamma[2,7,2] + delta[2,7,2]
   + epsilon[2,7] <= 0
 v[2,7,3]: 3480.538969584979 beta[2,7] + gamma[2,7,3] + delta[2,7,3]
   + epsilon[2,7] <= 0
 v[2,7,4]: 4399.649873225206 beta[2,7] + gamma[2,7,4] + delta[2,7,4]
   + epsilon[2,7] <= 0
 v[2,7,5]: 5144.053704330322 beta[2,7] + gamma[2,7,5] + delta[2,7,5]
   + epsilon[2,7] <= 0
 v[2,7,6]: 1737.416609389662 beta[2,7] + gamma[2,7,6] + delta[2,7,6]
   + epsilon[2,7] <= 0
 v[2,7,7]: gamma[2,7,7] + delta[2,7,7] + epsilon[2,7] <= 0
 v[2,7,8]: 3067.415131655416 beta[2,7] + gamma[2,7,8] + delta[2,7,8]
   + epsilon[2,7] <= 0
 v[2,7,9]: 2917.312873090447 beta[2,7] + gamma[2,7,9] + delta[2,7,9]
   + epsilon[2,7] <= 0
 v[2,8,0]: 1241.229498895382 beta[2,8] + gamma[2,8,0] + delta[2,8,0]
   + epsilon[2,8] <= 0
 v[2,8,1]: 951.219607005663 beta[2,8] + gamma[2,8,1] + delta[2,8,1]
   + epsilon[2,8] <= 0
 v[2,8,2]: 415.6965064850095 beta[2,8] + gamma[2,8,2] + delta[2,8,2]
   + epsilon[2,8] <= 0
 v[2,8,3]: 1270.501657247982 beta[2,8] + gamma[2,8,3] + delta[2,8,3]
   + epsilon[2,8] <= 0
 v[2,8,4]: 2885.897973701719 beta[2,8] + gamma[2,8,4] + delta[2,8,4]
   + epsilon[2,8] <= 0
 v[2,8,5]: 2946.354266367365 beta[2,8] + gamma[2,8,5] + delta[2,8,5]
   + epsilon[2,8] <= 0
 v[2,8,6]: 2427.813003442387 beta[2,8] + gamma[2,8,6] + delta[2,8,6]
   + epsilon[2,8] <= 0
 v[2,8,7]: 1552.102836288542 beta[2,8] + gamma[2,8,7] + delta[2,8,7]
   + epsilon[2,8] <= 0
 v[2,8,8]: gamma[2,8,8] + delta[2,8,8] + epsilon[2,8] <= 0
 v[2,8,9]: 1906.695457170305 beta[2,8] + gamma[2,8,9] + delta[2,8,9]
   + epsilon[2,8] <= 0
 v[2,9,0]: 5398.772473839408 beta[2,9] + gamma[2,9,0] + delta[2,9,0]
   + epsilon[2,9] <= 0
 v[2,9,1]: 1809.777510896833 beta[2,9] + gamma[2,9,1] + delta[2,9,1]
   + epsilon[2,9] <= 0
 v[2,9,2]: 3144.663744478095 beta[2,9] + gamma[2,9,2] + delta[2,9,2]
   + epsilon[2,9] <= 0
 v[2,9,3]: 5143.638340191013 beta[2,9] + gamma[2,9,3] + delta[2,9,3]
   + epsilon[2,9] <= 0
 v[2,9,4]: 3519.466259083943 beta[2,9] + gamma[2,9,4] + delta[2,9,4]
   + epsilon[2,9] <= 0
 v[2,9,5]: 3193.441025803417 beta[2,9] + gamma[2,9,5] + delta[2,9,5]
   + epsilon[2,9] <= 0
 v[2,9,6]: 17.80295006314239 beta[2,9] + gamma[2,9,6] + delta[2,9,6]
   + epsilon[2,9] <= 0
 v[2,9,7]: 1477.05384018801 beta[2,9] + gamma[2,9,7] + delta[2,9,7]
   + epsilon[2,9] <= 0
 v[2,9,8]: 1963.198063768448 beta[2,9] + gamma[2,9,8] + delta[2,9,8]
   + epsilon[2,9] <= 0
 v[2,9,9]: gamma[2,9,9] + delta[2,9,9] + epsilon[2,9] <= 0
 v[3,0,0]: gamma[3,0,0] + delta[3,0,0] + epsilon[3,0] <= 0
 v[3,0,1]: 3208.978648593075 beta[3,0] + gamma[3,0,1] + delta[3,0,1]
   + epsilon[3,0] <= 0
 v[3,0,2]: 2342.901038996515 beta[3,0] + gamma[3,0,2] + delta[3,0,2]
   + epsilon[3,0] <= 0
 v[3,0,3]: 615.957546999766 beta[3,0] + gamma[3,0,3] + delta[3,0,3]
   + epsilon[3,0] <= 0
 v[3,0,4]: 97.24247201388332 beta[3,0] + gamma[3,0,4] + delta[3,0,4]
   + epsilon[3,0] <= 0
 v[3,0,5]: 52.77535054264121 beta[3,0] + gamma[3,0,5] + delta[3,0,5]
   + epsilon[3,0] <= 0
 v[3,0,6]: 4499.578400287878 beta[3,0] + gamma[3,0,6] + delta[3,0,6]
   + epsilon[3,0] <= 0
 v[3,0,7]: 5266.42465442798 beta[3,0] + gamma[3,0,7] + delta[3,0,7]
   + epsilon[3,0] <= 0
 v[3,0,8]: 1240.822296051684 beta[3,0] + gamma[3,0,8] + delta[3,0,8]
   + epsilon[3,0] <= 0
 v[3,0,9]: 2341.17869378275 beta[3,0] + gamma[3,0,9] + delta[3,0,9]
   + epsilon[3,0] <= 0
 v[3,1,0]: 1925.964377499827 beta[3,1] + gamma[3,1,0] + delta[3,1,0]
   + epsilon[3,1] <= 0
 v[3,1,1]: gamma[3,1,1] + delta[3,1,1] + epsilon[3,1] <= 0
 v[3,1,2]: 930.9740704662021 beta[3,1] + gamma[3,1,2] + delta[3,1,2]
   + epsilon[3,1] <= 0
 v[3,1,3]: 778.9730457869714 beta[3,1] + gamma[3,1,3] + delta[3,1,3]
   + epsilon[3,1] <= 0
 v[3,1,4]: 2116.896152843372 beta[3,1] + gamma[3,1,4] + delta[3,1,4]
   + epsilon[3,1] <= 0
 v[3,1,5]: 2086.147499115753 beta[3,1] + gamma[3,1,5] + delta[3,1,5]
   + epsilon[3,1] <= 0
 v[3,1,6]: 1172.987993571106 beta[3,1] + gamma[3,1,6] + delta[3,1,6]
   + epsilon[3,1] <= 0
 v[3,1,7]: 1250.811614378712 beta[3,1] + gamma[3,1,7] + delta[3,1,7]
   + epsilon[3,1] <= 0
 v[3,1,8]: 861.2724371701593 beta[3,1] + gamma[3,1,8] + delta[3,1,8]
   + epsilon[3,1] <= 0
 v[3,1,9]: 902.1951178579441 beta[3,1] + gamma[3,1,9] + delta[3,1,9]
   + epsilon[3,1] <= 0
 v[3,2,0]: 1279.949553365379 beta[3,2] + gamma[3,2,0] + delta[3,2,0]
   + epsilon[3,2] <= 0
 v[3,2,1]: 1071.504337159917 beta[3,2] + gamma[3,2,1] + delta[3,2,1]
   + epsilon[3,2] <= 0
 v[3,2,2]: gamma[3,2,2] + delta[3,2,2] + epsilon[3,2] <= 0
 v[3,2,3]: 1907.420965260077 beta[3,2] + gamma[3,2,3] + delta[3,2,3]
   + epsilon[3,2] <= 0
 v[3,2,4]: 2123.052116151658 beta[3,2] + gamma[3,2,4] + delta[3,2,4]
   + epsilon[3,2] <= 0
 v[3,2,5]: 1389.047923139421 beta[3,2] + gamma[3,2,5] + delta[3,2,5]
   + epsilon[3,2] <= 0
 v[3,2,6]: 1704.71858922734 beta[3,2] + gamma[3,2,6] + delta[3,2,6]
   + epsilon[3,2] <= 0
 v[3,2,7]: 1035.715227648784 beta[3,2] + gamma[3,2,7] + delta[3,2,7]
   + epsilon[3,2] <= 0
 v[3,2,8]: 399.2571937974993 beta[3,2] + gamma[3,2,8] + delta[3,2,8]
   + epsilon[3,2] <= 0
 v[3,2,9]: 3009.230674273607 beta[3,2] + gamma[3,2,9] + delta[3,2,9]
   + epsilon[3,2] <= 0
 v[3,3,0]: 454.6066407840574 beta[3,3] + gamma[3,3,0] + delta[3,3,0]
   + epsilon[3,3] <= 0
 v[3,3,1]: 3855.030737517592 beta[3,3] + gamma[3,3,1] + delta[3,3,1]
   + epsilon[3,3] <= 0
 v[3,3,2]: 2446.775419617858 beta[3,3] + gamma[3,3,2] + delta[3,3,2]
   + epsilon[3,3] <= 0
 v[3,3,3]: gamma[3,3,3] + delta[3,3,3] + epsilon[3,3] <= 0
 v[3,3,4]: 464.7858443432719 beta[3,3] + gamma[3,3,4] + delta[3,3,4]
   + epsilon[3,3] <= 0
 v[3,3,5]: 372.3452744531189 beta[3,3] + gamma[3,3,5] + delta[3,3,5]
   + epsilon[3,3] <= 0
 v[3,3,6]: 3901.006482797561 beta[3,3] + gamma[3,3,6] + delta[3,3,6]
   + epsilon[3,3] <= 0
 v[3,3,7]: 3447.729010744882 beta[3,3] + gamma[3,3,7] + delta[3,3,7]
   + epsilon[3,3] <= 0
 v[3,3,8]: 1695.056694032141 beta[3,3] + gamma[3,3,8] + delta[3,3,8]
   + epsilon[3,3] <= 0
 v[3,3,9]: 4879.968772116327 beta[3,3] + gamma[3,3,9] + delta[3,3,9]
   + epsilon[3,3] <= 0
 v[3,4,0]: 73.25348996842068 beta[3,4] + gamma[3,4,0] + delta[3,4,0]
   + epsilon[3,4] <= 0
 v[3,4,1]: 1828.772807937364 beta[3,4] + gamma[3,4,1] + delta[3,4,1]
   + epsilon[3,4] <= 0
 v[3,4,2]: 1077.882583263636 beta[3,4] + gamma[3,4,2] + delta[3,4,2]
   + epsilon[3,4] <= 0
 v[3,4,3]: 482.1206888426964 beta[3,4] + gamma[3,4,3] + delta[3,4,3]
   + epsilon[3,4] <= 0
 v[3,4,4]: gamma[3,4,4] + delta[3,4,4] + epsilon[3,4] <= 0
 v[3,4,5]: 95.90624034486871 beta[3,4] + gamma[3,4,5] + delta[3,4,5]
   + epsilon[3,4] <= 0
 v[3,4,6]: 1990.128877350368 beta[3,4] + gamma[3,4,6] + delta[3,4,6]
   + epsilon[3,4] <= 0
 v[3,4,7]: 4104.436175606506 beta[3,4] + gamma[3,4,7] + delta[3,4,7]
   + epsilon[3,4] <= 0
 v[3,4,8]: 2204.957498176135 beta[3,4] + gamma[3,4,8] + delta[3,4,8]
   + epsilon[3,4] <= 0
 v[3,4,9]: 2211.916477704504 beta[3,4] + gamma[3,4,9] + delta[3,4,9]
   + epsilon[3,4] <= 0
 v[3,5,0]: 16.83003991641369 beta[3,5] + gamma[3,5,0] + delta[3,5,0]
   + epsilon[3,5] <= 0
 v[3,5,1]: 1070.580247383931 beta[3,5] + gamma[3,5,1] + delta[3,5,1]
   + epsilon[3,5] <= 0
 v[3,5,2]: 931.118449882295 beta[3,5] + gamma[3,5,2] + delta[3,5,2]
   + epsilon[3,5] <= 0
 v[3,5,3]: 147.6994440664921 beta[3,5] + gamma[3,5,3] + delta[3,5,3]
   + epsilon[3,5] <= 0
 v[3,5,4]: 54.69780994099824 beta[3,5] + gamma[3,5,4] + delta[3,5,4]
   + epsilon[3,5] <= 0
 v[3,5,5]: gamma[3,5,5] + delta[3,5,5] + epsilon[3,5] <= 0
 v[3,5,6]: 2475.074947090337 beta[3,5] + gamma[3,5,6] + delta[3,5,6]
   + epsilon[3,5] <= 0
 v[3,5,7]: 1701.607346276125 beta[3,5] + gamma[3,5,7] + delta[3,5,7]
   + epsilon[3,5] <= 0
 v[3,5,8]: 943.7577695088633 beta[3,5] + gamma[3,5,8] + delta[3,5,8]
   + epsilon[3,5] <= 0
 v[3,5,9]: 2046.982806347881 beta[3,5] + gamma[3,5,9] + delta[3,5,9]
   + epsilon[3,5] <= 0
 v[3,6,0]: 2168.512390160236 beta[3,6] + gamma[3,6,0] + delta[3,6,0]
   + epsilon[3,6] <= 0
 v[3,6,1]: 666.3419071067533 beta[3,6] + gamma[3,6,1] + delta[3,6,1]
   + epsilon[3,6] <= 0
 v[3,6,2]: 1527.380286504414 beta[3,6] + gamma[3,6,2] + delta[3,6,2]
   + epsilon[3,6] <= 0
 v[3,6,3]: 1061.894548203772 beta[3,6] + gamma[3,6,3] + delta[3,6,3]
   + epsilon[3,6] <= 0
 v[3,6,4]: 2217.717010487638 beta[3,6] + gamma[3,6,4] + delta[3,6,4]
   + epsilon[3,6] <= 0
 v[3,6,5]: 1567.02334441915 beta[3,6] + gamma[3,6,5] + delta[3,6,5]
   + epsilon[3,6] <= 0
 v[3,6,6]: gamma[3,6,6] + delta[3,6,6] + epsilon[3,6] <= 0
 v[3,6,7]: 570.7707514486364 beta[3,6] + gamma[3,6,7] + delta[3,6,7]
   + epsilon[3,6] <= 0
 v[3,6,8]: 524.4570999441429 beta[3,6] + gamma[3,6,8] + delta[3,6,8]
   + epsilon[3,6] <= 0
 v[3,6,9]: 7.051740607653295 beta[3,6] + gamma[3,6,9] + delta[3,6,9]
   + epsilon[3,6] <= 0
 v[3,7,0]: 3346.618628688321 beta[3,7] + gamma[3,7,0] + delta[3,7,0]
   + epsilon[3,7] <= 0
 v[3,7,1]: 1349.084968569255 beta[3,7] + gamma[3,7,1] + delta[3,7,1]
   + epsilon[3,7] <= 0
 v[3,7,2]: 1351.483736430492 beta[3,7] + gamma[3,7,2] + delta[3,7,2]
   + epsilon[3,7] <= 0
 v[3,7,3]: 1738.054967136183 beta[3,7] + gamma[3,7,3] + delta[3,7,3]
   + epsilon[3,7] <= 0
 v[3,7,4]: 4744.086091380238 beta[3,7] + gamma[3,7,4] + delta[3,7,4]
   + epsilon[3,7] <= 0
 v[3,7,5]: 3010.026784434425 beta[3,7] + gamma[3,7,5] + delta[3,7,5]
   + epsilon[3,7] <= 0
 v[3,7,6]: 1333.611209088748 beta[3,7] + gamma[3,7,6] + delta[3,7,6]
   + epsilon[3,7] <= 0
 v[3,7,7]: gamma[3,7,7] + delta[3,7,7] + epsilon[3,7] <= 0
 v[3,7,8]: 2185.281679182095 beta[3,7] + gamma[3,7,8] + delta[3,7,8]
   + epsilon[3,7] <= 0
 v[3,7,9]: 1922.584786079756 beta[3,7] + gamma[3,7,9] + delta[3,7,9]
   + epsilon[3,7] <= 0
 v[3,8,0]: 1081.142756980481 beta[3,8] + gamma[3,8,0] + delta[3,8,0]
   + epsilon[3,8] <= 0
 v[3,8,1]: 854.2533244409511 beta[3,8] + gamma[3,8,1] + delta[3,8,1]
   + epsilon[3,8] <= 0
 v[3,8,2]: 132.4288506664379 beta[3,8] + gamma[3,8,2] + delta[3,8,2]
   + epsilon[3,8] <= 0
 v[3,8,3]: 661.0659233326738 beta[3,8] + gamma[3,8,3] + delta[3,8,3]
   + epsilon[3,8] <= 0
 v[3,8,4]: 1299.724518986545 beta[3,8] + gamma[3,8,4] + delta[3,8,4]
   + epsilon[3,8] <= 0
 v[3,8,5]: 931.5969260671786 beta[3,8] + gamma[3,8,5] + delta[3,8,5]
   + epsilon[3,8] <= 0
 v[3,8,6]: 943.8352698816655 beta[3,8] + gamma[3,8,6] + delta[3,8,6]
   + epsilon[3,8] <= 0
 v[3,8,7]: 839.5872275154392 beta[3,8] + gamma[3,8,7] + delta[3,8,7]
   + epsilon[3,8] <= 0
 v[3,8,8]: gamma[3,8,8] + delta[3,8,8] + epsilon[3,8] <= 0
 v[3,8,9]: 2036.400400287301 beta[3,8] + gamma[3,8,9] + delta[3,8,9]
   + epsilon[3,8] <= 0
 v[3,9,0]: 3158.592906968413 beta[3,9] + gamma[3,9,0] + delta[3,9,0]
   + epsilon[3,9] <= 0
 v[3,9,1]: 1617.605452398856 beta[3,9] + gamma[3,9,1] + delta[3,9,1]
   + epsilon[3,9] <= 0
 v[3,9,2]: 3357.297089358304 beta[3,9] + gamma[3,9,2] + delta[3,9,2]
   + epsilon[3,9] <= 0
 v[3,9,3]: 2659.169133048421 beta[3,9] + gamma[3,9,3] + delta[3,9,3]
   + epsilon[3,9] <= 0
 v[3,9,4]: 6272.199045233398 beta[3,9] + gamma[3,9,4] + delta[3,9,4]
   + epsilon[3,9] <= 0
 v[3,9,5]: 4612.561016479242 beta[3,9] + gamma[3,9,5] + delta[3,9,5]
   + epsilon[3,9] <= 0
 v[3,9,6]: 15.88971347344664 beta[3,9] + gamma[3,9,6] + delta[3,9,6]
   + epsilon[3,9] <= 0
 v[3,9,7]: 1010.392834866365 beta[3,9] + gamma[3,9,7] + delta[3,9,7]
   + epsilon[3,9] <= 0
 v[3,9,8]: 2015.964339438871 beta[3,9] + gamma[3,9,8] + delta[3,9,8]
   + epsilon[3,9] <= 0
 v[3,9,9]: gamma[3,9,9] + delta[3,9,9] + epsilon[3,9] <= 0
 v[4,0,0]: gamma[4,0,0] + delta[4,0,0] + epsilon[4,0] <= 0
 v[4,0,1]: 4746.594829819967 beta[4,0] + gamma[4,0,1] + delta[4,0,1]
   + epsilon[4,0] <= 0
 v[4,0,2]: 2591.557696249116 beta[4,0] + gamma[4,0,2] + delta[4,0,2]
   + epsilon[4,0] <= 0
 v[4,0,3]: 862.3496082486221 beta[4,0] + gamma[4,0,3] + delta[4,0,3]
   + epsilon[4,0] <= 0
 v[4,0,4]: 92.37247651013308 beta[4,0] + gamma[4,0,4] + delta[4,0,4]
   + epsilon[4,0] <= 0
 v[4,0,5]: 29.49643160630173 beta[4,0] + gamma[4,0,5] + delta[4,0,5]
   + epsilon[4,0] <= 0
 v[4,0,6]: 2710.637338234039 beta[4,0] + gamma[4,0,6] + delta[4,0,6]
   + epsilon[4,0] <= 0
 v[4,0,7]: 2909.731602801275 beta[4,0] + gamma[4,0,7] + delta[4,0,7]
   + epsilon[4,0] <= 0
 v[4,0,8]: 2047.276354548051 beta[4,0] + gamma[4,0,8] + delta[4,0,8]
   + epsilon[4,0] <= 0
 v[4,0,9]: 6961.198183740335 beta[4,0] + gamma[4,0,9] + delta[4,0,9]
   + epsilon[4,0] <= 0
 v[4,1,0]: 2316.532092494641 beta[4,1] + gamma[4,1,0] + delta[4,1,0]
   + epsilon[4,1] <= 0
 v[4,1,1]: gamma[4,1,1] + delta[4,1,1] + epsilon[4,1] <= 0
 v[4,1,2]: 1114.38353986556 beta[4,1] + gamma[4,1,2] + delta[4,1,2]
   + epsilon[4,1] <= 0
 v[4,1,3]: 1109.750535610367 beta[4,1] + gamma[4,1,3] + delta[4,1,3]
   + epsilon[4,1] <= 0
 v[4,1,4]: 1050.476468414755 beta[4,1] + gamma[4,1,4] + delta[4,1,4]
   + epsilon[4,1] <= 0
 v[4,1,5]: 1633.425847091895 beta[4,1] + gamma[4,1,5] + delta[4,1,5]
   + epsilon[4,1] <= 0
 v[4,1,6]: 1271.899080832162 beta[4,1] + gamma[4,1,6] + delta[4,1,6]
   + epsilon[4,1] <= 0
 v[4,1,7]: 967.9824716157534 beta[4,1] + gamma[4,1,7] + delta[4,1,7]
   + epsilon[4,1] <= 0
 v[4,1,8]: 843.0833266233458 beta[4,1] + gamma[4,1,8] + delta[4,1,8]
   + epsilon[4,1] <= 0
 v[4,1,9]: 709.1321787209818 beta[4,1] + gamma[4,1,9] + delta[4,1,9]
   + epsilon[4,1] <= 0
 v[4,2,0]: 2326.145901274243 beta[4,2] + gamma[4,2,0] + delta[4,2,0]
   + epsilon[4,2] <= 0
 v[4,2,1]: 1082.964646932824 beta[4,2] + gamma[4,2,1] + delta[4,2,1]
   + epsilon[4,2] <= 0
 v[4,2,2]: gamma[4,2,2] + delta[4,2,2] + epsilon[4,2] <= 0
 v[4,2,3]: 1908.67620403395 beta[4,2] + gamma[4,2,3] + delta[4,2,3]
   + epsilon[4,2] <= 0
 v[4,2,4]: 1283.316867547302 beta[4,2] + gamma[4,2,4] + delta[4,2,4]
   + epsilon[4,2] <= 0
 v[4,2,5]: 2918.590720288229 beta[4,2] + gamma[4,2,5] + delta[4,2,5]
   + epsilon[4,2] <= 0
 v[4,2,6]: 2776.922886248433 beta[4,2] + gamma[4,2,6] + delta[4,2,6]
   + epsilon[4,2] <= 0
 v[4,2,7]: 737.955175227482 beta[4,2] + gamma[4,2,7] + delta[4,2,7]
   + epsilon[4,2] <= 0
 v[4,2,8]: 159.3754127435812 beta[4,2] + gamma[4,2,8] + delta[4,2,8]
   + epsilon[4,2] <= 0
 v[4,2,9]: 1973.863215874105 beta[4,2] + gamma[4,2,9] + delta[4,2,9]
   + epsilon[4,2] <= 0
 v[4,3,0]: 525.5322857471776 beta[4,3] + gamma[4,3,0] + delta[4,3,0]
   + epsilon[4,3] <= 0
 v[4,3,1]: 2468.484279331002 beta[4,3] + gamma[4,3,1] + delta[4,3,1]
   + epsilon[4,3] <= 0
 v[4,3,2]: 1503.919144596161 beta[4,3] + gamma[4,3,2] + delta[4,3,2]
   + epsilon[4,3] <= 0
 v[4,3,3]: gamma[4,3,3] + delta[4,3,3] + epsilon[4,3] <= 0
 v[4,3,4]: 587.3961209665413 beta[4,3] + gamma[4,3,4] + delta[4,3,4]
   + epsilon[4,3] <= 0
 v[4,3,5]: 588.9907243345557 beta[4,3] + gamma[4,3,5] + delta[4,3,5]
   + epsilon[4,3] <= 0
 v[4,3,6]: 2135.034143288857 beta[4,3] + gamma[4,3,6] + delta[4,3,6]
   + epsilon[4,3] <= 0
 v[4,3,7]: 3752.971943981723 beta[4,3] + gamma[4,3,7] + delta[4,3,7]
   + epsilon[4,3] <= 0
 v[4,3,8]: 913.2131670227029 beta[4,3] + gamma[4,3,8] + delta[4,3,8]
   + epsilon[4,3] <= 0
 v[4,3,9]: 4912.309183535601 beta[4,3] + gamma[4,3,9] + delta[4,3,9]
   + epsilon[4,3] <= 0
 v[4,4,0]: 90.66860423519273 beta[4,4] + gamma[4,4,0] + delta[4,4,0]
   + epsilon[4,4] <= 0
 v[4,4,1]: 3008.02641127578 beta[4,4] + gamma[4,4,1] + delta[4,4,1]
   + epsilon[4,4] <= 0
 v[4,4,2]: 2128.253324225433 beta[4,4] + gamma[4,4,2] + delta[4,4,2]
   + epsilon[4,4] <= 0
 v[4,4,3]: 600.1657266940094 beta[4,4] + gamma[4,4,3] + delta[4,4,3]
   + epsilon[4,4] <= 0
 v[4,4,4]: gamma[4,4,4] + delta[4,4,4] + epsilon[4,4] <= 0
 v[4,4,5]: 63.72798920026596 beta[4,4] + gamma[4,4,5] + delta[4,4,5]
   + epsilon[4,4] <= 0
 v[4,4,6]: 3394.821975103665 beta[4,4] + gamma[4,4,6] + delta[4,4,6]
   + epsilon[4,4] <= 0
 v[4,4,7]: 2389.4970767903 beta[4,4] + gamma[4,4,7] + delta[4,4,7]
   + epsilon[4,4] <= 0
 v[4,4,8]: 1456.033418173623 beta[4,4] + gamma[4,4,8] + delta[4,4,8]
   + epsilon[4,4] <= 0
 v[4,4,9]: 4835.690005492854 beta[4,4] + gamma[4,4,9] + delta[4,4,9]
   + epsilon[4,4] <= 0
 v[4,5,0]: 57.73677991006338 beta[4,5] + gamma[4,5,0] + delta[4,5,0]
   + epsilon[4,5] <= 0
 v[4,5,1]: 1850.159065954677 beta[4,5] + gamma[4,5,1] + delta[4,5,1]
   + epsilon[4,5] <= 0
 v[4,5,2]: 2342.855708834278 beta[4,5] + gamma[4,5,2] + delta[4,5,2]
   + epsilon[4,5] <= 0
 v[4,5,3]: 667.4216230786274 beta[4,5] + gamma[4,5,3] + delta[4,5,3]
   + epsilon[4,5] <= 0
 v[4,5,4]: 104.6853122635332 beta[4,5] + gamma[4,5,4] + delta[4,5,4]
   + epsilon[4,5] <= 0
 v[4,5,5]: gamma[4,5,5] + delta[4,5,5] + epsilon[4,5] <= 0
 v[4,5,6]: 3886.643547518208 beta[4,5] + gamma[4,5,6] + delta[4,5,6]
   + epsilon[4,5] <= 0
 v[4,5,7]: 2353.751834145549 beta[4,5] + gamma[4,5,7] + delta[4,5,7]
   + epsilon[4,5] <= 0
 v[4,5,8]: 1870.570970547994 beta[4,5] + gamma[4,5,8] + delta[4,5,8]
   + epsilon[4,5] <= 0
 v[4,5,9]: 4261.255701542335 beta[4,5] + gamma[4,5,9] + delta[4,5,9]
   + epsilon[4,5] <= 0
 v[4,6,0]: 1544.149927853291 beta[4,6] + gamma[4,6,0] + delta[4,6,0]
   + epsilon[4,6] <= 0
 v[4,6,1]: 710.555077279873 beta[4,6] + gamma[4,6,1] + delta[4,6,1]
   + epsilon[4,6] <= 0
 v[4,6,2]: 866.2503545094985 beta[4,6] + gamma[4,6,2] + delta[4,6,2]
   + epsilon[4,6] <= 0
 v[4,6,3]: 2006.074740588352 beta[4,6] + gamma[4,6,3] + delta[4,6,3]
   + epsilon[4,6] <= 0
 v[4,6,4]: 1697.911504706406 beta[4,6] + gamma[4,6,4] + delta[4,6,4]
   + epsilon[4,6] <= 0
 v[4,6,5]: 2332.68710821873 beta[4,6] + gamma[4,6,5] + delta[4,6,5]
   + epsilon[4,6] <= 0
 v[4,6,6]: gamma[4,6,6] + delta[4,6,6] + epsilon[4,6] <= 0
 v[4,6,7]: 1166.736987137108 beta[4,6] + gamma[4,6,7] + delta[4,6,7]
   + epsilon[4,6] <= 0
 v[4,6,8]: 1232.385217093077 beta[4,6] + gamma[4,6,8] + delta[4,6,8]
   + epsilon[4,6] <= 0
 v[4,6,9]: 7.305244040538343 beta[4,6] + gamma[4,6,9] + delta[4,6,9]
   + epsilon[4,6] <= 0
 v[4,7,0]: 3363.78095714536 beta[4,7] + gamma[4,7,0] + delta[4,7,0]
   + epsilon[4,7] <= 0
 v[4,7,1]: 1510.493252937328 beta[4,7] + gamma[4,7,1] + delta[4,7,1]
   + epsilon[4,7] <= 0
 v[4,7,2]: 1554.039190034829 beta[4,7] + gamma[4,7,2] + delta[4,7,2]
   + epsilon[4,7] <= 0
 v[4,7,3]: 4230.381314980858 beta[4,7] + gamma[4,7,3] + delta[4,7,3]
   + epsilon[4,7] <= 0
 v[4,7,4]: 3184.290336477887 beta[4,7] + gamma[4,7,4] + delta[4,7,4]
   + epsilon[4,7] <= 0
 v[4,7,5]: 3093.14927923781 beta[4,7] + gamma[4,7,5] + delta[4,7,5]
   + epsilon[4,7] <= 0
 v[4,7,6]: 2259.453010744211 beta[4,7] + gamma[4,7,6] + delta[4,7,6]
   + epsilon[4,7] <= 0
 v[4,7,7]: gamma[4,7,7] + delta[4,7,7] + epsilon[4,7] <= 0
 v[4,7,8]: 2469.341901168487 beta[4,7] + gamma[4,7,8] + delta[4,7,8]
   + epsilon[4,7] <= 0
 v[4,7,9]: 1978.447409743468 beta[4,7] + gamma[4,7,9] + delta[4,7,9]
   + epsilon[4,7] <= 0
 v[4,8,0]: 906.1119161202895 beta[4,8] + gamma[4,8,0] + delta[4,8,0]
   + epsilon[4,8] <= 0
 v[4,8,1]: 585.8922008680116 beta[4,8] + gamma[4,8,1] + delta[4,8,1]
   + epsilon[4,8] <= 0
 v[4,8,2]: 274.1700126295809 beta[4,8] + gamma[4,8,2] + delta[4,8,2]
   + epsilon[4,8] <= 0
 v[4,8,3]: 875.0968403870539 beta[4,8] + gamma[4,8,3] + delta[4,8,3]
   + epsilon[4,8] <= 0
 v[4,8,4]: 2009.078613030907 beta[4,8] + gamma[4,8,4] + delta[4,8,4]
   + epsilon[4,8] <= 0
 v[4,8,5]: 1702.671362661393 beta[4,8] + gamma[4,8,5] + delta[4,8,5]
   + epsilon[4,8] <= 0
 v[4,8,6]: 2052.840771738036 beta[4,8] + gamma[4,8,6] + delta[4,8,6]
   + epsilon[4,8] <= 0
 v[4,8,7]: 764.6919414397835 beta[4,8] + gamma[4,8,7] + delta[4,8,7]
   + epsilon[4,8] <= 0
 v[4,8,8]: gamma[4,8,8] + delta[4,8,8] + epsilon[4,8] <= 0
 v[4,8,9]: 2124.916297762169 beta[4,8] + gamma[4,8,9] + delta[4,8,9]
   + epsilon[4,8] <= 0
 v[4,9,0]: 2340.302044550013 beta[4,9] + gamma[4,9,0] + delta[4,9,0]
   + epsilon[4,9] <= 0
 v[4,9,1]: 1113.049402702473 beta[4,9] + gamma[4,9,1] + delta[4,9,1]
   + epsilon[4,9] <= 0
 v[4,9,2]: 1319.039512064695 beta[4,9] + gamma[4,9,2] + delta[4,9,2]
   + epsilon[4,9] <= 0
 v[4,9,3]: 2513.531787076773 beta[4,9] + gamma[4,9,3] + delta[4,9,3]
   + epsilon[4,9] <= 0
 v[4,9,4]: 2767.247686510893 beta[4,9] + gamma[4,9,4] + delta[4,9,4]
   + epsilon[4,9] <= 0
 v[4,9,5]: 2526.462460782454 beta[4,9] + gamma[4,9,5] + delta[4,9,5]
   + epsilon[4,9] <= 0
 v[4,9,6]: 7.552541475499328 beta[4,9] + gamma[4,9,6] + delta[4,9,6]
   + epsilon[4,9] <= 0
 v[4,9,7]: 820.5329507835795 beta[4,9] + gamma[4,9,7] + delta[4,9,7]
   + epsilon[4,9] <= 0
 v[4,9,8]: 991.4527540687974 beta[4,9] + gamma[4,9,8] + delta[4,9,8]
   + epsilon[4,9] <= 0
 v[4,9,9]: gamma[4,9,9] + delta[4,9,9] + epsilon[4,9] <= 0
 L3[0]: - beta[0,0] - beta[0,1] - beta[0,2] - beta[0,3] - beta[0,4]
   - beta[0,5] - beta[0,6] - beta[0,7] - beta[0,8] - beta[0,9] <= 1
 L3[1]: - beta[1,0] - beta[1,1] - beta[1,2] - beta[1,3] - beta[1,4]
   - beta[1,5] - beta[1,6] - beta[1,7] - beta[1,8] - beta[1,9] <= 1
 L3[2]: - beta[2,0] - beta[2,1] - beta[2,2] - beta[2,3] - beta[2,4]
   - beta[2,5] - beta[2,6] - beta[2,7] - beta[2,8] - beta[2,9] <= 1
 L3[3]: - beta[3,0] - beta[3,1] - beta[3,2] - beta[3,3] - beta[3,4]
   - beta[3,5] - beta[3,6] - beta[3,7] - beta[3,8] - beta[3,9] <= 1
 L3[4]: - beta[4,0] - beta[4,1] - beta[4,2] - beta[4,3] - beta[4,4]
   - beta[4,5] - beta[4,6] - beta[4,7] - beta[4,8] - beta[4,9] <= 1
Bounds
 -infinity <= beta[0,0] <= 0
 -infinity <= beta[0,1] <= 0
 -infinity <= beta[0,2] <= 0
 -infinity <= beta[0,3] <= 0
 -infinity <= beta[0,4] <= 0
 -infinity <= beta[0,5] <= 0
 -infinity <= beta[0,6] <= 0
 -infinity <= beta[0,7] <= 0
 -infinity <= beta[0,8] <= 0
 -infinity <= beta[0,9] <= 0
 -infinity <= beta[1,0] <= 0
 -infinity <= beta[1,1] <= 0
 -infinity <= beta[1,2] <= 0
 -infinity <= beta[1,3] <= 0
 -infinity <= beta[1,4] <= 0
 -infinity <= beta[1,5] <= 0
 -infinity <= beta[1,6] <= 0
 -infinity <= beta[1,7] <= 0
 -infinity <= beta[1,8] <= 0
 -infinity <= beta[1,9] <= 0
 -infinity <= beta[2,0] <= 0
 -infinity <= beta[2,1] <= 0
 -infinity <= beta[2,2] <= 0
 -infinity <= beta[2,3] <= 0
 -infinity <= beta[2,4] <= 0
 -infinity <= beta[2,5] <= 0
 -infinity <= beta[2,6] <= 0
 -infinity <= beta[2,7] <= 0
 -infinity <= beta[2,8] <= 0
 -infinity <= beta[2,9] <= 0
 -infinity <= beta[3,0] <= 0
 -infinity <= beta[3,1] <= 0
 -infinity <= beta[3,2] <= 0
 -infinity <= beta[3,3] <= 0
 -infinity <= beta[3,4] <= 0
 -infinity <= beta[3,5] <= 0
 -infinity <= beta[3,6] <= 0
 -infinity <= beta[3,7] <= 0
 -infinity <= beta[3,8] <= 0
 -infinity <= beta[3,9] <= 0
 -infinity <= beta[4,0] <= 0
 -infinity <= beta[4,1] <= 0
 -infinity <= beta[4,2] <= 0
 -infinity <= beta[4,3] <= 0
 -infinity <= beta[4,4] <= 0
 -infinity <= beta[4,5] <= 0
 -infinity <= beta[4,6] <= 0
 -infinity <= beta[4,7] <= 0
 -infinity <= beta[4,8] <= 0
 -infinity <= beta[4,9] <= 0
 -infinity <= gamma[0,0,0] <= 0
 -infinity <= gamma[0,0,1] <= 0
 -infinity <= gamma[0,0,2] <= 0
 -infinity <= gamma[0,0,3] <= 0
 -infinity <= gamma[0,0,4] <= 0
 -infinity <= gamma[0,0,5] <= 0
 -infinity <= gamma[0,0,6] <= 0
 -infinity <= gamma[0,0,7] <= 0
 -infinity <= gamma[0,0,8] <= 0
 -infinity <= gamma[0,0,9] <= 0
 -infinity <= gamma[0,1,0] <= 0
 -infinity <= gamma[0,1,1] <= 0
 -infinity <= gamma[0,1,2] <= 0
 -infinity <= gamma[0,1,3] <= 0
 -infinity <= gamma[0,1,4] <= 0
 -infinity <= gamma[0,1,5] <= 0
 -infinity <= gamma[0,1,6] <= 0
 -infinity <= gamma[0,1,7] <= 0
 -infinity <= gamma[0,1,8] <= 0
 -infinity <= gamma[0,1,9] <= 0
 -infinity <= gamma[0,2,0] <= 0
 -infinity <= gamma[0,2,1] <= 0
 -infinity <= gamma[0,2,2] <= 0
 -infinity <= gamma[0,2,3] <= 0
 -infinity <= gamma[0,2,4] <= 0
 -infinity <= gamma[0,2,5] <= 0
 -infinity <= gamma[0,2,6] <= 0
 -infinity <= gamma[0,2,7] <= 0
 -infinity <= gamma[0,2,8] <= 0
 -infinity <= gamma[0,2,9] <= 0
 -infinity <= gamma[0,3,0] <= 0
 -infinity <= gamma[0,3,1] <= 0
 -infinity <= gamma[0,3,2] <= 0
 -infinity <= gamma[0,3,3] <= 0
 -infinity <= gamma[0,3,4] <= 0
 -infinity <= gamma[0,3,5] <= 0
 -infinity <= gamma[0,3,6] <= 0
 -infinity <= gamma[0,3,7] <= 0
 -infinity <= gamma[0,3,8] <= 0
 -infinity <= gamma[0,3,9] <= 0
 -infinity <= gamma[0,4,0] <= 0
 -infinity <= gamma[0,4,1] <= 0
 -infinity <= gamma[0,4,2] <= 0
 -infinity <= gamma[0,4,3] <= 0
 -infinity <= gamma[0,4,4] <= 0
 -infinity <= gamma[0,4,5] <= 0
 -infinity <= gamma[0,4,6] <= 0
 -infinity <= gamma[0,4,7] <= 0
 -infinity <= gamma[0,4,8] <= 0
 -infinity <= gamma[0,4,9] <= 0
 -infinity <= gamma[0,5,0] <= 0
 -infinity <= gamma[0,5,1] <= 0
 -infinity <= gamma[0,5,2] <= 0
 -infinity <= gamma[0,5,3] <= 0
 -infinity <= gamma[0,5,4] <= 0
 -infinity <= gamma[0,5,5] <= 0
 -infinity <= gamma[0,5,6] <= 0
 -infinity <= gamma[0,5,7] <= 0
 -infinity <= gamma[0,5,8] <= 0
 -infinity <= gamma[0,5,9] <= 0
 -infinity <= gamma[0,6,0] <= 0
 -infinity <= gamma[0,6,1] <= 0
 -infinity <= gamma[0,6,2] <= 0
 -infinity <= gamma[0,6,3] <= 0
 -infinity <= gamma[0,6,4] <= 0
 -infinity <= gamma[0,6,5] <= 0
 -infinity <= gamma[0,6,6] <= 0
 -infinity <= gamma[0,6,7] <= 0
 -infinity <= gamma[0,6,8] <= 0
 -infinity <= gamma[0,6,9] <= 0
 -infinity <= gamma[0,7,0] <= 0
 -infinity <= gamma[0,7,1] <= 0
 -infinity <= gamma[0,7,2] <= 0
 -infinity <= gamma[0,7,3] <= 0
 -infinity <= gamma[0,7,4] <= 0
 -infinity <= gamma[0,7,5] <= 0
 -infinity <= gamma[0,7,6] <= 0
 -infinity <= gamma[0,7,7] <= 0
 -infinity <= gamma[0,7,8] <= 0
 -infinity <= gamma[0,7,9] <= 0
 -infinity <= gamma[0,8,0] <= 0
 -infinity <= gamma[0,8,1] <= 0
 -infinity <= gamma[0,8,2] <= 0
 -infinity <= gamma[0,8,3] <= 0
 -infinity <= gamma[0,8,4] <= 0
 -infinity <= gamma[0,8,5] <= 0
 -infinity <= gamma[0,8,6] <= 0
 -infinity <= gamma[0,8,7] <= 0
 -infinity <= gamma[0,8,8] <= 0
 -infinity <= gamma[0,8,9] <= 0
 -infinity <= gamma[0,9,0] <= 0
 -infinity <= gamma[0,9,1] <= 0
 -infinity <= gamma[0,9,2] <= 0
 -infinity <= gamma[0,9,3] <= 0
 -infinity <= gamma[0,9,4] <= 0
 -infinity <= gamma[0,9,5] <= 0
 -infinity <= gamma[0,9,6] <= 0
 -infinity <= gamma[0,9,7] <= 0
 -infinity <= gamma[0,9,8] <= 0
 -infinity <= gamma[0,9,9] <= 0
 -infinity <= gamma[1,0,0] <= 0
 -infinity <= gamma[1,0,1] <= 0
 -infinity <= gamma[1,0,2] <= 0
 -infinity <= gamma[1,0,3] <= 0
 -infinity <= gamma[1,0,4] <= 0
 -infinity <= gamma[1,0,5] <= 0
 -infinity <= gamma[1,0,6] <= 0
 -infinity <= gamma[1,0,7] <= 0
 -infinity <= gamma[1,0,8] <= 0
 -infinity <= gamma[1,0,9] <= 0
 -infinity <= gamma[1,1,0] <= 0
 -infinity <= gamma[1,1,1] <= 0
 -infinity <= gamma[1,1,2] <= 0
 -infinity <= gamma[1,1,3] <= 0
 -infinity <= gamma[1,1,4] <= 0
 -infinity <= gamma[1,1,5] <= 0
 -infinity <= gamma[1,1,6] <= 0
 -infinity <= gamma[1,1,7] <= 0
 -infinity <= gamma[1,1,8] <= 0
 -infinity <= gamma[1,1,9] <= 0
 -infinity <= gamma[1,2,0] <= 0
 -infinity <= gamma[1,2,1] <= 0
 -infinity <= gamma[1,2,2] <= 0
 -infinity <= gamma[1,2,3] <= 0
 -infinity <= gamma[1,2,4] <= 0
 -infinity <= gamma[1,2,5] <= 0
 -infinity <= gamma[1,2,6] <= 0
 -infinity <= gamma[1,2,7] <= 0
 -infinity <= gamma[1,2,8] <= 0
 -infinity <= gamma[1,2,9] <= 0
 -infinity <= gamma[1,3,0] <= 0
 -infinity <= gamma[1,3,1] <= 0
 -infinity <= gamma[1,3,2] <= 0
 -infinity <= gamma[1,3,3] <= 0
 -infinity <= gamma[1,3,4] <= 0
 -infinity <= gamma[1,3,5] <= 0
 -infinity <= gamma[1,3,6] <= 0
 -infinity <= gamma[1,3,7] <= 0
 -infinity <= gamma[1,3,8] <= 0
 -infinity <= gamma[1,3,9] <= 0
 -infinity <= gamma[1,4,0] <= 0
 -infinity <= gamma[1,4,1] <= 0
 -infinity <= gamma[1,4,2] <= 0
 -infinity <= gamma[1,4,3] <= 0
 -infinity <= gamma[1,4,4] <= 0
 -infinity <= gamma[1,4,5] <= 0
 -infinity <= gamma[1,4,6] <= 0
 -infinity <= gamma[1,4,7] <= 0
 -infinity <= gamma[1,4,8] <= 0
 -infinity <= gamma[1,4,9] <= 0
 -infinity <= gamma[1,5,0] <= 0
 -infinity <= gamma[1,5,1] <= 0
 -infinity <= gamma[1,5,2] <= 0
 -infinity <= gamma[1,5,3] <= 0
 -infinity <= gamma[1,5,4] <= 0
 -infinity <= gamma[1,5,5] <= 0
 -infinity <= gamma[1,5,6] <= 0
 -infinity <= gamma[1,5,7] <= 0
 -infinity <= gamma[1,5,8] <= 0
 -infinity <= gamma[1,5,9] <= 0
 -infinity <= gamma[1,6,0] <= 0
 -infinity <= gamma[1,6,1] <= 0
 -infinity <= gamma[1,6,2] <= 0
 -infinity <= gamma[1,6,3] <= 0
 -infinity <= gamma[1,6,4] <= 0
 -infinity <= gamma[1,6,5] <= 0
 -infinity <= gamma[1,6,6] <= 0
 -infinity <= gamma[1,6,7] <= 0
 -infinity <= gamma[1,6,8] <= 0
 -infinity <= gamma[1,6,9] <= 0
 -infinity <= gamma[1,7,0] <= 0
 -infinity <= gamma[1,7,1] <= 0
 -infinity <= gamma[1,7,2] <= 0
 -infinity <= gamma[1,7,3] <= 0
 -infinity <= gamma[1,7,4] <= 0
 -infinity <= gamma[1,7,5] <= 0
 -infinity <= gamma[1,7,6] <= 0
 -infinity <= gamma[1,7,7] <= 0
 -infinity <= gamma[1,7,8] <= 0
 -infinity <= gamma[1,7,9] <= 0
 -infinity <= gamma[1,8,0] <= 0
 -infinity <= gamma[1,8,1] <= 0
 -infinity <= gamma[1,8,2] <= 0
 -infinity <= gamma[1,8,3] <= 0
 -infinity <= gamma[1,8,4] <= 0
 -infinity <= gamma[1,8,5] <= 0
 -infinity <= gamma[1,8,6] <= 0
 -infinity <= gamma[1,8,7] <= 0
 -infinity <= gamma[1,8,8] <= 0
 -infinity <= gamma[1,8,9] <= 0
 -infinity <= gamma[1,9,0] <= 0
 -infinity <= gamma[1,9,1] <= 0
 -infinity <= gamma[1,9,2] <= 0
 -infinity <= gamma[1,9,3] <= 0
 -infinity <= gamma[1,9,4] <= 0
 -infinity <= gamma[1,9,5] <= 0
 -infinity <= gamma[1,9,6] <= 0
 -infinity <= gamma[1,9,7] <= 0
 -infinity <= gamma[1,9,8] <= 0
 -infinity <= gamma[1,9,9] <= 0
 -infinity <= gamma[2,0,0] <= 0
 -infinity <= gamma[2,0,1] <= 0
 -infinity <= gamma[2,0,2] <= 0
 -infinity <= gamma[2,0,3] <= 0
 -infinity <= gamma[2,0,4] <= 0
 -infinity <= gamma[2,0,5] <= 0
 -infinity <= gamma[2,0,6] <= 0
 -infinity <= gamma[2,0,7] <= 0
 -infinity <= gamma[2,0,8] <= 0
 -infinity <= gamma[2,0,9] <= 0
 -infinity <= gamma[2,1,0] <= 0
 -infinity <= gamma[2,1,1] <= 0
 -infinity <= gamma[2,1,2] <= 0
 -infinity <= gamma[2,1,3] <= 0
 -infinity <= gamma[2,1,4] <= 0
 -infinity <= gamma[2,1,5] <= 0
 -infinity <= gamma[2,1,6] <= 0
 -infinity <= gamma[2,1,7] <= 0
 -infinity <= gamma[2,1,8] <= 0
 -infinity <= gamma[2,1,9] <= 0
 -infinity <= gamma[2,2,0] <= 0
 -infinity <= gamma[2,2,1] <= 0
 -infinity <= gamma[2,2,2] <= 0
 -infinity <= gamma[2,2,3] <= 0
 -infinity <= gamma[2,2,4] <= 0
 -infinity <= gamma[2,2,5] <= 0
 -infinity <= gamma[2,2,6] <= 0
 -infinity <= gamma[2,2,7] <= 0
 -infinity <= gamma[2,2,8] <= 0
 -infinity <= gamma[2,2,9] <= 0
 -infinity <= gamma[2,3,0] <= 0
 -infinity <= gamma[2,3,1] <= 0
 -infinity <= gamma[2,3,2] <= 0
 -infinity <= gamma[2,3,3] <= 0
 -infinity <= gamma[2,3,4] <= 0
 -infinity <= gamma[2,3,5] <= 0
 -infinity <= gamma[2,3,6] <= 0
 -infinity <= gamma[2,3,7] <= 0
 -infinity <= gamma[2,3,8] <= 0
 -infinity <= gamma[2,3,9] <= 0
 -infinity <= gamma[2,4,0] <= 0
 -infinity <= gamma[2,4,1] <= 0
 -infinity <= gamma[2,4,2] <= 0
 -infinity <= gamma[2,4,3] <= 0
 -infinity <= gamma[2,4,4] <= 0
 -infinity <= gamma[2,4,5] <= 0
 -infinity <= gamma[2,4,6] <= 0
 -infinity <= gamma[2,4,7] <= 0
 -infinity <= gamma[2,4,8] <= 0
 -infinity <= gamma[2,4,9] <= 0
 -infinity <= gamma[2,5,0] <= 0
 -infinity <= gamma[2,5,1] <= 0
 -infinity <= gamma[2,5,2] <= 0
 -infinity <= gamma[2,5,3] <= 0
 -infinity <= gamma[2,5,4] <= 0
 -infinity <= gamma[2,5,5] <= 0
 -infinity <= gamma[2,5,6] <= 0
 -infinity <= gamma[2,5,7] <= 0
 -infinity <= gamma[2,5,8] <= 0
 -infinity <= gamma[2,5,9] <= 0
 -infinity <= gamma[2,6,0] <= 0
 -infinity <= gamma[2,6,1] <= 0
 -infinity <= gamma[2,6,2] <= 0
 -infinity <= gamma[2,6,3] <= 0
 -infinity <= gamma[2,6,4] <= 0
 -infinity <= gamma[2,6,5] <= 0
 -infinity <= gamma[2,6,6] <= 0
 -infinity <= gamma[2,6,7] <= 0
 -infinity <= gamma[2,6,8] <= 0
 -infinity <= gamma[2,6,9] <= 0
 -infinity <= gamma[2,7,0] <= 0
 -infinity <= gamma[2,7,1] <= 0
 -infinity <= gamma[2,7,2] <= 0
 -infinity <= gamma[2,7,3] <= 0
 -infinity <= gamma[2,7,4] <= 0
 -infinity <= gamma[2,7,5] <= 0
 -infinity <= gamma[2,7,6] <= 0
 -infinity <= gamma[2,7,7] <= 0
 -infinity <= gamma[2,7,8] <= 0
 -infinity <= gamma[2,7,9] <= 0
 -infinity <= gamma[2,8,0] <= 0
 -infinity <= gamma[2,8,1] <= 0
 -infinity <= gamma[2,8,2] <= 0
 -infinity <= gamma[2,8,3] <= 0
 -infinity <= gamma[2,8,4] <= 0
 -infinity <= gamma[2,8,5] <= 0
 -infinity <= gamma[2,8,6] <= 0
 -infinity <= gamma[2,8,7] <= 0
 -infinity <= gamma[2,8,8] <= 0
 -infinity <= gamma[2,8,9] <= 0
 -infinity <= gamma[2,9,0] <= 0
 -infinity <= gamma[2,9,1] <= 0
 -infinity <= gamma[2,9,2] <= 0
 -infinity <= gamma[2,9,3] <= 0
 -infinity <= gamma[2,9,4] <= 0
 -infinity <= gamma[2,9,5] <= 0
 -infinity <= gamma[2,9,6] <= 0
 -infinity <= gamma[2,9,7] <= 0
 -infinity <= gamma[2,9,8] <= 0
 -infinity <= gamma[2,9,9] <= 0
 -infinity <= gamma[3,0,0] <= 0
 -infinity <= gamma[3,0,1] <= 0
 -infinity <= gamma[3,0,2] <= 0
 -infinity <= gamma[3,0,3] <= 0
 -infinity <= gamma[3,0,4] <= 0
 -infinity <= gamma[3,0,5] <= 0
 -infinity <= gamma[3,0,6] <= 0
 -infinity <= gamma[3,0,7] <= 0
 -infinity <= gamma[3,0,8] <= 0
 -infinity <= gamma[3,0,9] <= 0
 -infinity <= gamma[3,1,0] <= 0
 -infinity <= gamma[3,1,1] <= 0
 -infinity <= gamma[3,1,2] <= 0
 -infinity <= gamma[3,1,3] <= 0
 -infinity <= gamma[3,1,4] <= 0
 -infinity <= gamma[3,1,5] <= 0
 -infinity <= gamma[3,1,6] <= 0
 -infinity <= gamma[3,1,7] <= 0
 -infinity <= gamma[3,1,8] <= 0
 -infinity <= gamma[3,1,9] <= 0
 -infinity <= gamma[3,2,0] <= 0
 -infinity <= gamma[3,2,1] <= 0
 -infinity <= gamma[3,2,2] <= 0
 -infinity <= gamma[3,2,3] <= 0
 -infinity <= gamma[3,2,4] <= 0
 -infinity <= gamma[3,2,5] <= 0
 -infinity <= gamma[3,2,6] <= 0
 -infinity <= gamma[3,2,7] <= 0
 -infinity <= gamma[3,2,8] <= 0
 -infinity <= gamma[3,2,9] <= 0
 -infinity <= gamma[3,3,0] <= 0
 -infinity <= gamma[3,3,1] <= 0
 -infinity <= gamma[3,3,2] <= 0
 -infinity <= gamma[3,3,3] <= 0
 -infinity <= gamma[3,3,4] <= 0
 -infinity <= gamma[3,3,5] <= 0
 -infinity <= gamma[3,3,6] <= 0
 -infinity <= gamma[3,3,7] <= 0
 -infinity <= gamma[3,3,8] <= 0
 -infinity <= gamma[3,3,9] <= 0
 -infinity <= gamma[3,4,0] <= 0
 -infinity <= gamma[3,4,1] <= 0
 -infinity <= gamma[3,4,2] <= 0
 -infinity <= gamma[3,4,3] <= 0
 -infinity <= gamma[3,4,4] <= 0
 -infinity <= gamma[3,4,5] <= 0
 -infinity <= gamma[3,4,6] <= 0
 -infinity <= gamma[3,4,7] <= 0
 -infinity <= gamma[3,4,8] <= 0
 -infinity <= gamma[3,4,9] <= 0
 -infinity <= gamma[3,5,0] <= 0
 -infinity <= gamma[3,5,1] <= 0
 -infinity <= gamma[3,5,2] <= 0
 -infinity <= gamma[3,5,3] <= 0
 -infinity <= gamma[3,5,4] <= 0
 -infinity <= gamma[3,5,5] <= 0
 -infinity <= gamma[3,5,6] <= 0
 -infinity <= gamma[3,5,7] <= 0
 -infinity <= gamma[3,5,8] <= 0
 -infinity <= gamma[3,5,9] <= 0
 -infinity <= gamma[3,6,0] <= 0
 -infinity <= gamma[3,6,1] <= 0
 -infinity <= gamma[3,6,2] <= 0
 -infinity <= gamma[3,6,3] <= 0
 -infinity <= gamma[3,6,4] <= 0
 -infinity <= gamma[3,6,5] <= 0
 -infinity <= gamma[3,6,6] <= 0
 -infinity <= gamma[3,6,7] <= 0
 -infinity <= gamma[3,6,8] <= 0
 -infinity <= gamma[3,6,9] <= 0
 -infinity <= gamma[3,7,0] <= 0
 -infinity <= gamma[3,7,1] <= 0
 -infinity <= gamma[3,7,2] <= 0
 -infinity <= gamma[3,7,3] <= 0
 -infinity <= gamma[3,7,4] <= 0
 -infinity <= gamma[3,7,5] <= 0
 -infinity <= gamma[3,7,6] <= 0
 -infinity <= gamma[3,7,7] <= 0
 -infinity <= gamma[3,7,8] <= 0
 -infinity <= gamma[3,7,9] <= 0
 -infinity <= gamma[3,8,0] <= 0
 -infinity <= gamma[3,8,1] <= 0
 -infinity <= gamma[3,8,2] <= 0
 -infinity <= gamma[3,8,3] <= 0
 -infinity <= gamma[3,8,4] <= 0
 -infinity <= gamma[3,8,5] <= 0
 -infinity <= gamma[3,8,6] <= 0
 -infinity <= gamma[3,8,7] <= 0
 -infinity <= gamma[3,8,8] <= 0
 -infinity <= gamma[3,8,9] <= 0
 -infinity <= gamma[3,9,0] <= 0
 -infinity <= gamma[3,9,1] <= 0
 -infinity <= gamma[3,9,2] <= 0
 -infinity <= gamma[3,9,3] <= 0
 -infinity <= gamma[3,9,4] <= 0
 -infinity <= gamma[3,9,5] <= 0
 -infinity <= gamma[3,9,6] <= 0
 -infinity <= gamma[3,9,7] <= 0
 -infinity <= gamma[3,9,8] <= 0
 -infinity <= gamma[3,9,9] <= 0
 -infinity <= gamma[4,0,0] <= 0
 -infinity <= gamma[4,0,1] <= 0
 -infinity <= gamma[4,0,2] <= 0
 -infinity <= gamma[4,0,3] <= 0
 -infinity <= gamma[4,0,4] <= 0
 -infinity <= gamma[4,0,5] <= 0
 -infinity <= gamma[4,0,6] <= 0
 -infinity <= gamma[4,0,7] <= 0
 -infinity <= gamma[4,0,8] <= 0
 -infinity <= gamma[4,0,9] <= 0
 -infinity <= gamma[4,1,0] <= 0
 -infinity <= gamma[4,1,1] <= 0
 -infinity <= gamma[4,1,2] <= 0
 -infinity <= gamma[4,1,3] <= 0
 -infinity <= gamma[4,1,4] <= 0
 -infinity <= gamma[4,1,5] <= 0
 -infinity <= gamma[4,1,6] <= 0
 -infinity <= gamma[4,1,7] <= 0
 -infinity <= gamma[4,1,8] <= 0
 -infinity <= gamma[4,1,9] <= 0
 -infinity <= gamma[4,2,0] <= 0
 -infinity <= gamma[4,2,1] <= 0
 -infinity <= gamma[4,2,2] <= 0
 -infinity <= gamma[4,2,3] <= 0
 -infinity <= gamma[4,2,4] <= 0
 -infinity <= gamma[4,2,5] <= 0
 -infinity <= gamma[4,2,6] <= 0
 -infinity <= gamma[4,2,7] <= 0
 -infinity <= gamma[4,2,8] <= 0
 -infinity <= gamma[4,2,9] <= 0
 -infinity <= gamma[4,3,0] <= 0
 -infinity <= gamma[4,3,1] <= 0
 -infinity <= gamma[4,3,2] <= 0
 -infinity <= gamma[4,3,3] <= 0
 -infinity <= gamma[4,3,4] <= 0
 -infinity <= gamma[4,3,5] <= 0
 -infinity <= gamma[4,3,6] <= 0
 -infinity <= gamma[4,3,7] <= 0
 -infinity <= gamma[4,3,8] <= 0
 -infinity <= gamma[4,3,9] <= 0
 -infinity <= gamma[4,4,0] <= 0
 -infinity <= gamma[4,4,1] <= 0
 -infinity <= gamma[4,4,2] <= 0
 -infinity <= gamma[4,4,3] <= 0
 -infinity <= gamma[4,4,4] <= 0
 -infinity <= gamma[4,4,5] <= 0
 -infinity <= gamma[4,4,6] <= 0
 -infinity <= gamma[4,4,7] <= 0
 -infinity <= gamma[4,4,8] <= 0
 -infinity <= gamma[4,4,9] <= 0
 -infinity <= gamma[4,5,0] <= 0
 -infinity <= gamma[4,5,1] <= 0
 -infinity <= gamma[4,5,2] <= 0
 -infinity <= gamma[4,5,3] <= 0
 -infinity <= gamma[4,5,4] <= 0
 -infinity <= gamma[4,5,5] <= 0
 -infinity <= gamma[4,5,6] <= 0
 -infinity <= gamma[4,5,7] <= 0
 -infinity <= gamma[4,5,8] <= 0
 -infinity <= gamma[4,5,9] <= 0
 -infinity <= gamma[4,6,0] <= 0
 -infinity <= gamma[4,6,1] <= 0
 -infinity <= gamma[4,6,2] <= 0
 -infinity <= gamma[4,6,3] <= 0
 -infinity <= gamma[4,6,4] <= 0
 -infinity <= gamma[4,6,5] <= 0
 -infinity <= gamma[4,6,6] <= 0
 -infinity <= gamma[4,6,7] <= 0
 -infinity <= gamma[4,6,8] <= 0
 -infinity <= gamma[4,6,9] <= 0
 -infinity <= gamma[4,7,0] <= 0
 -infinity <= gamma[4,7,1] <= 0
 -infinity <= gamma[4,7,2] <= 0
 -infinity <= gamma[4,7,3] <= 0
 -infinity <= gamma[4,7,4] <= 0
 -infinity <= gamma[4,7,5] <= 0
 -infinity <= gamma[4,7,6] <= 0
 -infinity <= gamma[4,7,7] <= 0
 -infinity <= gamma[4,7,8] <= 0
 -infinity <= gamma[4,7,9] <= 0
 -infinity <= gamma[4,8,0] <= 0
 -infinity <= gamma[4,8,1] <= 0
 -infinity <= gamma[4,8,2] <= 0
 -infinity <= gamma[4,8,3] <= 0
 -infinity <= gamma[4,8,4] <= 0
 -infinity <= gamma[4,8,5] <= 0
 -infinity <= gamma[4,8,6] <= 0
 -infinity <= gamma[4,8,7] <= 0
 -infinity <= gamma[4,8,8] <= 0
 -infinity <= gamma[4,8,9] <= 0
 -infinity <= gamma[4,9,0] <= 0
 -infinity <= gamma[4,9,1] <= 0
 -infinity <= gamma[4,9,2] <= 0
 -infinity <= gamma[4,9,3] <= 0
 -infinity <= gamma[4,9,4] <= 0
 -infinity <= gamma[4,9,5] <= 0
 -infinity <= gamma[4,9,6] <= 0
 -infinity <= gamma[4,9,7] <= 0
 -infinity <= gamma[4,9,8] <= 0
 -infinity <= gamma[4,9,9] <= 0
 -infinity <= delta[0,0,0] <= 0
 -infinity <= delta[0,0,1] <= 0
 -infinity <= delta[0,0,2] <= 0
 -infinity <= delta[0,0,3] <= 0
 -infinity <= delta[0,0,4] <= 0
 -infinity <= delta[0,0,5] <= 0
 -infinity <= delta[0,0,6] <= 0
 -infinity <= delta[0,0,7] <= 0
 -infinity <= delta[0,0,8] <= 0
 -infinity <= delta[0,0,9] <= 0
 -infinity <= delta[0,1,0] <= 0
 -infinity <= delta[0,1,1] <= 0
 -infinity <= delta[0,1,2] <= 0
 -infinity <= delta[0,1,3] <= 0
 -infinity <= delta[0,1,4] <= 0
 -infinity <= delta[0,1,5] <= 0
 -infinity <= delta[0,1,6] <= 0
 -infinity <= delta[0,1,7] <= 0
 -infinity <= delta[0,1,8] <= 0
 -infinity <= delta[0,1,9] <= 0
 -infinity <= delta[0,2,0] <= 0
 -infinity <= delta[0,2,1] <= 0
 -infinity <= delta[0,2,2] <= 0
 -infinity <= delta[0,2,3] <= 0
 -infinity <= delta[0,2,4] <= 0
 -infinity <= delta[0,2,5] <= 0
 -infinity <= delta[0,2,6] <= 0
 -infinity <= delta[0,2,7] <= 0
 -infinity <= delta[0,2,8] <= 0
 -infinity <= delta[0,2,9] <= 0
 -infinity <= delta[0,3,0] <= 0
 -infinity <= delta[0,3,1] <= 0
 -infinity <= delta[0,3,2] <= 0
 -infinity <= delta[0,3,3] <= 0
 -infinity <= delta[0,3,4] <= 0
 -infinity <= delta[0,3,5] <= 0
 -infinity <= delta[0,3,6] <= 0
 -infinity <= delta[0,3,7] <= 0
 -infinity <= delta[0,3,8] <= 0
 -infinity <= delta[0,3,9] <= 0
 -infinity <= delta[0,4,0] <= 0
 -infinity <= delta[0,4,1] <= 0
 -infinity <= delta[0,4,2] <= 0
 -infinity <= delta[0,4,3] <= 0
 -infinity <= delta[0,4,4] <= 0
 -infinity <= delta[0,4,5] <= 0
 -infinity <= delta[0,4,6] <= 0
 -infinity <= delta[0,4,7] <= 0
 -infinity <= delta[0,4,8] <= 0
 -infinity <= delta[0,4,9] <= 0
 -infinity <= delta[0,5,0] <= 0
 -infinity <= delta[0,5,1] <= 0
 -infinity <= delta[0,5,2] <= 0
 -infinity <= delta[0,5,3] <= 0
 -infinity <= delta[0,5,4] <= 0
 -infinity <= delta[0,5,5] <= 0
 -infinity <= delta[0,5,6] <= 0
 -infinity <= delta[0,5,7] <= 0
 -infinity <= delta[0,5,8] <= 0
 -infinity <= delta[0,5,9] <= 0
 -infinity <= delta[0,6,0] <= 0
 -infinity <= delta[0,6,1] <= 0
 -infinity <= delta[0,6,2] <= 0
 -infinity <= delta[0,6,3] <= 0
 -infinity <= delta[0,6,4] <= 0
 -infinity <= delta[0,6,5] <= 0
 -infinity <= delta[0,6,6] <= 0
 -infinity <= delta[0,6,7] <= 0
 -infinity <= delta[0,6,8] <= 0
 -infinity <= delta[0,6,9] <= 0
 -infinity <= delta[0,7,0] <= 0
 -infinity <= delta[0,7,1] <= 0
 -infinity <= delta[0,7,2] <= 0
 -infinity <= delta[0,7,3] <= 0
 -infinity <= delta[0,7,4] <= 0
 -infinity <= delta[0,7,5] <= 0
 -infinity <= delta[0,7,6] <= 0
 -infinity <= delta[0,7,7] <= 0
 -infinity <= delta[0,7,8] <= 0
 -infinity <= delta[0,7,9] <= 0
 -infinity <= delta[0,8,0] <= 0
 -infinity <= delta[0,8,1] <= 0
 -infinity <= delta[0,8,2] <= 0
 -infinity <= delta[0,8,3] <= 0
 -infinity <= delta[0,8,4] <= 0
 -infinity <= delta[0,8,5] <= 0
 -infinity <= delta[0,8,6] <= 0
 -infinity <= delta[0,8,7] <= 0
 -infinity <= delta[0,8,8] <= 0
 -infinity <= delta[0,8,9] <= 0
 -infinity <= delta[0,9,0] <= 0
 -infinity <= delta[0,9,1] <= 0
 -infinity <= delta[0,9,2] <= 0
 -infinity <= delta[0,9,3] <= 0
 -infinity <= delta[0,9,4] <= 0
 -infinity <= delta[0,9,5] <= 0
 -infinity <= delta[0,9,6] <= 0
 -infinity <= delta[0,9,7] <= 0
 -infinity <= delta[0,9,8] <= 0
 -infinity <= delta[0,9,9] <= 0
 -infinity <= delta[1,0,0] <= 0
 -infinity <= delta[1,0,1] <= 0
 -infinity <= delta[1,0,2] <= 0
 -infinity <= delta[1,0,3] <= 0
 -infinity <= delta[1,0,4] <= 0
 -infinity <= delta[1,0,5] <= 0
 -infinity <= delta[1,0,6] <= 0
 -infinity <= delta[1,0,7] <= 0
 -infinity <= delta[1,0,8] <= 0
 -infinity <= delta[1,0,9] <= 0
 -infinity <= delta[1,1,0] <= 0
 -infinity <= delta[1,1,1] <= 0
 -infinity <= delta[1,1,2] <= 0
 -infinity <= delta[1,1,3] <= 0
 -infinity <= delta[1,1,4] <= 0
 -infinity <= delta[1,1,5] <= 0
 -infinity <= delta[1,1,6] <= 0
 -infinity <= delta[1,1,7] <= 0
 -infinity <= delta[1,1,8] <= 0
 -infinity <= delta[1,1,9] <= 0
 -infinity <= delta[1,2,0] <= 0
 -infinity <= delta[1,2,1] <= 0
 -infinity <= delta[1,2,2] <= 0
 -infinity <= delta[1,2,3] <= 0
 -infinity <= delta[1,2,4] <= 0
 -infinity <= delta[1,2,5] <= 0
 -infinity <= delta[1,2,6] <= 0
 -infinity <= delta[1,2,7] <= 0
 -infinity <= delta[1,2,8] <= 0
 -infinity <= delta[1,2,9] <= 0
 -infinity <= delta[1,3,0] <= 0
 -infinity <= delta[1,3,1] <= 0
 -infinity <= delta[1,3,2] <= 0
 -infinity <= delta[1,3,3] <= 0
 -infinity <= delta[1,3,4] <= 0
 -infinity <= delta[1,3,5] <= 0
 -infinity <= delta[1,3,6] <= 0
 -infinity <= delta[1,3,7] <= 0
 -infinity <= delta[1,3,8] <= 0
 -infinity <= delta[1,3,9] <= 0
 -infinity <= delta[1,4,0] <= 0
 -infinity <= delta[1,4,1] <= 0
 -infinity <= delta[1,4,2] <= 0
 -infinity <= delta[1,4,3] <= 0
 -infinity <= delta[1,4,4] <= 0
 -infinity <= delta[1,4,5] <= 0
 -infinity <= delta[1,4,6] <= 0
 -infinity <= delta[1,4,7] <= 0
 -infinity <= delta[1,4,8] <= 0
 -infinity <= delta[1,4,9] <= 0
 -infinity <= delta[1,5,0] <= 0
 -infinity <= delta[1,5,1] <= 0
 -infinity <= delta[1,5,2] <= 0
 -infinity <= delta[1,5,3] <= 0
 -infinity <= delta[1,5,4] <= 0
 -infinity <= delta[1,5,5] <= 0
 -infinity <= delta[1,5,6] <= 0
 -infinity <= delta[1,5,7] <= 0
 -infinity <= delta[1,5,8] <= 0
 -infinity <= delta[1,5,9] <= 0
 -infinity <= delta[1,6,0] <= 0
 -infinity <= delta[1,6,1] <= 0
 -infinity <= delta[1,6,2] <= 0
 -infinity <= delta[1,6,3] <= 0
 -infinity <= delta[1,6,4] <= 0
 -infinity <= delta[1,6,5] <= 0
 -infinity <= delta[1,6,6] <= 0
 -infinity <= delta[1,6,7] <= 0
 -infinity <= delta[1,6,8] <= 0
 -infinity <= delta[1,6,9] <= 0
 -infinity <= delta[1,7,0] <= 0
 -infinity <= delta[1,7,1] <= 0
 -infinity <= delta[1,7,2] <= 0
 -infinity <= delta[1,7,3] <= 0
 -infinity <= delta[1,7,4] <= 0
 -infinity <= delta[1,7,5] <= 0
 -infinity <= delta[1,7,6] <= 0
 -infinity <= delta[1,7,7] <= 0
 -infinity <= delta[1,7,8] <= 0
 -infinity <= delta[1,7,9] <= 0
 -infinity <= delta[1,8,0] <= 0
 -infinity <= delta[1,8,1] <= 0
 -infinity <= delta[1,8,2] <= 0
 -infinity <= delta[1,8,3] <= 0
 -infinity <= delta[1,8,4] <= 0
 -infinity <= delta[1,8,5] <= 0
 -infinity <= delta[1,8,6] <= 0
 -infinity <= delta[1,8,7] <= 0
 -infinity <= delta[1,8,8] <= 0
 -infinity <= delta[1,8,9] <= 0
 -infinity <= delta[1,9,0] <= 0
 -infinity <= delta[1,9,1] <= 0
 -infinity <= delta[1,9,2] <= 0
 -infinity <= delta[1,9,3] <= 0
 -infinity <= delta[1,9,4] <= 0
 -infinity <= delta[1,9,5] <= 0
 -infinity <= delta[1,9,6] <= 0
 -infinity <= delta[1,9,7] <= 0
 -infinity <= delta[1,9,8] <= 0
 -infinity <= delta[1,9,9] <= 0
 -infinity <= delta[2,0,0] <= 0
 -infinity <= delta[2,0,1] <= 0
 -infinity <= delta[2,0,2] <= 0
 -infinity <= delta[2,0,3] <= 0
 -infinity <= delta[2,0,4] <= 0
 -infinity <= delta[2,0,5] <= 0
 -infinity <= delta[2,0,6] <= 0
 -infinity <= delta[2,0,7] <= 0
 -infinity <= delta[2,0,8] <= 0
 -infinity <= delta[2,0,9] <= 0
 -infinity <= delta[2,1,0] <= 0
 -infinity <= delta[2,1,1] <= 0
 -infinity <= delta[2,1,2] <= 0
 -infinity <= delta[2,1,3] <= 0
 -infinity <= delta[2,1,4] <= 0
 -infinity <= delta[2,1,5] <= 0
 -infinity <= delta[2,1,6] <= 0
 -infinity <= delta[2,1,7] <= 0
 -infinity <= delta[2,1,8] <= 0
 -infinity <= delta[2,1,9] <= 0
 -infinity <= delta[2,2,0] <= 0
 -infinity <= delta[2,2,1] <= 0
 -infinity <= delta[2,2,2] <= 0
 -infinity <= delta[2,2,3] <= 0
 -infinity <= delta[2,2,4] <= 0
 -infinity <= delta[2,2,5] <= 0
 -infinity <= delta[2,2,6] <= 0
 -infinity <= delta[2,2,7] <= 0
 -infinity <= delta[2,2,8] <= 0
 -infinity <= delta[2,2,9] <= 0
 -infinity <= delta[2,3,0] <= 0
 -infinity <= delta[2,3,1] <= 0
 -infinity <= delta[2,3,2] <= 0
 -infinity <= delta[2,3,3] <= 0
 -infinity <= delta[2,3,4] <= 0
 -infinity <= delta[2,3,5] <= 0
 -infinity <= delta[2,3,6] <= 0
 -infinity <= delta[2,3,7] <= 0
 -infinity <= delta[2,3,8] <= 0
 -infinity <= delta[2,3,9] <= 0
 -infinity <= delta[2,4,0] <= 0
 -infinity <= delta[2,4,1] <= 0
 -infinity <= delta[2,4,2] <= 0
 -infinity <= delta[2,4,3] <= 0
 -infinity <= delta[2,4,4] <= 0
 -infinity <= delta[2,4,5] <= 0
 -infinity <= delta[2,4,6] <= 0
 -infinity <= delta[2,4,7] <= 0
 -infinity <= delta[2,4,8] <= 0
 -infinity <= delta[2,4,9] <= 0
 -infinity <= delta[2,5,0] <= 0
 -infinity <= delta[2,5,1] <= 0
 -infinity <= delta[2,5,2] <= 0
 -infinity <= delta[2,5,3] <= 0
 -infinity <= delta[2,5,4] <= 0
 -infinity <= delta[2,5,5] <= 0
 -infinity <= delta[2,5,6] <= 0
 -infinity <= delta[2,5,7] <= 0
 -infinity <= delta[2,5,8] <= 0
 -infinity <= delta[2,5,9] <= 0
 -infinity <= delta[2,6,0] <= 0
 -infinity <= delta[2,6,1] <= 0
 -infinity <= delta[2,6,2] <= 0
 -infinity <= delta[2,6,3] <= 0
 -infinity <= delta[2,6,4] <= 0
 -infinity <= delta[2,6,5] <= 0
 -infinity <= delta[2,6,6] <= 0
 -infinity <= delta[2,6,7] <= 0
 -infinity <= delta[2,6,8] <= 0
 -infinity <= delta[2,6,9] <= 0
 -infinity <= delta[2,7,0] <= 0
 -infinity <= delta[2,7,1] <= 0
 -infinity <= delta[2,7,2] <= 0
 -infinity <= delta[2,7,3] <= 0
 -infinity <= delta[2,7,4] <= 0
 -infinity <= delta[2,7,5] <= 0
 -infinity <= delta[2,7,6] <= 0
 -infinity <= delta[2,7,7] <= 0
 -infinity <= delta[2,7,8] <= 0
 -infinity <= delta[2,7,9] <= 0
 -infinity <= delta[2,8,0] <= 0
 -infinity <= delta[2,8,1] <= 0
 -infinity <= delta[2,8,2] <= 0
 -infinity <= delta[2,8,3] <= 0
 -infinity <= delta[2,8,4] <= 0
 -infinity <= delta[2,8,5] <= 0
 -infinity <= delta[2,8,6] <= 0
 -infinity <= delta[2,8,7] <= 0
 -infinity <= delta[2,8,8] <= 0
 -infinity <= delta[2,8,9] <= 0
 -infinity <= delta[2,9,0] <= 0
 -infinity <= delta[2,9,1] <= 0
 -infinity <= delta[2,9,2] <= 0
 -infinity <= delta[2,9,3] <= 0
 -infinity <= delta[2,9,4] <= 0
 -infinity <= delta[2,9,5] <= 0
 -infinity <= delta[2,9,6] <= 0
 -infinity <= delta[2,9,7] <= 0
 -infinity <= delta[2,9,8] <= 0
 -infinity <= delta[2,9,9] <= 0
 -infinity <= delta[3,0,0] <= 0
 -infinity <= delta[3,0,1] <= 0
 -infinity <= delta[3,0,2] <= 0
 -infinity <= delta[3,0,3] <= 0
 -infinity <= delta[3,0,4] <= 0
 -infinity <= delta[3,0,5] <= 0
 -infinity <= delta[3,0,6] <= 0
 -infinity <= delta[3,0,7] <= 0
 -infinity <= delta[3,0,8] <= 0
 -infinity <= delta[3,0,9] <= 0
 -infinity <= delta[3,1,0] <= 0
 -infinity <= delta[3,1,1] <= 0
 -infinity <= delta[3,1,2] <= 0
 -infinity <= delta[3,1,3] <= 0
 -infinity <= delta[3,1,4] <= 0
 -infinity <= delta[3,1,5] <= 0
 -infinity <= delta[3,1,6] <= 0
 -infinity <= delta[3,1,7] <= 0
 -infinity <= delta[3,1,8] <= 0
 -infinity <= delta[3,1,9] <= 0
 -infinity <= delta[3,2,0] <= 0
 -infinity <= delta[3,2,1] <= 0
 -infinity <= delta[3,2,2] <= 0
 -infinity <= delta[3,2,3] <= 0
 -infinity <= delta[3,2,4] <= 0
 -infinity <= delta[3,2,5] <= 0
 -infinity <= delta[3,2,6] <= 0
 -infinity <= delta[3,2,7] <= 0
 -infinity <= delta[3,2,8] <= 0
 -infinity <= delta[3,2,9] <= 0
 -infinity <= delta[3,3,0] <= 0
 -infinity <= delta[3,3,1] <= 0
 -infinity <= delta[3,3,2] <= 0
 -infinity <= delta[3,3,3] <= 0
 -infinity <= delta[3,3,4] <= 0
 -infinity <= delta[3,3,5] <= 0
 -infinity <= delta[3,3,6] <= 0
 -infinity <= delta[3,3,7] <= 0
 -infinity <= delta[3,3,8] <= 0
 -infinity <= delta[3,3,9] <= 0
 -infinity <= delta[3,4,0] <= 0
 -infinity <= delta[3,4,1] <= 0
 -infinity <= delta[3,4,2] <= 0
 -infinity <= delta[3,4,3] <= 0
 -infinity <= delta[3,4,4] <= 0
 -infinity <= delta[3,4,5] <= 0
 -infinity <= delta[3,4,6] <= 0
 -infinity <= delta[3,4,7] <= 0
 -infinity <= delta[3,4,8] <= 0
 -infinity <= delta[3,4,9] <= 0
 -infinity <= delta[3,5,0] <= 0
 -infinity <= delta[3,5,1] <= 0
 -infinity <= delta[3,5,2] <= 0
 -infinity <= delta[3,5,3] <= 0
 -infinity <= delta[3,5,4] <= 0
 -infinity <= delta[3,5,5] <= 0
 -infinity <= delta[3,5,6] <= 0
 -infinity <= delta[3,5,7] <= 0
 -infinity <= delta[3,5,8] <= 0
 -infinity <= delta[3,5,9] <= 0
 -infinity <= delta[3,6,0] <= 0
 -infinity <= delta[3,6,1] <= 0
 -infinity <= delta[3,6,2] <= 0
 -infinity <= delta[3,6,3] <= 0
 -infinity <= delta[3,6,4] <= 0
 -infinity <= delta[3,6,5] <= 0
 -infinity <= delta[3,6,6] <= 0
 -infinity <= delta[3,6,7] <= 0
 -infinity <= delta[3,6,8] <= 0
 -infinity <= delta[3,6,9] <= 0
 -infinity <= delta[3,7,0] <= 0
 -infinity <= delta[3,7,1] <= 0
 -infinity <= delta[3,7,2] <= 0
 -infinity <= delta[3,7,3] <= 0
 -infinity <= delta[3,7,4] <= 0
 -infinity <= delta[3,7,5] <= 0
 -infinity <= delta[3,7,6] <= 0
 -infinity <= delta[3,7,7] <= 0
 -infinity <= delta[3,7,8] <= 0
 -infinity <= delta[3,7,9] <= 0
 -infinity <= delta[3,8,0] <= 0
 -infinity <= delta[3,8,1] <= 0
 -infinity <= delta[3,8,2] <= 0
 -infinity <= delta[3,8,3] <= 0
 -infinity <= delta[3,8,4] <= 0
 -infinity <= delta[3,8,5] <= 0
 -infinity <= delta[3,8,6] <= 0
 -infinity <= delta[3,8,7] <= 0
 -infinity <= delta[3,8,8] <= 0
 -infinity <= delta[3,8,9] <= 0
 -infinity <= delta[3,9,0] <= 0
 -infinity <= delta[3,9,1] <= 0
 -infinity <= delta[3,9,2] <= 0
 -infinity <= delta[3,9,3] <= 0
 -infinity <= delta[3,9,4] <= 0
 -infinity <= delta[3,9,5] <= 0
 -infinity <= delta[3,9,6] <= 0
 -infinity <= delta[3,9,7] <= 0
 -infinity <= delta[3,9,8] <= 0
 -infinity <= delta[3,9,9] <= 0
 -infinity <= delta[4,0,0] <= 0
 -infinity <= delta[4,0,1] <= 0
 -infinity <= delta[4,0,2] <= 0
 -infinity <= delta[4,0,3] <= 0
 -infinity <= delta[4,0,4] <= 0
 -infinity <= delta[4,0,5] <= 0
 -infinity <= delta[4,0,6] <= 0
 -infinity <= delta[4,0,7] <= 0
 -infinity <= delta[4,0,8] <= 0
 -infinity <= delta[4,0,9] <= 0
 -infinity <= delta[4,1,0] <= 0
 -infinity <= delta[4,1,1] <= 0
 -infinity <= delta[4,1,2] <= 0
 -infinity <= delta[4,1,3] <= 0
 -infinity <= delta[4,1,4] <= 0
 -infinity <= delta[4,1,5] <= 0
 -infinity <= delta[4,1,6] <= 0
 -infinity <= delta[4,1,7] <= 0
 -infinity <= delta[4,1,8] <= 0
 -infinity <= delta[4,1,9] <= 0
 -infinity <= delta[4,2,0] <= 0
 -infinity <= delta[4,2,1] <= 0
 -infinity <= delta[4,2,2] <= 0
 -infinity <= delta[4,2,3] <= 0
 -infinity <= delta[4,2,4] <= 0
 -infinity <= delta[4,2,5] <= 0
 -infinity <= delta[4,2,6] <= 0
 -infinity <= delta[4,2,7] <= 0
 -infinity <= delta[4,2,8] <= 0
 -infinity <= delta[4,2,9] <= 0
 -infinity <= delta[4,3,0] <= 0
 -infinity <= delta[4,3,1] <= 0
 -infinity <= delta[4,3,2] <= 0
 -infinity <= delta[4,3,3] <= 0
 -infinity <= delta[4,3,4] <= 0
 -infinity <= delta[4,3,5] <= 0
 -infinity <= delta[4,3,6] <= 0
 -infinity <= delta[4,3,7] <= 0
 -infinity <= delta[4,3,8] <= 0
 -infinity <= delta[4,3,9] <= 0
 -infinity <= delta[4,4,0] <= 0
 -infinity <= delta[4,4,1] <= 0
 -infinity <= delta[4,4,2] <= 0
 -infinity <= delta[4,4,3] <= 0
 -infinity <= delta[4,4,4] <= 0
 -infinity <= delta[4,4,5] <= 0
 -infinity <= delta[4,4,6] <= 0
 -infinity <= delta[4,4,7] <= 0
 -infinity <= delta[4,4,8] <= 0
 -infinity <= delta[4,4,9] <= 0
 -infinity <= delta[4,5,0] <= 0
 -infinity <= delta[4,5,1] <= 0
 -infinity <= delta[4,5,2] <= 0
 -infinity <= delta[4,5,3] <= 0
 -infinity <= delta[4,5,4] <= 0
 -infinity <= delta[4,5,5] <= 0
 -infinity <= delta[4,5,6] <= 0
 -infinity <= delta[4,5,7] <= 0
 -infinity <= delta[4,5,8] <= 0
 -infinity <= delta[4,5,9] <= 0
 -infinity <= delta[4,6,0] <= 0
 -infinity <= delta[4,6,1] <= 0
 -infinity <= delta[4,6,2] <= 0
 -infinity <= delta[4,6,3] <= 0
 -infinity <= delta[4,6,4] <= 0
 -infinity <= delta[4,6,5] <= 0
 -infinity <= delta[4,6,6] <= 0
 -infinity <= delta[4,6,7] <= 0
 -infinity <= delta[4,6,8] <= 0
 -infinity <= delta[4,6,9] <= 0
 -infinity <= delta[4,7,0] <= 0
 -infinity <= delta[4,7,1] <= 0
 -infinity <= delta[4,7,2] <= 0
 -infinity <= delta[4,7,3] <= 0
 -infinity <= delta[4,7,4] <= 0
 -infinity <= delta[4,7,5] <= 0
 -infinity <= delta[4,7,6] <= 0
 -infinity <= delta[4,7,7] <= 0
 -infinity <= delta[4,7,8] <= 0
 -infinity <= delta[4,7,9] <= 0
 -infinity <= delta[4,8,0] <= 0
 -infinity <= delta[4,8,1] <= 0
 -infinity <= delta[4,8,2] <= 0
 -infinity <= delta[4,8,3] <= 0
 -infinity <= delta[4,8,4] <= 0
 -infinity <= delta[4,8,5] <= 0
 -infinity <= delta[4,8,6] <= 0
 -infinity <= delta[4,8,7] <= 0
 -infinity <= delta[4,8,8] <= 0
 -infinity <= delta[4,8,9] <= 0
 -infinity <= delta[4,9,0] <= 0
 -infinity <= delta[4,9,1] <= 0
 -infinity <= delta[4,9,2] <= 0
 -infinity <= delta[4,9,3] <= 0
 -infinity <= delta[4,9,4] <= 0
 -infinity <= delta[4,9,5] <= 0
 -infinity <= delta[4,9,6] <= 0
 -infinity <= delta[4,9,7] <= 0
 -infinity <= delta[4,9,8] <= 0
 -infinity <= delta[4,9,9] <= 0
 -infinity <= lambda[0,0] <= 0
 -infinity <= lambda[0,1] <= 0
 -infinity <= lambda[0,2] <= 0
 -infinity <= lambda[0,3] <= 0
 -infinity <= lambda[0,4] <= 0
 -infinity <= lambda[0,5] <= 0
 -infinity <= lambda[0,6] <= 0
 -infinity <= lambda[0,7] <= 0
 -infinity <= lambda[0,8] <= 0
 -infinity <= lambda[0,9] <= 0
 -infinity <= lambda[1,0] <= 0
 -infinity <= lambda[1,1] <= 0
 -infinity <= lambda[1,2] <= 0
 -infinity <= lambda[1,3] <= 0
 -infinity <= lambda[1,4] <= 0
 -infinity <= lambda[1,5] <= 0
 -infinity <= lambda[1,6] <= 0
 -infinity <= lambda[1,7] <= 0
 -infinity <= lambda[1,8] <= 0
 -infinity <= lambda[1,9] <= 0
 -infinity <= lambda[2,0] <= 0
 -infinity <= lambda[2,1] <= 0
 -infinity <= lambda[2,2] <= 0
 -infinity <= lambda[2,3] <= 0
 -infinity <= lambda[2,4] <= 0
 -infinity <= lambda[2,5] <= 0
 -infinity <= lambda[2,6] <= 0
 -infinity <= lambda[2,7] <= 0
 -infinity <= lambda[2,8] <= 0
 -infinity <= lambda[2,9] <= 0
 -infinity <= lambda[3,0] <= 0
 -infinity <= lambda[3,1] <= 0
 -infinity <= lambda[3,2] <= 0
 -infinity <= lambda[3,3] <= 0
 -infinity <= lambda[3,4] <= 0
 -infinity <= lambda[3,5] <= 0
 -infinity <= lambda[3,6] <= 0
 -infinity <= lambda[3,7] <= 0
 -infinity <= lambda[3,8] <= 0
 -infinity <= lambda[3,9] <= 0
 -infinity <= lambda[4,0] <= 0
 -infinity <= lambda[4,1] <= 0
 -infinity <= lambda[4,2] <= 0
 -infinity <= lambda[4,3] <= 0
 -infinity <= lambda[4,4] <= 0
 -infinity <= lambda[4,5] <= 0
 -infinity <= lambda[4,6] <= 0
 -infinity <= lambda[4,7] <= 0
 -infinity <= lambda[4,8] <= 0
 -infinity <= lambda[4,9] <= 0
 -infinity <= mu[0,0] <= 0
 -infinity <= mu[0,1] <= 0
 -infinity <= mu[0,2] <= 0
 -infinity <= mu[0,3] <= 0
 -infinity <= mu[0,4] <= 0
 -infinity <= mu[0,5] <= 0
 -infinity <= mu[0,6] <= 0
 -infinity <= mu[0,7] <= 0
 -infinity <= mu[0,8] <= 0
 -infinity <= mu[0,9] <= 0
 -infinity <= mu[1,0] <= 0
 -infinity <= mu[1,1] <= 0
 -infinity <= mu[1,2] <= 0
 -infinity <= mu[1,3] <= 0
 -infinity <= mu[1,4] <= 0
 -infinity <= mu[1,5] <= 0
 -infinity <= mu[1,6] <= 0
 -infinity <= mu[1,7] <= 0
 -infinity <= mu[1,8] <= 0
 -infinity <= mu[1,9] <= 0
 -infinity <= mu[2,0] <= 0
 -infinity <= mu[2,1] <= 0
 -infinity <= mu[2,2] <= 0
 -infinity <= mu[2,3] <= 0
 -infinity <= mu[2,4] <= 0
 -infinity <= mu[2,5] <= 0
 -infinity <= mu[2,6] <= 0
 -infinity <= mu[2,7] <= 0
 -infinity <= mu[2,8] <= 0
 -infinity <= mu[2,9] <= 0
 -infinity <= mu[3,0] <= 0
 -infinity <= mu[3,1] <= 0
 -infinity <= mu[3,2] <= 0
 -infinity <= mu[3,3] <= 0
 -infinity <= mu[3,4] <= 0
 -infinity <= mu[3,5] <= 0
 -infinity <= mu[3,6] <= 0
 -infinity <= mu[3,7] <= 0
 -infinity <= mu[3,8] <= 0
 -infinity <= mu[3,9] <= 0
 -infinity <= mu[4,0] <= 0
 -infinity <= mu[4,1] <= 0
 -infinity <= mu[4,2] <= 0
 -infinity <= mu[4,3] <= 0
 -infinity <= mu[4,4] <= 0
 -infinity <= mu[4,5] <= 0
 -infinity <= mu[4,6] <= 0
 -infinity <= mu[4,7] <= 0
 -infinity <= mu[4,8] <= 0
 -infinity <= mu[4,9] <= 0
End
