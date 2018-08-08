\ Model sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0] + Qk[1] + Qk[2] + Qk[3] + Qk[4] + Qk[5] + Qk[6] + Qk[7] + Qk[8]
   + Qk[9]
Subject To
 Q(k)[0]: - gamma[0,0,0] - gamma[0,0,1] - gamma[0,0,7] - gamma[0,1,0]
   - gamma[0,1,1] - gamma[0,1,7] - gamma[0,2,0] - gamma[0,2,1]
   - gamma[0,2,7] - gamma[0,3,0] - gamma[0,3,1] - gamma[0,3,7]
   - gamma[0,4,0] - gamma[0,4,1] - gamma[0,4,7] - gamma[0,5,0]
   - gamma[0,5,1] - gamma[0,5,7] - gamma[0,6,0] - gamma[0,6,1]
   - gamma[0,6,7] - gamma[0,7,0] - gamma[0,7,1] - gamma[0,7,7]
   - gamma[0,8,0] - gamma[0,8,1] - gamma[0,8,7] - gamma[0,9,0]
   - gamma[0,9,1] - gamma[0,9,7] - delta[0,0,0] - delta[0,0,3]
   - delta[0,0,4] - delta[0,0,5] - delta[0,0,6] - delta[0,0,7]
   - delta[0,0,8] - delta[0,0,9] - delta[0,1,0] - delta[0,1,3]
   - delta[0,1,4] - delta[0,1,5] - delta[0,1,6] - delta[0,1,7]
   - delta[0,1,8] - delta[0,1,9] - delta[0,2,0] - delta[0,2,3]
   - delta[0,2,4] - delta[0,2,5] - delta[0,2,6] - delta[0,2,7]
   - delta[0,2,8] - delta[0,2,9] - delta[0,3,0] - delta[0,3,3]
   - delta[0,3,4] - delta[0,3,5] - delta[0,3,6] - delta[0,3,7]
   - delta[0,3,8] - delta[0,3,9] - delta[0,4,0] - delta[0,4,3]
   - delta[0,4,4] - delta[0,4,5] - delta[0,4,6] - delta[0,4,7]
   - delta[0,4,8] - delta[0,4,9] - delta[0,5,0] - delta[0,5,3]
   - delta[0,5,4] - delta[0,5,5] - delta[0,5,6] - delta[0,5,7]
   - delta[0,5,8] - delta[0,5,9] - delta[0,6,0] - delta[0,6,3]
   - delta[0,6,4] - delta[0,6,5] - delta[0,6,6] - delta[0,6,7]
   - delta[0,6,8] - delta[0,6,9] - delta[0,7,0] - delta[0,7,3]
   - delta[0,7,4] - delta[0,7,5] - delta[0,7,6] - delta[0,7,7]
   - delta[0,7,8] - delta[0,7,9] - delta[0,8,0] - delta[0,8,3]
   - delta[0,8,4] - delta[0,8,5] - delta[0,8,6] - delta[0,8,7]
   - delta[0,8,8] - delta[0,8,9] - delta[0,9,0] - delta[0,9,3]
   - delta[0,9,4] - delta[0,9,5] - delta[0,9,6] - delta[0,9,7]
   - delta[0,9,8] - delta[0,9,9] - epsilon[0,0] - epsilon[0,1]
   - epsilon[0,2] - epsilon[0,3] - epsilon[0,4] - epsilon[0,5]
   - epsilon[0,6] - epsilon[0,7] - epsilon[0,8] - epsilon[0,9]
   - lambda[0,2] - lambda[0,3] - lambda[0,4] - lambda[0,5] - lambda[0,6]
   - lambda[0,8] - lambda[0,9] - mu[0,0] - mu[0,3] - mu[0,4] - mu[0,5]
   - mu[0,6] - mu[0,7] - mu[0,8] - mu[0,9] - nu[0] + Qk[0] = 0
 Q(k)[1]: - gamma[1,0,0] - gamma[1,0,1] - gamma[1,0,7] - gamma[1,1,0]
   - gamma[1,1,1] - gamma[1,1,7] - gamma[1,2,0] - gamma[1,2,1]
   - gamma[1,2,7] - gamma[1,3,0] - gamma[1,3,1] - gamma[1,3,7]
   - gamma[1,4,0] - gamma[1,4,1] - gamma[1,4,7] - gamma[1,5,0]
   - gamma[1,5,1] - gamma[1,5,7] - gamma[1,6,0] - gamma[1,6,1]
   - gamma[1,6,7] - gamma[1,7,0] - gamma[1,7,1] - gamma[1,7,7]
   - gamma[1,8,0] - gamma[1,8,1] - gamma[1,8,7] - gamma[1,9,0]
   - gamma[1,9,1] - gamma[1,9,7] - delta[1,0,0] - delta[1,0,2]
   - delta[1,0,3] - delta[1,0,5] - delta[1,0,6] - delta[1,0,7]
   - delta[1,0,8] - delta[1,0,9] - delta[1,1,0] - delta[1,1,2]
   - delta[1,1,3] - delta[1,1,5] - delta[1,1,6] - delta[1,1,7]
   - delta[1,1,8] - delta[1,1,9] - delta[1,2,0] - delta[1,2,2]
   - delta[1,2,3] - delta[1,2,5] - delta[1,2,6] - delta[1,2,7]
   - delta[1,2,8] - delta[1,2,9] - delta[1,3,0] - delta[1,3,2]
   - delta[1,3,3] - delta[1,3,5] - delta[1,3,6] - delta[1,3,7]
   - delta[1,3,8] - delta[1,3,9] - delta[1,4,0] - delta[1,4,2]
   - delta[1,4,3] - delta[1,4,5] - delta[1,4,6] - delta[1,4,7]
   - delta[1,4,8] - delta[1,4,9] - delta[1,5,0] - delta[1,5,2]
   - delta[1,5,3] - delta[1,5,5] - delta[1,5,6] - delta[1,5,7]
   - delta[1,5,8] - delta[1,5,9] - delta[1,6,0] - delta[1,6,2]
   - delta[1,6,3] - delta[1,6,5] - delta[1,6,6] - delta[1,6,7]
   - delta[1,6,8] - delta[1,6,9] - delta[1,7,0] - delta[1,7,2]
   - delta[1,7,3] - delta[1,7,5] - delta[1,7,6] - delta[1,7,7]
   - delta[1,7,8] - delta[1,7,9] - delta[1,8,0] - delta[1,8,2]
   - delta[1,8,3] - delta[1,8,5] - delta[1,8,6] - delta[1,8,7]
   - delta[1,8,8] - delta[1,8,9] - delta[1,9,0] - delta[1,9,2]
   - delta[1,9,3] - delta[1,9,5] - delta[1,9,6] - delta[1,9,7]
   - delta[1,9,8] - delta[1,9,9] - epsilon[1,0] - epsilon[1,1]
   - epsilon[1,2] - epsilon[1,3] - epsilon[1,4] - epsilon[1,5]
   - epsilon[1,6] - epsilon[1,7] - epsilon[1,8] - epsilon[1,9]
   - lambda[1,2] - lambda[1,3] - lambda[1,4] - lambda[1,5] - lambda[1,6]
   - lambda[1,8] - lambda[1,9] - mu[1,0] - mu[1,2] - mu[1,3] - mu[1,5]
   - mu[1,6] - mu[1,7] - mu[1,8] - mu[1,9] - nu[1] + Qk[1] = 0
 Q(k)[2]: - gamma[2,0,0] - gamma[2,0,1] - gamma[2,0,7] - gamma[2,1,0]
   - gamma[2,1,1] - gamma[2,1,7] - gamma[2,2,0] - gamma[2,2,1]
   - gamma[2,2,7] - gamma[2,3,0] - gamma[2,3,1] - gamma[2,3,7]
   - gamma[2,4,0] - gamma[2,4,1] - gamma[2,4,7] - gamma[2,5,0]
   - gamma[2,5,1] - gamma[2,5,7] - gamma[2,6,0] - gamma[2,6,1]
   - gamma[2,6,7] - gamma[2,7,0] - gamma[2,7,1] - gamma[2,7,7]
   - gamma[2,8,0] - gamma[2,8,1] - gamma[2,8,7] - gamma[2,9,0]
   - gamma[2,9,1] - gamma[2,9,7] - delta[2,0,0] - delta[2,0,1]
   - delta[2,0,2] - delta[2,0,3] - delta[2,0,4] - delta[2,0,5]
   - delta[2,0,6] - delta[2,0,8] - delta[2,1,0] - delta[2,1,1]
   - delta[2,1,2] - delta[2,1,3] - delta[2,1,4] - delta[2,1,5]
   - delta[2,1,6] - delta[2,1,8] - delta[2,2,0] - delta[2,2,1]
   - delta[2,2,2] - delta[2,2,3] - delta[2,2,4] - delta[2,2,5]
   - delta[2,2,6] - delta[2,2,8] - delta[2,3,0] - delta[2,3,1]
   - delta[2,3,2] - delta[2,3,3] - delta[2,3,4] - delta[2,3,5]
   - delta[2,3,6] - delta[2,3,8] - delta[2,4,0] - delta[2,4,1]
   - delta[2,4,2] - delta[2,4,3] - delta[2,4,4] - delta[2,4,5]
   - delta[2,4,6] - delta[2,4,8] - delta[2,5,0] - delta[2,5,1]
   - delta[2,5,2] - delta[2,5,3] - delta[2,5,4] - delta[2,5,5]
   - delta[2,5,6] - delta[2,5,8] - delta[2,6,0] - delta[2,6,1]
   - delta[2,6,2] - delta[2,6,3] - delta[2,6,4] - delta[2,6,5]
   - delta[2,6,6] - delta[2,6,8] - delta[2,7,0] - delta[2,7,1]
   - delta[2,7,2] - delta[2,7,3] - delta[2,7,4] - delta[2,7,5]
   - delta[2,7,6] - delta[2,7,8] - delta[2,8,0] - delta[2,8,1]
   - delta[2,8,2] - delta[2,8,3] - delta[2,8,4] - delta[2,8,5]
   - delta[2,8,6] - delta[2,8,8] - delta[2,9,0] - delta[2,9,1]
   - delta[2,9,2] - delta[2,9,3] - delta[2,9,4] - delta[2,9,5]
   - delta[2,9,6] - delta[2,9,8] - epsilon[2,0] - epsilon[2,1]
   - epsilon[2,2] - epsilon[2,3] - epsilon[2,4] - epsilon[2,5]
   - epsilon[2,6] - epsilon[2,7] - epsilon[2,8] - epsilon[2,9]
   - lambda[2,2] - lambda[2,3] - lambda[2,4] - lambda[2,5] - lambda[2,6]
   - lambda[2,8] - lambda[2,9] - mu[2,0] - mu[2,1] - mu[2,2] - mu[2,3]
   - mu[2,4] - mu[2,5] - mu[2,6] - mu[2,8] - nu[2] + Qk[2] = 0
 Q(k)[3]: - gamma[3,0,0] - gamma[3,0,1] - gamma[3,0,7] - gamma[3,1,0]
   - gamma[3,1,1] - gamma[3,1,7] - gamma[3,2,0] - gamma[3,2,1]
   - gamma[3,2,7] - gamma[3,3,0] - gamma[3,3,1] - gamma[3,3,7]
   - gamma[3,4,0] - gamma[3,4,1] - gamma[3,4,7] - gamma[3,5,0]
   - gamma[3,5,1] - gamma[3,5,7] - gamma[3,6,0] - gamma[3,6,1]
   - gamma[3,6,7] - gamma[3,7,0] - gamma[3,7,1] - gamma[3,7,7]
   - gamma[3,8,0] - gamma[3,8,1] - gamma[3,8,7] - gamma[3,9,0]
   - gamma[3,9,1] - gamma[3,9,7] - delta[3,0,0] - delta[3,0,1]
   - delta[3,0,2] - delta[3,0,3] - delta[3,0,4] - delta[3,0,5]
   - delta[3,0,8] - delta[3,0,9] - delta[3,1,0] - delta[3,1,1]
   - delta[3,1,2] - delta[3,1,3] - delta[3,1,4] - delta[3,1,5]
   - delta[3,1,8] - delta[3,1,9] - delta[3,2,0] - delta[3,2,1]
   - delta[3,2,2] - delta[3,2,3] - delta[3,2,4] - delta[3,2,5]
   - delta[3,2,8] - delta[3,2,9] - delta[3,3,0] - delta[3,3,1]
   - delta[3,3,2] - delta[3,3,3] - delta[3,3,4] - delta[3,3,5]
   - delta[3,3,8] - delta[3,3,9] - delta[3,4,0] - delta[3,4,1]
   - delta[3,4,2] - delta[3,4,3] - delta[3,4,4] - delta[3,4,5]
   - delta[3,4,8] - delta[3,4,9] - delta[3,5,0] - delta[3,5,1]
   - delta[3,5,2] - delta[3,5,3] - delta[3,5,4] - delta[3,5,5]
   - delta[3,5,8] - delta[3,5,9] - delta[3,6,0] - delta[3,6,1]
   - delta[3,6,2] - delta[3,6,3] - delta[3,6,4] - delta[3,6,5]
   - delta[3,6,8] - delta[3,6,9] - delta[3,7,0] - delta[3,7,1]
   - delta[3,7,2] - delta[3,7,3] - delta[3,7,4] - delta[3,7,5]
   - delta[3,7,8] - delta[3,7,9] - delta[3,8,0] - delta[3,8,1]
   - delta[3,8,2] - delta[3,8,3] - delta[3,8,4] - delta[3,8,5]
   - delta[3,8,8] - delta[3,8,9] - delta[3,9,0] - delta[3,9,1]
   - delta[3,9,2] - delta[3,9,3] - delta[3,9,4] - delta[3,9,5]
   - delta[3,9,8] - delta[3,9,9] - epsilon[3,0] - epsilon[3,1]
   - epsilon[3,2] - epsilon[3,3] - epsilon[3,4] - epsilon[3,5]
   - epsilon[3,6] - epsilon[3,7] - epsilon[3,8] - epsilon[3,9]
   - lambda[3,2] - lambda[3,3] - lambda[3,4] - lambda[3,5] - lambda[3,6]
   - lambda[3,8] - lambda[3,9] - mu[3,0] - mu[3,1] - mu[3,2] - mu[3,3]
   - mu[3,4] - mu[3,5] - mu[3,8] - mu[3,9] - nu[3] + Qk[3] = 0
 Q(k)[4]: - gamma[4,0,0] - gamma[4,0,1] - gamma[4,0,7] - gamma[4,1,0]
   - gamma[4,1,1] - gamma[4,1,7] - gamma[4,2,0] - gamma[4,2,1]
   - gamma[4,2,7] - gamma[4,3,0] - gamma[4,3,1] - gamma[4,3,7]
   - gamma[4,4,0] - gamma[4,4,1] - gamma[4,4,7] - gamma[4,5,0]
   - gamma[4,5,1] - gamma[4,5,7] - gamma[4,6,0] - gamma[4,6,1]
   - gamma[4,6,7] - gamma[4,7,0] - gamma[4,7,1] - gamma[4,7,7]
   - gamma[4,8,0] - gamma[4,8,1] - gamma[4,8,7] - gamma[4,9,0]
   - gamma[4,9,1] - gamma[4,9,7] - delta[4,0,0] - delta[4,0,2]
   - delta[4,0,4] - delta[4,0,5] - delta[4,0,6] - delta[4,0,7]
   - delta[4,0,8] - delta[4,0,9] - delta[4,1,0] - delta[4,1,2]
   - delta[4,1,4] - delta[4,1,5] - delta[4,1,6] - delta[4,1,7]
   - delta[4,1,8] - delta[4,1,9] - delta[4,2,0] - delta[4,2,2]
   - delta[4,2,4] - delta[4,2,5] - delta[4,2,6] - delta[4,2,7]
   - delta[4,2,8] - delta[4,2,9] - delta[4,3,0] - delta[4,3,2]
   - delta[4,3,4] - delta[4,3,5] - delta[4,3,6] - delta[4,3,7]
   - delta[4,3,8] - delta[4,3,9] - delta[4,4,0] - delta[4,4,2]
   - delta[4,4,4] - delta[4,4,5] - delta[4,4,6] - delta[4,4,7]
   - delta[4,4,8] - delta[4,4,9] - delta[4,5,0] - delta[4,5,2]
   - delta[4,5,4] - delta[4,5,5] - delta[4,5,6] - delta[4,5,7]
   - delta[4,5,8] - delta[4,5,9] - delta[4,6,0] - delta[4,6,2]
   - delta[4,6,4] - delta[4,6,5] - delta[4,6,6] - delta[4,6,7]
   - delta[4,6,8] - delta[4,6,9] - delta[4,7,0] - delta[4,7,2]
   - delta[4,7,4] - delta[4,7,5] - delta[4,7,6] - delta[4,7,7]
   - delta[4,7,8] - delta[4,7,9] - delta[4,8,0] - delta[4,8,2]
   - delta[4,8,4] - delta[4,8,5] - delta[4,8,6] - delta[4,8,7]
   - delta[4,8,8] - delta[4,8,9] - delta[4,9,0] - delta[4,9,2]
   - delta[4,9,4] - delta[4,9,5] - delta[4,9,6] - delta[4,9,7]
   - delta[4,9,8] - delta[4,9,9] - epsilon[4,0] - epsilon[4,1]
   - epsilon[4,2] - epsilon[4,3] - epsilon[4,4] - epsilon[4,5]
   - epsilon[4,6] - epsilon[4,7] - epsilon[4,8] - epsilon[4,9]
   - lambda[4,2] - lambda[4,3] - lambda[4,4] - lambda[4,5] - lambda[4,6]
   - lambda[4,8] - lambda[4,9] - mu[4,0] - mu[4,2] - mu[4,4] - mu[4,5]
   - mu[4,6] - mu[4,7] - mu[4,8] - mu[4,9] - nu[4] + Qk[4] = 0
 Q(k)[5]: - gamma[5,0,0] - gamma[5,0,1] - gamma[5,0,7] - gamma[5,1,0]
   - gamma[5,1,1] - gamma[5,1,7] - gamma[5,2,0] - gamma[5,2,1]
   - gamma[5,2,7] - gamma[5,3,0] - gamma[5,3,1] - gamma[5,3,7]
   - gamma[5,4,0] - gamma[5,4,1] - gamma[5,4,7] - gamma[5,5,0]
   - gamma[5,5,1] - gamma[5,5,7] - gamma[5,6,0] - gamma[5,6,1]
   - gamma[5,6,7] - gamma[5,7,0] - gamma[5,7,1] - gamma[5,7,7]
   - gamma[5,8,0] - gamma[5,8,1] - gamma[5,8,7] - gamma[5,9,0]
   - gamma[5,9,1] - gamma[5,9,7] - delta[5,0,1] - delta[5,0,2]
   - delta[5,0,3] - delta[5,0,4] - delta[5,0,5] - delta[5,0,6]
   - delta[5,0,8] - delta[5,0,9] - delta[5,1,1] - delta[5,1,2]
   - delta[5,1,3] - delta[5,1,4] - delta[5,1,5] - delta[5,1,6]
   - delta[5,1,8] - delta[5,1,9] - delta[5,2,1] - delta[5,2,2]
   - delta[5,2,3] - delta[5,2,4] - delta[5,2,5] - delta[5,2,6]
   - delta[5,2,8] - delta[5,2,9] - delta[5,3,1] - delta[5,3,2]
   - delta[5,3,3] - delta[5,3,4] - delta[5,3,5] - delta[5,3,6]
   - delta[5,3,8] - delta[5,3,9] - delta[5,4,1] - delta[5,4,2]
   - delta[5,4,3] - delta[5,4,4] - delta[5,4,5] - delta[5,4,6]
   - delta[5,4,8] - delta[5,4,9] - delta[5,5,1] - delta[5,5,2]
   - delta[5,5,3] - delta[5,5,4] - delta[5,5,5] - delta[5,5,6]
   - delta[5,5,8] - delta[5,5,9] - delta[5,6,1] - delta[5,6,2]
   - delta[5,6,3] - delta[5,6,4] - delta[5,6,5] - delta[5,6,6]
   - delta[5,6,8] - delta[5,6,9] - delta[5,7,1] - delta[5,7,2]
   - delta[5,7,3] - delta[5,7,4] - delta[5,7,5] - delta[5,7,6]
   - delta[5,7,8] - delta[5,7,9] - delta[5,8,1] - delta[5,8,2]
   - delta[5,8,3] - delta[5,8,4] - delta[5,8,5] - delta[5,8,6]
   - delta[5,8,8] - delta[5,8,9] - delta[5,9,1] - delta[5,9,2]
   - delta[5,9,3] - delta[5,9,4] - delta[5,9,5] - delta[5,9,6]
   - delta[5,9,8] - delta[5,9,9] - epsilon[5,0] - epsilon[5,1]
   - epsilon[5,2] - epsilon[5,3] - epsilon[5,4] - epsilon[5,5]
   - epsilon[5,6] - epsilon[5,7] - epsilon[5,8] - epsilon[5,9]
   - lambda[5,2] - lambda[5,3] - lambda[5,4] - lambda[5,5] - lambda[5,6]
   - lambda[5,8] - lambda[5,9] - mu[5,1] - mu[5,2] - mu[5,3] - mu[5,4]
   - mu[5,5] - mu[5,6] - mu[5,8] - mu[5,9] - 2 nu[5] + Qk[5] = 0
 Q(k)[6]: - gamma[6,0,0] - gamma[6,0,1] - gamma[6,0,7] - gamma[6,1,0]
   - gamma[6,1,1] - gamma[6,1,7] - gamma[6,2,0] - gamma[6,2,1]
   - gamma[6,2,7] - gamma[6,3,0] - gamma[6,3,1] - gamma[6,3,7]
   - gamma[6,4,0] - gamma[6,4,1] - gamma[6,4,7] - gamma[6,5,0]
   - gamma[6,5,1] - gamma[6,5,7] - gamma[6,6,0] - gamma[6,6,1]
   - gamma[6,6,7] - gamma[6,7,0] - gamma[6,7,1] - gamma[6,7,7]
   - gamma[6,8,0] - gamma[6,8,1] - gamma[6,8,7] - gamma[6,9,0]
   - gamma[6,9,1] - gamma[6,9,7] - delta[6,0,0] - delta[6,0,1]
   - delta[6,0,2] - delta[6,0,3] - delta[6,0,4] - delta[6,0,5]
   - delta[6,0,7] - delta[6,0,8] - delta[6,1,0] - delta[6,1,1]
   - delta[6,1,2] - delta[6,1,3] - delta[6,1,4] - delta[6,1,5]
   - delta[6,1,7] - delta[6,1,8] - delta[6,2,0] - delta[6,2,1]
   - delta[6,2,2] - delta[6,2,3] - delta[6,2,4] - delta[6,2,5]
   - delta[6,2,7] - delta[6,2,8] - delta[6,3,0] - delta[6,3,1]
   - delta[6,3,2] - delta[6,3,3] - delta[6,3,4] - delta[6,3,5]
   - delta[6,3,7] - delta[6,3,8] - delta[6,4,0] - delta[6,4,1]
   - delta[6,4,2] - delta[6,4,3] - delta[6,4,4] - delta[6,4,5]
   - delta[6,4,7] - delta[6,4,8] - delta[6,5,0] - delta[6,5,1]
   - delta[6,5,2] - delta[6,5,3] - delta[6,5,4] - delta[6,5,5]
   - delta[6,5,7] - delta[6,5,8] - delta[6,6,0] - delta[6,6,1]
   - delta[6,6,2] - delta[6,6,3] - delta[6,6,4] - delta[6,6,5]
   - delta[6,6,7] - delta[6,6,8] - delta[6,7,0] - delta[6,7,1]
   - delta[6,7,2] - delta[6,7,3] - delta[6,7,4] - delta[6,7,5]
   - delta[6,7,7] - delta[6,7,8] - delta[6,8,0] - delta[6,8,1]
   - delta[6,8,2] - delta[6,8,3] - delta[6,8,4] - delta[6,8,5]
   - delta[6,8,7] - delta[6,8,8] - delta[6,9,0] - delta[6,9,1]
   - delta[6,9,2] - delta[6,9,3] - delta[6,9,4] - delta[6,9,5]
   - delta[6,9,7] - delta[6,9,8] - epsilon[6,0] - epsilon[6,1]
   - epsilon[6,2] - epsilon[6,3] - epsilon[6,4] - epsilon[6,5]
   - epsilon[6,6] - epsilon[6,7] - epsilon[6,8] - epsilon[6,9]
   - lambda[6,2] - lambda[6,3] - lambda[6,4] - lambda[6,5] - lambda[6,6]
   - lambda[6,8] - lambda[6,9] - mu[6,0] - mu[6,1] - mu[6,2] - mu[6,3]
   - mu[6,4] - mu[6,5] - mu[6,7] - mu[6,8] + Qk[6] = 0
 Q(k)[7]: - gamma[7,0,0] - gamma[7,0,1] - gamma[7,0,7] - gamma[7,1,0]
   - gamma[7,1,1] - gamma[7,1,7] - gamma[7,2,0] - gamma[7,2,1]
   - gamma[7,2,7] - gamma[7,3,0] - gamma[7,3,1] - gamma[7,3,7]
   - gamma[7,4,0] - gamma[7,4,1] - gamma[7,4,7] - gamma[7,5,0]
   - gamma[7,5,1] - gamma[7,5,7] - gamma[7,6,0] - gamma[7,6,1]
   - gamma[7,6,7] - gamma[7,7,0] - gamma[7,7,1] - gamma[7,7,7]
   - gamma[7,8,0] - gamma[7,8,1] - gamma[7,8,7] - gamma[7,9,0]
   - gamma[7,9,1] - gamma[7,9,7] - delta[7,0,1] - delta[7,0,2]
   - delta[7,0,3] - delta[7,0,4] - delta[7,0,5] - delta[7,0,6]
   - delta[7,0,7] - delta[7,0,8] - delta[7,1,1] - delta[7,1,2]
   - delta[7,1,3] - delta[7,1,4] - delta[7,1,5] - delta[7,1,6]
   - delta[7,1,7] - delta[7,1,8] - delta[7,2,1] - delta[7,2,2]
   - delta[7,2,3] - delta[7,2,4] - delta[7,2,5] - delta[7,2,6]
   - delta[7,2,7] - delta[7,2,8] - delta[7,3,1] - delta[7,3,2]
   - delta[7,3,3] - delta[7,3,4] - delta[7,3,5] - delta[7,3,6]
   - delta[7,3,7] - delta[7,3,8] - delta[7,4,1] - delta[7,4,2]
   - delta[7,4,3] - delta[7,4,4] - delta[7,4,5] - delta[7,4,6]
   - delta[7,4,7] - delta[7,4,8] - delta[7,5,1] - delta[7,5,2]
   - delta[7,5,3] - delta[7,5,4] - delta[7,5,5] - delta[7,5,6]
   - delta[7,5,7] - delta[7,5,8] - delta[7,6,1] - delta[7,6,2]
   - delta[7,6,3] - delta[7,6,4] - delta[7,6,5] - delta[7,6,6]
   - delta[7,6,7] - delta[7,6,8] - delta[7,7,1] - delta[7,7,2]
   - delta[7,7,3] - delta[7,7,4] - delta[7,7,5] - delta[7,7,6]
   - delta[7,7,7] - delta[7,7,8] - delta[7,8,1] - delta[7,8,2]
   - delta[7,8,3] - delta[7,8,4] - delta[7,8,5] - delta[7,8,6]
   - delta[7,8,7] - delta[7,8,8] - delta[7,9,1] - delta[7,9,2]
   - delta[7,9,3] - delta[7,9,4] - delta[7,9,5] - delta[7,9,6]
   - delta[7,9,7] - delta[7,9,8] - epsilon[7,0] - epsilon[7,1]
   - epsilon[7,2] - epsilon[7,3] - epsilon[7,4] - epsilon[7,5]
   - epsilon[7,6] - epsilon[7,7] - epsilon[7,8] - epsilon[7,9]
   - lambda[7,2] - lambda[7,3] - lambda[7,4] - lambda[7,5] - lambda[7,6]
   - lambda[7,8] - lambda[7,9] - mu[7,1] - mu[7,2] - mu[7,3] - mu[7,4]
   - mu[7,5] - mu[7,6] - mu[7,7] - mu[7,8] - nu[7] + Qk[7] = 0
 Q(k)[8]: - gamma[8,0,0] - gamma[8,0,1] - gamma[8,0,7] - gamma[8,1,0]
   - gamma[8,1,1] - gamma[8,1,7] - gamma[8,2,0] - gamma[8,2,1]
   - gamma[8,2,7] - gamma[8,3,0] - gamma[8,3,1] - gamma[8,3,7]
   - gamma[8,4,0] - gamma[8,4,1] - gamma[8,4,7] - gamma[8,5,0]
   - gamma[8,5,1] - gamma[8,5,7] - gamma[8,6,0] - gamma[8,6,1]
   - gamma[8,6,7] - gamma[8,7,0] - gamma[8,7,1] - gamma[8,7,7]
   - gamma[8,8,0] - gamma[8,8,1] - gamma[8,8,7] - gamma[8,9,0]
   - gamma[8,9,1] - gamma[8,9,7] - delta[8,0,0] - delta[8,0,1]
   - delta[8,0,2] - delta[8,0,3] - delta[8,0,5] - delta[8,0,6]
   - delta[8,0,8] - delta[8,0,9] - delta[8,1,0] - delta[8,1,1]
   - delta[8,1,2] - delta[8,1,3] - delta[8,1,5] - delta[8,1,6]
   - delta[8,1,8] - delta[8,1,9] - delta[8,2,0] - delta[8,2,1]
   - delta[8,2,2] - delta[8,2,3] - delta[8,2,5] - delta[8,2,6]
   - delta[8,2,8] - delta[8,2,9] - delta[8,3,0] - delta[8,3,1]
   - delta[8,3,2] - delta[8,3,3] - delta[8,3,5] - delta[8,3,6]
   - delta[8,3,8] - delta[8,3,9] - delta[8,4,0] - delta[8,4,1]
   - delta[8,4,2] - delta[8,4,3] - delta[8,4,5] - delta[8,4,6]
   - delta[8,4,8] - delta[8,4,9] - delta[8,5,0] - delta[8,5,1]
   - delta[8,5,2] - delta[8,5,3] - delta[8,5,5] - delta[8,5,6]
   - delta[8,5,8] - delta[8,5,9] - delta[8,6,0] - delta[8,6,1]
   - delta[8,6,2] - delta[8,6,3] - delta[8,6,5] - delta[8,6,6]
   - delta[8,6,8] - delta[8,6,9] - delta[8,7,0] - delta[8,7,1]
   - delta[8,7,2] - delta[8,7,3] - delta[8,7,5] - delta[8,7,6]
   - delta[8,7,8] - delta[8,7,9] - delta[8,8,0] - delta[8,8,1]
   - delta[8,8,2] - delta[8,8,3] - delta[8,8,5] - delta[8,8,6]
   - delta[8,8,8] - delta[8,8,9] - delta[8,9,0] - delta[8,9,1]
   - delta[8,9,2] - delta[8,9,3] - delta[8,9,5] - delta[8,9,6]
   - delta[8,9,8] - delta[8,9,9] - epsilon[8,0] - epsilon[8,1]
   - epsilon[8,2] - epsilon[8,3] - epsilon[8,4] - epsilon[8,5]
   - epsilon[8,6] - epsilon[8,7] - epsilon[8,8] - epsilon[8,9]
   - lambda[8,2] - lambda[8,3] - lambda[8,4] - lambda[8,5] - lambda[8,6]
   - lambda[8,8] - lambda[8,9] - mu[8,0] - mu[8,1] - mu[8,2] - mu[8,3]
   - mu[8,5] - mu[8,6] - mu[8,8] - mu[8,9] - nu[8] + Qk[8] = 0
 Q(k)[9]: - gamma[9,0,0] - gamma[9,0,1] - gamma[9,0,7] - gamma[9,1,0]
   - gamma[9,1,1] - gamma[9,1,7] - gamma[9,2,0] - gamma[9,2,1]
   - gamma[9,2,7] - gamma[9,3,0] - gamma[9,3,1] - gamma[9,3,7]
   - gamma[9,4,0] - gamma[9,4,1] - gamma[9,4,7] - gamma[9,5,0]
   - gamma[9,5,1] - gamma[9,5,7] - gamma[9,6,0] - gamma[9,6,1]
   - gamma[9,6,7] - gamma[9,7,0] - gamma[9,7,1] - gamma[9,7,7]
   - gamma[9,8,0] - gamma[9,8,1] - gamma[9,8,7] - gamma[9,9,0]
   - gamma[9,9,1] - gamma[9,9,7] - delta[9,0,0] - delta[9,0,2]
   - delta[9,0,4] - delta[9,0,5] - delta[9,0,6] - delta[9,0,7]
   - delta[9,0,8] - delta[9,0,9] - delta[9,1,0] - delta[9,1,2]
   - delta[9,1,4] - delta[9,1,5] - delta[9,1,6] - delta[9,1,7]
   - delta[9,1,8] - delta[9,1,9] - delta[9,2,0] - delta[9,2,2]
   - delta[9,2,4] - delta[9,2,5] - delta[9,2,6] - delta[9,2,7]
   - delta[9,2,8] - delta[9,2,9] - delta[9,3,0] - delta[9,3,2]
   - delta[9,3,4] - delta[9,3,5] - delta[9,3,6] - delta[9,3,7]
   - delta[9,3,8] - delta[9,3,9] - delta[9,4,0] - delta[9,4,2]
   - delta[9,4,4] - delta[9,4,5] - delta[9,4,6] - delta[9,4,7]
   - delta[9,4,8] - delta[9,4,9] - delta[9,5,0] - delta[9,5,2]
   - delta[9,5,4] - delta[9,5,5] - delta[9,5,6] - delta[9,5,7]
   - delta[9,5,8] - delta[9,5,9] - delta[9,6,0] - delta[9,6,2]
   - delta[9,6,4] - delta[9,6,5] - delta[9,6,6] - delta[9,6,7]
   - delta[9,6,8] - delta[9,6,9] - delta[9,7,0] - delta[9,7,2]
   - delta[9,7,4] - delta[9,7,5] - delta[9,7,6] - delta[9,7,7]
   - delta[9,7,8] - delta[9,7,9] - delta[9,8,0] - delta[9,8,2]
   - delta[9,8,4] - delta[9,8,5] - delta[9,8,6] - delta[9,8,7]
   - delta[9,8,8] - delta[9,8,9] - delta[9,9,0] - delta[9,9,2]
   - delta[9,9,4] - delta[9,9,5] - delta[9,9,6] - delta[9,9,7]
   - delta[9,9,8] - delta[9,9,9] - epsilon[9,0] - epsilon[9,1]
   - epsilon[9,2] - epsilon[9,3] - epsilon[9,4] - epsilon[9,5]
   - epsilon[9,6] - epsilon[9,7] - epsilon[9,8] - epsilon[9,9]
   - lambda[9,2] - lambda[9,3] - lambda[9,4] - lambda[9,5] - lambda[9,6]
   - lambda[9,8] - lambda[9,9] - mu[9,0] - mu[9,2] - mu[9,4] - mu[9,5]
   - mu[9,6] - mu[9,7] - mu[9,8] - mu[9,9] - nu[9] + Qk[9] = 0
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
 u[5,0]: - gamma[5,0,0] - gamma[5,1,0] - gamma[5,2,0] - gamma[5,3,0]
   - gamma[5,4,0] - gamma[5,5,0] - gamma[5,6,0] - gamma[5,7,0]
   - gamma[5,8,0] - gamma[5,9,0] + lambda[5,0] + mu[5,0] + nu[5] <= 0
 u[5,1]: - gamma[5,0,1] - gamma[5,1,1] - gamma[5,2,1] - gamma[5,3,1]
   - gamma[5,4,1] - gamma[5,5,1] - gamma[5,6,1] - gamma[5,7,1]
   - gamma[5,8,1] - gamma[5,9,1] + lambda[5,1] + mu[5,1] + nu[5] <= 0
 u[5,2]: - gamma[5,0,2] - gamma[5,1,2] - gamma[5,2,2] - gamma[5,3,2]
   - gamma[5,4,2] - gamma[5,5,2] - gamma[5,6,2] - gamma[5,7,2]
   - gamma[5,8,2] - gamma[5,9,2] + lambda[5,2] + mu[5,2] + nu[5] <= 0
 u[5,3]: - gamma[5,0,3] - gamma[5,1,3] - gamma[5,2,3] - gamma[5,3,3]
   - gamma[5,4,3] - gamma[5,5,3] - gamma[5,6,3] - gamma[5,7,3]
   - gamma[5,8,3] - gamma[5,9,3] + lambda[5,3] + mu[5,3] + nu[5] <= 0
 u[5,4]: - gamma[5,0,4] - gamma[5,1,4] - gamma[5,2,4] - gamma[5,3,4]
   - gamma[5,4,4] - gamma[5,5,4] - gamma[5,6,4] - gamma[5,7,4]
   - gamma[5,8,4] - gamma[5,9,4] + lambda[5,4] + mu[5,4] + nu[5] <= 0
 u[5,5]: - gamma[5,0,5] - gamma[5,1,5] - gamma[5,2,5] - gamma[5,3,5]
   - gamma[5,4,5] - gamma[5,5,5] - gamma[5,6,5] - gamma[5,7,5]
   - gamma[5,8,5] - gamma[5,9,5] + lambda[5,5] + mu[5,5] + nu[5] <= 0
 u[5,6]: - gamma[5,0,6] - gamma[5,1,6] - gamma[5,2,6] - gamma[5,3,6]
   - gamma[5,4,6] - gamma[5,5,6] - gamma[5,6,6] - gamma[5,7,6]
   - gamma[5,8,6] - gamma[5,9,6] + lambda[5,6] + mu[5,6] + nu[5] <= 0
 u[5,7]: - gamma[5,0,7] - gamma[5,1,7] - gamma[5,2,7] - gamma[5,3,7]
   - gamma[5,4,7] - gamma[5,5,7] - gamma[5,6,7] - gamma[5,7,7]
   - gamma[5,8,7] - gamma[5,9,7] + lambda[5,7] + mu[5,7] + nu[5] <= 0
 u[5,8]: - gamma[5,0,8] - gamma[5,1,8] - gamma[5,2,8] - gamma[5,3,8]
   - gamma[5,4,8] - gamma[5,5,8] - gamma[5,6,8] - gamma[5,7,8]
   - gamma[5,8,8] - gamma[5,9,8] + lambda[5,8] + mu[5,8] + nu[5] <= 0
 u[5,9]: - gamma[5,0,9] - gamma[5,1,9] - gamma[5,2,9] - gamma[5,3,9]
   - gamma[5,4,9] - gamma[5,5,9] - gamma[5,6,9] - gamma[5,7,9]
   - gamma[5,8,9] - gamma[5,9,9] + lambda[5,9] + mu[5,9] + nu[5] <= 0
 u[6,0]: - gamma[6,0,0] - gamma[6,1,0] - gamma[6,2,0] - gamma[6,3,0]
   - gamma[6,4,0] - gamma[6,5,0] - gamma[6,6,0] - gamma[6,7,0]
   - gamma[6,8,0] - gamma[6,9,0] + lambda[6,0] + mu[6,0] + nu[6] <= 0
 u[6,1]: - gamma[6,0,1] - gamma[6,1,1] - gamma[6,2,1] - gamma[6,3,1]
   - gamma[6,4,1] - gamma[6,5,1] - gamma[6,6,1] - gamma[6,7,1]
   - gamma[6,8,1] - gamma[6,9,1] + lambda[6,1] + mu[6,1] + nu[6] <= 0
 u[6,2]: - gamma[6,0,2] - gamma[6,1,2] - gamma[6,2,2] - gamma[6,3,2]
   - gamma[6,4,2] - gamma[6,5,2] - gamma[6,6,2] - gamma[6,7,2]
   - gamma[6,8,2] - gamma[6,9,2] + lambda[6,2] + mu[6,2] + nu[6] <= 0
 u[6,3]: - gamma[6,0,3] - gamma[6,1,3] - gamma[6,2,3] - gamma[6,3,3]
   - gamma[6,4,3] - gamma[6,5,3] - gamma[6,6,3] - gamma[6,7,3]
   - gamma[6,8,3] - gamma[6,9,3] + lambda[6,3] + mu[6,3] + nu[6] <= 0
 u[6,4]: - gamma[6,0,4] - gamma[6,1,4] - gamma[6,2,4] - gamma[6,3,4]
   - gamma[6,4,4] - gamma[6,5,4] - gamma[6,6,4] - gamma[6,7,4]
   - gamma[6,8,4] - gamma[6,9,4] + lambda[6,4] + mu[6,4] + nu[6] <= 0
 u[6,5]: - gamma[6,0,5] - gamma[6,1,5] - gamma[6,2,5] - gamma[6,3,5]
   - gamma[6,4,5] - gamma[6,5,5] - gamma[6,6,5] - gamma[6,7,5]
   - gamma[6,8,5] - gamma[6,9,5] + lambda[6,5] + mu[6,5] + nu[6] <= 0
 u[6,6]: - gamma[6,0,6] - gamma[6,1,6] - gamma[6,2,6] - gamma[6,3,6]
   - gamma[6,4,6] - gamma[6,5,6] - gamma[6,6,6] - gamma[6,7,6]
   - gamma[6,8,6] - gamma[6,9,6] + lambda[6,6] + mu[6,6] + nu[6] <= 0
 u[6,7]: - gamma[6,0,7] - gamma[6,1,7] - gamma[6,2,7] - gamma[6,3,7]
   - gamma[6,4,7] - gamma[6,5,7] - gamma[6,6,7] - gamma[6,7,7]
   - gamma[6,8,7] - gamma[6,9,7] + lambda[6,7] + mu[6,7] + nu[6] <= 0
 u[6,8]: - gamma[6,0,8] - gamma[6,1,8] - gamma[6,2,8] - gamma[6,3,8]
   - gamma[6,4,8] - gamma[6,5,8] - gamma[6,6,8] - gamma[6,7,8]
   - gamma[6,8,8] - gamma[6,9,8] + lambda[6,8] + mu[6,8] + nu[6] <= 0
 u[6,9]: - gamma[6,0,9] - gamma[6,1,9] - gamma[6,2,9] - gamma[6,3,9]
   - gamma[6,4,9] - gamma[6,5,9] - gamma[6,6,9] - gamma[6,7,9]
   - gamma[6,8,9] - gamma[6,9,9] + lambda[6,9] + mu[6,9] + nu[6] <= 0
 u[7,0]: - gamma[7,0,0] - gamma[7,1,0] - gamma[7,2,0] - gamma[7,3,0]
   - gamma[7,4,0] - gamma[7,5,0] - gamma[7,6,0] - gamma[7,7,0]
   - gamma[7,8,0] - gamma[7,9,0] + lambda[7,0] + mu[7,0] + nu[7] <= 0
 u[7,1]: - gamma[7,0,1] - gamma[7,1,1] - gamma[7,2,1] - gamma[7,3,1]
   - gamma[7,4,1] - gamma[7,5,1] - gamma[7,6,1] - gamma[7,7,1]
   - gamma[7,8,1] - gamma[7,9,1] + lambda[7,1] + mu[7,1] + nu[7] <= 0
 u[7,2]: - gamma[7,0,2] - gamma[7,1,2] - gamma[7,2,2] - gamma[7,3,2]
   - gamma[7,4,2] - gamma[7,5,2] - gamma[7,6,2] - gamma[7,7,2]
   - gamma[7,8,2] - gamma[7,9,2] + lambda[7,2] + mu[7,2] + nu[7] <= 0
 u[7,3]: - gamma[7,0,3] - gamma[7,1,3] - gamma[7,2,3] - gamma[7,3,3]
   - gamma[7,4,3] - gamma[7,5,3] - gamma[7,6,3] - gamma[7,7,3]
   - gamma[7,8,3] - gamma[7,9,3] + lambda[7,3] + mu[7,3] + nu[7] <= 0
 u[7,4]: - gamma[7,0,4] - gamma[7,1,4] - gamma[7,2,4] - gamma[7,3,4]
   - gamma[7,4,4] - gamma[7,5,4] - gamma[7,6,4] - gamma[7,7,4]
   - gamma[7,8,4] - gamma[7,9,4] + lambda[7,4] + mu[7,4] + nu[7] <= 0
 u[7,5]: - gamma[7,0,5] - gamma[7,1,5] - gamma[7,2,5] - gamma[7,3,5]
   - gamma[7,4,5] - gamma[7,5,5] - gamma[7,6,5] - gamma[7,7,5]
   - gamma[7,8,5] - gamma[7,9,5] + lambda[7,5] + mu[7,5] + nu[7] <= 0
 u[7,6]: - gamma[7,0,6] - gamma[7,1,6] - gamma[7,2,6] - gamma[7,3,6]
   - gamma[7,4,6] - gamma[7,5,6] - gamma[7,6,6] - gamma[7,7,6]
   - gamma[7,8,6] - gamma[7,9,6] + lambda[7,6] + mu[7,6] + nu[7] <= 0
 u[7,7]: - gamma[7,0,7] - gamma[7,1,7] - gamma[7,2,7] - gamma[7,3,7]
   - gamma[7,4,7] - gamma[7,5,7] - gamma[7,6,7] - gamma[7,7,7]
   - gamma[7,8,7] - gamma[7,9,7] + lambda[7,7] + mu[7,7] + nu[7] <= 0
 u[7,8]: - gamma[7,0,8] - gamma[7,1,8] - gamma[7,2,8] - gamma[7,3,8]
   - gamma[7,4,8] - gamma[7,5,8] - gamma[7,6,8] - gamma[7,7,8]
   - gamma[7,8,8] - gamma[7,9,8] + lambda[7,8] + mu[7,8] + nu[7] <= 0
 u[7,9]: - gamma[7,0,9] - gamma[7,1,9] - gamma[7,2,9] - gamma[7,3,9]
   - gamma[7,4,9] - gamma[7,5,9] - gamma[7,6,9] - gamma[7,7,9]
   - gamma[7,8,9] - gamma[7,9,9] + lambda[7,9] + mu[7,9] + nu[7] <= 0
 u[8,0]: - gamma[8,0,0] - gamma[8,1,0] - gamma[8,2,0] - gamma[8,3,0]
   - gamma[8,4,0] - gamma[8,5,0] - gamma[8,6,0] - gamma[8,7,0]
   - gamma[8,8,0] - gamma[8,9,0] + lambda[8,0] + mu[8,0] + nu[8] <= 0
 u[8,1]: - gamma[8,0,1] - gamma[8,1,1] - gamma[8,2,1] - gamma[8,3,1]
   - gamma[8,4,1] - gamma[8,5,1] - gamma[8,6,1] - gamma[8,7,1]
   - gamma[8,8,1] - gamma[8,9,1] + lambda[8,1] + mu[8,1] + nu[8] <= 0
 u[8,2]: - gamma[8,0,2] - gamma[8,1,2] - gamma[8,2,2] - gamma[8,3,2]
   - gamma[8,4,2] - gamma[8,5,2] - gamma[8,6,2] - gamma[8,7,2]
   - gamma[8,8,2] - gamma[8,9,2] + lambda[8,2] + mu[8,2] + nu[8] <= 0
 u[8,3]: - gamma[8,0,3] - gamma[8,1,3] - gamma[8,2,3] - gamma[8,3,3]
   - gamma[8,4,3] - gamma[8,5,3] - gamma[8,6,3] - gamma[8,7,3]
   - gamma[8,8,3] - gamma[8,9,3] + lambda[8,3] + mu[8,3] + nu[8] <= 0
 u[8,4]: - gamma[8,0,4] - gamma[8,1,4] - gamma[8,2,4] - gamma[8,3,4]
   - gamma[8,4,4] - gamma[8,5,4] - gamma[8,6,4] - gamma[8,7,4]
   - gamma[8,8,4] - gamma[8,9,4] + lambda[8,4] + mu[8,4] + nu[8] <= 0
 u[8,5]: - gamma[8,0,5] - gamma[8,1,5] - gamma[8,2,5] - gamma[8,3,5]
   - gamma[8,4,5] - gamma[8,5,5] - gamma[8,6,5] - gamma[8,7,5]
   - gamma[8,8,5] - gamma[8,9,5] + lambda[8,5] + mu[8,5] + nu[8] <= 0
 u[8,6]: - gamma[8,0,6] - gamma[8,1,6] - gamma[8,2,6] - gamma[8,3,6]
   - gamma[8,4,6] - gamma[8,5,6] - gamma[8,6,6] - gamma[8,7,6]
   - gamma[8,8,6] - gamma[8,9,6] + lambda[8,6] + mu[8,6] + nu[8] <= 0
 u[8,7]: - gamma[8,0,7] - gamma[8,1,7] - gamma[8,2,7] - gamma[8,3,7]
   - gamma[8,4,7] - gamma[8,5,7] - gamma[8,6,7] - gamma[8,7,7]
   - gamma[8,8,7] - gamma[8,9,7] + lambda[8,7] + mu[8,7] + nu[8] <= 0
 u[8,8]: - gamma[8,0,8] - gamma[8,1,8] - gamma[8,2,8] - gamma[8,3,8]
   - gamma[8,4,8] - gamma[8,5,8] - gamma[8,6,8] - gamma[8,7,8]
   - gamma[8,8,8] - gamma[8,9,8] + lambda[8,8] + mu[8,8] + nu[8] <= 0
 u[8,9]: - gamma[8,0,9] - gamma[8,1,9] - gamma[8,2,9] - gamma[8,3,9]
   - gamma[8,4,9] - gamma[8,5,9] - gamma[8,6,9] - gamma[8,7,9]
   - gamma[8,8,9] - gamma[8,9,9] + lambda[8,9] + mu[8,9] + nu[8] <= 0
 u[9,0]: - gamma[9,0,0] - gamma[9,1,0] - gamma[9,2,0] - gamma[9,3,0]
   - gamma[9,4,0] - gamma[9,5,0] - gamma[9,6,0] - gamma[9,7,0]
   - gamma[9,8,0] - gamma[9,9,0] + lambda[9,0] + mu[9,0] + nu[9] <= 0
 u[9,1]: - gamma[9,0,1] - gamma[9,1,1] - gamma[9,2,1] - gamma[9,3,1]
   - gamma[9,4,1] - gamma[9,5,1] - gamma[9,6,1] - gamma[9,7,1]
   - gamma[9,8,1] - gamma[9,9,1] + lambda[9,1] + mu[9,1] + nu[9] <= 0
 u[9,2]: - gamma[9,0,2] - gamma[9,1,2] - gamma[9,2,2] - gamma[9,3,2]
   - gamma[9,4,2] - gamma[9,5,2] - gamma[9,6,2] - gamma[9,7,2]
   - gamma[9,8,2] - gamma[9,9,2] + lambda[9,2] + mu[9,2] + nu[9] <= 0
 u[9,3]: - gamma[9,0,3] - gamma[9,1,3] - gamma[9,2,3] - gamma[9,3,3]
   - gamma[9,4,3] - gamma[9,5,3] - gamma[9,6,3] - gamma[9,7,3]
   - gamma[9,8,3] - gamma[9,9,3] + lambda[9,3] + mu[9,3] + nu[9] <= 0
 u[9,4]: - gamma[9,0,4] - gamma[9,1,4] - gamma[9,2,4] - gamma[9,3,4]
   - gamma[9,4,4] - gamma[9,5,4] - gamma[9,6,4] - gamma[9,7,4]
   - gamma[9,8,4] - gamma[9,9,4] + lambda[9,4] + mu[9,4] + nu[9] <= 0
 u[9,5]: - gamma[9,0,5] - gamma[9,1,5] - gamma[9,2,5] - gamma[9,3,5]
   - gamma[9,4,5] - gamma[9,5,5] - gamma[9,6,5] - gamma[9,7,5]
   - gamma[9,8,5] - gamma[9,9,5] + lambda[9,5] + mu[9,5] + nu[9] <= 0
 u[9,6]: - gamma[9,0,6] - gamma[9,1,6] - gamma[9,2,6] - gamma[9,3,6]
   - gamma[9,4,6] - gamma[9,5,6] - gamma[9,6,6] - gamma[9,7,6]
   - gamma[9,8,6] - gamma[9,9,6] + lambda[9,6] + mu[9,6] + nu[9] <= 0
 u[9,7]: - gamma[9,0,7] - gamma[9,1,7] - gamma[9,2,7] - gamma[9,3,7]
   - gamma[9,4,7] - gamma[9,5,7] - gamma[9,6,7] - gamma[9,7,7]
   - gamma[9,8,7] - gamma[9,9,7] + lambda[9,7] + mu[9,7] + nu[9] <= 0
 u[9,8]: - gamma[9,0,8] - gamma[9,1,8] - gamma[9,2,8] - gamma[9,3,8]
   - gamma[9,4,8] - gamma[9,5,8] - gamma[9,6,8] - gamma[9,7,8]
   - gamma[9,8,8] - gamma[9,9,8] + lambda[9,8] + mu[9,8] + nu[9] <= 0
 u[9,9]: - gamma[9,0,9] - gamma[9,1,9] - gamma[9,2,9] - gamma[9,3,9]
   - gamma[9,4,9] - gamma[9,5,9] - gamma[9,6,9] - gamma[9,7,9]
   - gamma[9,8,9] - gamma[9,9,9] + lambda[9,9] + mu[9,9] + nu[9] <= 0
 v[0,0,0]: gamma[0,0,0] + delta[0,0,0] + epsilon[0,0] <= 0
 v[0,0,1]: 3870.026035080948 beta[0,0] + gamma[0,0,1] + delta[0,0,1]
   + epsilon[0,0] <= 0
 v[0,0,2]: 1170.30715372198 beta[0,0] + gamma[0,0,2] + delta[0,0,2]
   + epsilon[0,0] <= 0
 v[0,0,3]: 503.6616792821138 beta[0,0] + gamma[0,0,3] + delta[0,0,3]
   + epsilon[0,0] <= 0
 v[0,0,4]: 106.6339244195825 beta[0,0] + gamma[0,0,4] + delta[0,0,4]
   + epsilon[0,0] <= 0
 v[0,0,5]: 30.15230660255738 beta[0,0] + gamma[0,0,5] + delta[0,0,5]
   + epsilon[0,0] <= 0
 v[0,0,6]: 3096.684809812803 beta[0,0] + gamma[0,0,6] + delta[0,0,6]
   + epsilon[0,0] <= 0
 v[0,0,7]: 3309.704453467217 beta[0,0] + gamma[0,0,7] + delta[0,0,7]
   + epsilon[0,0] <= 0
 v[0,0,8]: 1942.370401490486 beta[0,0] + gamma[0,0,8] + delta[0,0,8]
   + epsilon[0,0] <= 0
 v[0,0,9]: 4701.22268139104 beta[0,0] + gamma[0,0,9] + delta[0,0,9]
   + epsilon[0,0] <= 0
 v[0,1,0]: 2375.688345596691 beta[0,1] + gamma[0,1,0] + delta[0,1,0]
   + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 826.5192834345412 beta[0,1] + gamma[0,1,2] + delta[0,1,2]
   + epsilon[0,1] <= 0
 v[0,1,3]: 2907.335119441443 beta[0,1] + gamma[0,1,3] + delta[0,1,3]
   + epsilon[0,1] <= 0
 v[0,1,4]: 1329.706164684987 beta[0,1] + gamma[0,1,4] + delta[0,1,4]
   + epsilon[0,1] <= 0
 v[0,1,5]: 2997.756489908481 beta[0,1] + gamma[0,1,5] + delta[0,1,5]
   + epsilon[0,1] <= 0
 v[0,1,6]: 1043.913721068225 beta[0,1] + gamma[0,1,6] + delta[0,1,6]
   + epsilon[0,1] <= 0
 v[0,1,7]: 1542.155515239143 beta[0,1] + gamma[0,1,7] + delta[0,1,7]
   + epsilon[0,1] <= 0
 v[0,1,8]: 672.6806278056603 beta[0,1] + gamma[0,1,8] + delta[0,1,8]
   + epsilon[0,1] <= 0
 v[0,1,9]: 801.0443444637759 beta[0,1] + gamma[0,1,9] + delta[0,1,9]
   + epsilon[0,1] <= 0
 v[0,2,0]: 3406.817354620825 beta[0,2] + gamma[0,2,0] + delta[0,2,0]
   + epsilon[0,2] <= 0
 v[0,2,1]: 2365.852617650562 beta[0,2] + gamma[0,2,1] + delta[0,2,1]
   + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 v[0,2,3]: 2290.002396478753 beta[0,2] + gamma[0,2,3] + delta[0,2,3]
   + epsilon[0,2] <= 0
 v[0,2,4]: 3401.879700088519 beta[0,2] + gamma[0,2,4] + delta[0,2,4]
   + epsilon[0,2] <= 0
 v[0,2,5]: 3584.356584485162 beta[0,2] + gamma[0,2,5] + delta[0,2,5]
   + epsilon[0,2] <= 0
 v[0,2,6]: 1586.792403762411 beta[0,2] + gamma[0,2,6] + delta[0,2,6]
   + epsilon[0,2] <= 0
 v[0,2,7]: 951.3847145806224 beta[0,2] + gamma[0,2,7] + delta[0,2,7]
   + epsilon[0,2] <= 0
 v[0,2,8]: 249.299309351511 beta[0,2] + gamma[0,2,8] + delta[0,2,8]
   + epsilon[0,2] <= 0
 v[0,2,9]: 3759.680745064636 beta[0,2] + gamma[0,2,9] + delta[0,2,9]
   + epsilon[0,2] <= 0
 v[0,3,0]: 244.8157681500021 beta[0,3] + gamma[0,3,0] + delta[0,3,0]
   + epsilon[0,3] <= 0
 v[0,3,1]: 1554.095358464138 beta[0,3] + gamma[0,3,1] + delta[0,3,1]
   + epsilon[0,3] <= 0
 v[0,3,2]: 1630.580966103643 beta[0,3] + gamma[0,3,2] + delta[0,3,2]
   + epsilon[0,3] <= 0
 v[0,3,3]: gamma[0,3,3] + delta[0,3,3] + epsilon[0,3] <= 0
 v[0,3,4]: 399.903783476488 beta[0,3] + gamma[0,3,4] + delta[0,3,4]
   + epsilon[0,3] <= 0
 v[0,3,5]: 312.1013448363041 beta[0,3] + gamma[0,3,5] + delta[0,3,5]
   + epsilon[0,3] <= 0
 v[0,3,6]: 3020.216076807229 beta[0,3] + gamma[0,3,6] + delta[0,3,6]
   + epsilon[0,3] <= 0
 v[0,3,7]: 2863.768379395586 beta[0,3] + gamma[0,3,7] + delta[0,3,7]
   + epsilon[0,3] <= 0
 v[0,3,8]: 519.7946528445602 beta[0,3] + gamma[0,3,8] + delta[0,3,8]
   + epsilon[0,3] <= 0
 v[0,3,9]: 3193.088784489257 beta[0,3] + gamma[0,3,9] + delta[0,3,9]
   + epsilon[0,3] <= 0
 v[0,4,0]: 254.5487918224523 beta[0,4] + gamma[0,4,0] + delta[0,4,0]
   + epsilon[0,4] <= 0
 v[0,4,1]: 4004.205678203949 beta[0,4] + gamma[0,4,1] + delta[0,4,1]
   + epsilon[0,4] <= 0
 v[0,4,2]: 1787.12530827881 beta[0,4] + gamma[0,4,2] + delta[0,4,2]
   + epsilon[0,4] <= 0
 v[0,4,3]: 688.2258938388384 beta[0,4] + gamma[0,4,3] + delta[0,4,3]
   + epsilon[0,4] <= 0
 v[0,4,4]: gamma[0,4,4] + delta[0,4,4] + epsilon[0,4] <= 0
 v[0,4,5]: 112.3318746974663 beta[0,4] + gamma[0,4,5] + delta[0,4,5]
   + epsilon[0,4] <= 0
 v[0,4,6]: 6454.951868509551 beta[0,4] + gamma[0,4,6] + delta[0,4,6]
   + epsilon[0,4] <= 0
 v[0,4,7]: 3109.529630568012 beta[0,4] + gamma[0,4,7] + delta[0,4,7]
   + epsilon[0,4] <= 0
 v[0,4,8]: 1673.220524034328 beta[0,4] + gamma[0,4,8] + delta[0,4,8]
   + epsilon[0,4] <= 0
 v[0,4,9]: 2895.825007946899 beta[0,4] + gamma[0,4,9] + delta[0,4,9]
   + epsilon[0,4] <= 0
 v[0,5,0]: 24.64166164689316 beta[0,5] + gamma[0,5,0] + delta[0,5,0]
   + epsilon[0,5] <= 0
 v[0,5,1]: 3452.962691130367 beta[0,5] + gamma[0,5,1] + delta[0,5,1]
   + epsilon[0,5] <= 0
 v[0,5,2]: 1523.110552042076 beta[0,5] + gamma[0,5,2] + delta[0,5,2]
   + epsilon[0,5] <= 0
 v[0,5,3]: 418.9655134141246 beta[0,5] + gamma[0,5,3] + delta[0,5,3]
   + epsilon[0,5] <= 0
 v[0,5,4]: 105.6213016536175 beta[0,5] + gamma[0,5,4] + delta[0,5,4]
   + epsilon[0,5] <= 0
 v[0,5,5]: gamma[0,5,5] + delta[0,5,5] + epsilon[0,5] <= 0
 v[0,5,6]: 4484.599220685022 beta[0,5] + gamma[0,5,6] + delta[0,5,6]
   + epsilon[0,5] <= 0
 v[0,5,7]: 2332.574298284483 beta[0,5] + gamma[0,5,7] + delta[0,5,7]
   + epsilon[0,5] <= 0
 v[0,5,8]: 2158.408636422978 beta[0,5] + gamma[0,5,8] + delta[0,5,8]
   + epsilon[0,5] <= 0
 v[0,5,9]: 5023.922077900266 beta[0,5] + gamma[0,5,9] + delta[0,5,9]
   + epsilon[0,5] <= 0
 v[0,6,0]: 1148.184607074712 beta[0,6] + gamma[0,6,0] + delta[0,6,0]
   + epsilon[0,6] <= 0
 v[0,6,1]: 539.1953031284435 beta[0,6] + gamma[0,6,1] + delta[0,6,1]
   + epsilon[0,6] <= 0
 v[0,6,2]: 1222.70255295892 beta[0,6] + gamma[0,6,2] + delta[0,6,2]
   + epsilon[0,6] <= 0
 v[0,6,3]: 1508.221912652303 beta[0,6] + gamma[0,6,3] + delta[0,6,3]
   + epsilon[0,6] <= 0
 v[0,6,4]: 1026.683627876996 beta[0,6] + gamma[0,6,4] + delta[0,6,4]
   + epsilon[0,6] <= 0
 v[0,6,5]: 1960.238041354198 beta[0,6] + gamma[0,6,5] + delta[0,6,5]
   + epsilon[0,6] <= 0
 v[0,6,6]: gamma[0,6,6] + delta[0,6,6] + epsilon[0,6] <= 0
 v[0,6,7]: 698.7382523118695 beta[0,6] + gamma[0,6,7] + delta[0,6,7]
   + epsilon[0,6] <= 0
 v[0,6,8]: 1484.743661941874 beta[0,6] + gamma[0,6,8] + delta[0,6,8]
   + epsilon[0,6] <= 0
 v[0,6,9]: 6.498537008678615 beta[0,6] + gamma[0,6,9] + delta[0,6,9]
   + epsilon[0,6] <= 0
 v[0,7,0]: 6365.782654543784 beta[0,7] + gamma[0,7,0] + delta[0,7,0]
   + epsilon[0,7] <= 0
 v[0,7,1]: 1320.396100037685 beta[0,7] + gamma[0,7,1] + delta[0,7,1]
   + epsilon[0,7] <= 0
 v[0,7,2]: 1168.183115371031 beta[0,7] + gamma[0,7,2] + delta[0,7,2]
   + epsilon[0,7] <= 0
 v[0,7,3]: 4987.781601658012 beta[0,7] + gamma[0,7,3] + delta[0,7,3]
   + epsilon[0,7] <= 0
 v[0,7,4]: 3930.689939510291 beta[0,7] + gamma[0,7,4] + delta[0,7,4]
   + epsilon[0,7] <= 0
 v[0,7,5]: 2983.64056877798 beta[0,7] + gamma[0,7,5] + delta[0,7,5]
   + epsilon[0,7] <= 0
 v[0,7,6]: 2696.088448209302 beta[0,7] + gamma[0,7,6] + delta[0,7,6]
   + epsilon[0,7] <= 0
 v[0,7,7]: gamma[0,7,7] + delta[0,7,7] + epsilon[0,7] <= 0
 v[0,7,8]: 2618.70531704454 beta[0,7] + gamma[0,7,8] + delta[0,7,8]
   + epsilon[0,7] <= 0
 v[0,7,9]: 2331.669867637379 beta[0,7] + gamma[0,7,9] + delta[0,7,9]
   + epsilon[0,7] <= 0
 v[0,8,0]: 3128.173302833535 beta[0,8] + gamma[0,8,0] + delta[0,8,0]
   + epsilon[0,8] <= 0
 v[0,8,1]: 1512.108386875047 beta[0,8] + gamma[0,8,1] + delta[0,8,1]
   + epsilon[0,8] <= 0
 v[0,8,2]: 434.2977088405181 beta[0,8] + gamma[0,8,2] + delta[0,8,2]
   + epsilon[0,8] <= 0
 v[0,8,3]: 1363.749220689497 beta[0,8] + gamma[0,8,3] + delta[0,8,3]
   + epsilon[0,8] <= 0
 v[0,8,4]: 1646.490794025562 beta[0,8] + gamma[0,8,4] + delta[0,8,4]
   + epsilon[0,8] <= 0
 v[0,8,5]: 3098.197698462424 beta[0,8] + gamma[0,8,5] + delta[0,8,5]
   + epsilon[0,8] <= 0
 v[0,8,6]: 2357.493804196921 beta[0,8] + gamma[0,8,6] + delta[0,8,6]
   + epsilon[0,8] <= 0
 v[0,8,7]: 2762.992724441164 beta[0,8] + gamma[0,8,7] + delta[0,8,7]
   + epsilon[0,8] <= 0
 v[0,8,8]: gamma[0,8,8] + delta[0,8,8] + epsilon[0,8] <= 0
 v[0,8,9]: 2865.540029431045 beta[0,8] + gamma[0,8,9] + delta[0,8,9]
   + epsilon[0,8] <= 0
 v[0,9,0]: 3681.432329055701 beta[0,9] + gamma[0,9,0] + delta[0,9,0]
   + epsilon[0,9] <= 0
 v[0,9,1]: 2700.590616515457 beta[0,9] + gamma[0,9,1] + delta[0,9,1]
   + epsilon[0,9] <= 0
 v[0,9,2]: 3063.81975459874 beta[0,9] + gamma[0,9,2] + delta[0,9,2]
   + epsilon[0,9] <= 0
 v[0,9,3]: 5592.0864978291 beta[0,9] + gamma[0,9,3] + delta[0,9,3]
   + epsilon[0,9] <= 0
 v[0,9,4]: 5326.607888987719 beta[0,9] + gamma[0,9,4] + delta[0,9,4]
   + epsilon[0,9] <= 0
 v[0,9,5]: 4381.638823440833 beta[0,9] + gamma[0,9,5] + delta[0,9,5]
   + epsilon[0,9] <= 0
 v[0,9,6]: 26.29192564252354 beta[0,9] + gamma[0,9,6] + delta[0,9,6]
   + epsilon[0,9] <= 0
 v[0,9,7]: 2331.717178175971 beta[0,9] + gamma[0,9,7] + delta[0,9,7]
   + epsilon[0,9] <= 0
 v[0,9,8]: 1628.23226974173 beta[0,9] + gamma[0,9,8] + delta[0,9,8]
   + epsilon[0,9] <= 0
 v[0,9,9]: gamma[0,9,9] + delta[0,9,9] + epsilon[0,9] <= 0
 v[1,0,0]: gamma[1,0,0] + delta[1,0,0] + epsilon[1,0] <= 0
 v[1,0,1]: 2187.32951769123 beta[1,0] + gamma[1,0,1] + delta[1,0,1]
   + epsilon[1,0] <= 0
 v[1,0,2]: 3629.090301559848 beta[1,0] + gamma[1,0,2] + delta[1,0,2]
   + epsilon[1,0] <= 0
 v[1,0,3]: 813.3104219964154 beta[1,0] + gamma[1,0,3] + delta[1,0,3]
   + epsilon[1,0] <= 0
 v[1,0,4]: 101.0218330650952 beta[1,0] + gamma[1,0,4] + delta[1,0,4]
   + epsilon[1,0] <= 0
 v[1,0,5]: 69.18279878155902 beta[1,0] + gamma[1,0,5] + delta[1,0,5]
   + epsilon[1,0] <= 0
 v[1,0,6]: 6125.006026945133 beta[1,0] + gamma[1,0,6] + delta[1,0,6]
   + epsilon[1,0] <= 0
 v[1,0,7]: 6029.766300814853 beta[1,0] + gamma[1,0,7] + delta[1,0,7]
   + epsilon[1,0] <= 0
 v[1,0,8]: 2840.71836225207 beta[1,0] + gamma[1,0,8] + delta[1,0,8]
   + epsilon[1,0] <= 0
 v[1,0,9]: 3058.360678017436 beta[1,0] + gamma[1,0,9] + delta[1,0,9]
   + epsilon[1,0] <= 0
 v[1,1,0]: 1914.311038315645 beta[1,1] + gamma[1,1,0] + delta[1,1,0]
   + epsilon[1,1] <= 0
 v[1,1,1]: gamma[1,1,1] + delta[1,1,1] + epsilon[1,1] <= 0
 v[1,1,2]: 883.1096135015304 beta[1,1] + gamma[1,1,2] + delta[1,1,2]
   + epsilon[1,1] <= 0
 v[1,1,3]: 2242.837408656146 beta[1,1] + gamma[1,1,3] + delta[1,1,3]
   + epsilon[1,1] <= 0
 v[1,1,4]: 4885.408063513409 beta[1,1] + gamma[1,1,4] + delta[1,1,4]
   + epsilon[1,1] <= 0
 v[1,1,5]: 3790.610445007482 beta[1,1] + gamma[1,1,5] + delta[1,1,5]
   + epsilon[1,1] <= 0
 v[1,1,6]: 1706.99842212319 beta[1,1] + gamma[1,1,6] + delta[1,1,6]
   + epsilon[1,1] <= 0
 v[1,1,7]: 2785.153797365277 beta[1,1] + gamma[1,1,7] + delta[1,1,7]
   + epsilon[1,1] <= 0
 v[1,1,8]: 934.1118709673776 beta[1,1] + gamma[1,1,8] + delta[1,1,8]
   + epsilon[1,1] <= 0
 v[1,1,9]: 1273.811248163906 beta[1,1] + gamma[1,1,9] + delta[1,1,9]
   + epsilon[1,1] <= 0
 v[1,2,0]: 1382.854917762706 beta[1,2] + gamma[1,2,0] + delta[1,2,0]
   + epsilon[1,2] <= 0
 v[1,2,1]: 1151.511771607619 beta[1,2] + gamma[1,2,1] + delta[1,2,1]
   + epsilon[1,2] <= 0
 v[1,2,2]: gamma[1,2,2] + delta[1,2,2] + epsilon[1,2] <= 0
 v[1,2,3]: 484.3959572166975 beta[1,2] + gamma[1,2,3] + delta[1,2,3]
   + epsilon[1,2] <= 0
 v[1,2,4]: 1563.214634773698 beta[1,2] + gamma[1,2,4] + delta[1,2,4]
   + epsilon[1,2] <= 0
 v[1,2,5]: 916.0822037637515 beta[1,2] + gamma[1,2,5] + delta[1,2,5]
   + epsilon[1,2] <= 0
 v[1,2,6]: 1723.964066045872 beta[1,2] + gamma[1,2,6] + delta[1,2,6]
   + epsilon[1,2] <= 0
 v[1,2,7]: 761.9988453330432 beta[1,2] + gamma[1,2,7] + delta[1,2,7]
   + epsilon[1,2] <= 0
 v[1,2,8]: 246.7805479838349 beta[1,2] + gamma[1,2,8] + delta[1,2,8]
   + epsilon[1,2] <= 0
 v[1,2,9]: 1654.443343915749 beta[1,2] + gamma[1,2,9] + delta[1,2,9]
   + epsilon[1,2] <= 0
 v[1,3,0]: 315.7435388293446 beta[1,3] + gamma[1,3,0] + delta[1,3,0]
   + epsilon[1,3] <= 0
 v[1,3,1]: 784.5498054612057 beta[1,3] + gamma[1,3,1] + delta[1,3,1]
   + epsilon[1,3] <= 0
 v[1,3,2]: 457.5381546873854 beta[1,3] + gamma[1,3,2] + delta[1,3,2]
   + epsilon[1,3] <= 0
 v[1,3,3]: gamma[1,3,3] + delta[1,3,3] + epsilon[1,3] <= 0
 v[1,3,4]: 133.4899804303137 beta[1,3] + gamma[1,3,4] + delta[1,3,4]
   + epsilon[1,3] <= 0
 v[1,3,5]: 169.8598796950947 beta[1,3] + gamma[1,3,5] + delta[1,3,5]
   + epsilon[1,3] <= 0
 v[1,3,6]: 1349.640248398168 beta[1,3] + gamma[1,3,6] + delta[1,3,6]
   + epsilon[1,3] <= 0
 v[1,3,7]: 1901.606473008364 beta[1,3] + gamma[1,3,7] + delta[1,3,7]
   + epsilon[1,3] <= 0
 v[1,3,8]: 776.4896870128167 beta[1,3] + gamma[1,3,8] + delta[1,3,8]
   + epsilon[1,3] <= 0
 v[1,3,9]: 956.4835946537349 beta[1,3] + gamma[1,3,9] + delta[1,3,9]
   + epsilon[1,3] <= 0
 v[1,4,0]: 44.09136367244794 beta[1,4] + gamma[1,4,0] + delta[1,4,0]
   + epsilon[1,4] <= 0
 v[1,4,1]: 2122.346217083808 beta[1,4] + gamma[1,4,1] + delta[1,4,1]
   + epsilon[1,4] <= 0
 v[1,4,2]: 1102.771508548173 beta[1,4] + gamma[1,4,2] + delta[1,4,2]
   + epsilon[1,4] <= 0
 v[1,4,3]: 169.2702359295827 beta[1,4] + gamma[1,4,3] + delta[1,4,3]
   + epsilon[1,4] <= 0
 v[1,4,4]: gamma[1,4,4] + delta[1,4,4] + epsilon[1,4] <= 0
 v[1,4,5]: 66.46043817169006 beta[1,4] + gamma[1,4,5] + delta[1,4,5]
   + epsilon[1,4] <= 0
 v[1,4,6]: 1436.878277151181 beta[1,4] + gamma[1,4,6] + delta[1,4,6]
   + epsilon[1,4] <= 0
 v[1,4,7]: 1910.3605387677 beta[1,4] + gamma[1,4,7] + delta[1,4,7]
   + epsilon[1,4] <= 0
 v[1,4,8]: 1407.778126679765 beta[1,4] + gamma[1,4,8] + delta[1,4,8]
   + epsilon[1,4] <= 0
 v[1,4,9]: 2791.159371072006 beta[1,4] + gamma[1,4,9] + delta[1,4,9]
   + epsilon[1,4] <= 0
 v[1,5,0]: 31.42470053051888 beta[1,5] + gamma[1,5,0] + delta[1,5,0]
   + epsilon[1,5] <= 0
 v[1,5,1]: 2605.898742423949 beta[1,5] + gamma[1,5,1] + delta[1,5,1]
   + epsilon[1,5] <= 0
 v[1,5,2]: 2145.913891961862 beta[1,5] + gamma[1,5,2] + delta[1,5,2]
   + epsilon[1,5] <= 0
 v[1,5,3]: 651.1161399471181 beta[1,5] + gamma[1,5,3] + delta[1,5,3]
   + epsilon[1,5] <= 0
 v[1,5,4]: 62.63070756731946 beta[1,5] + gamma[1,5,4] + delta[1,5,4]
   + epsilon[1,5] <= 0
 v[1,5,5]: gamma[1,5,5] + delta[1,5,5] + epsilon[1,5] <= 0
 v[1,5,6]: 2130.574086105224 beta[1,5] + gamma[1,5,6] + delta[1,5,6]
   + epsilon[1,5] <= 0
 v[1,5,7]: 3184.750711754107 beta[1,5] + gamma[1,5,7] + delta[1,5,7]
   + epsilon[1,5] <= 0
 v[1,5,8]: 1962.463319299846 beta[1,5] + gamma[1,5,8] + delta[1,5,8]
   + epsilon[1,5] <= 0
 v[1,5,9]: 4006.81754428061 beta[1,5] + gamma[1,5,9] + delta[1,5,9]
   + epsilon[1,5] <= 0
 v[1,6,0]: 1575.304229957168 beta[1,6] + gamma[1,6,0] + delta[1,6,0]
   + epsilon[1,6] <= 0
 v[1,6,1]: 887.8340074584326 beta[1,6] + gamma[1,6,1] + delta[1,6,1]
   + epsilon[1,6] <= 0
 v[1,6,2]: 1117.552675959074 beta[1,6] + gamma[1,6,2] + delta[1,6,2]
   + epsilon[1,6] <= 0
 v[1,6,3]: 1468.12137289096 beta[1,6] + gamma[1,6,3] + delta[1,6,3]
   + epsilon[1,6] <= 0
 v[1,6,4]: 1983.662408612589 beta[1,6] + gamma[1,6,4] + delta[1,6,4]
   + epsilon[1,6] <= 0
 v[1,6,5]: 2382.616957830182 beta[1,6] + gamma[1,6,5] + delta[1,6,5]
   + epsilon[1,6] <= 0
 v[1,6,6]: gamma[1,6,6] + delta[1,6,6] + epsilon[1,6] <= 0
 v[1,6,7]: 532.4605509886147 beta[1,6] + gamma[1,6,7] + delta[1,6,7]
   + epsilon[1,6] <= 0
 v[1,6,8]: 588.5108409266858 beta[1,6] + gamma[1,6,8] + delta[1,6,8]
   + epsilon[1,6] <= 0
 v[1,6,9]: 5.261394822037821 beta[1,6] + gamma[1,6,9] + delta[1,6,9]
   + epsilon[1,6] <= 0
 v[1,7,0]: 5190.798106934419 beta[1,7] + gamma[1,7,0] + delta[1,7,0]
   + epsilon[1,7] <= 0
 v[1,7,1]: 1490.049130198875 beta[1,7] + gamma[1,7,1] + delta[1,7,1]
   + epsilon[1,7] <= 0
 v[1,7,2]: 2318.906096272298 beta[1,7] + gamma[1,7,2] + delta[1,7,2]
   + epsilon[1,7] <= 0
 v[1,7,3]: 2189.395695499139 beta[1,7] + gamma[1,7,3] + delta[1,7,3]
   + epsilon[1,7] <= 0
 v[1,7,4]: 3372.383221971231 beta[1,7] + gamma[1,7,4] + delta[1,7,4]
   + epsilon[1,7] <= 0
 v[1,7,5]: 5926.754858866276 beta[1,7] + gamma[1,7,5] + delta[1,7,5]
   + epsilon[1,7] <= 0
 v[1,7,6]: 1326.646617651896 beta[1,7] + gamma[1,7,6] + delta[1,7,6]
   + epsilon[1,7] <= 0
 v[1,7,7]: gamma[1,7,7] + delta[1,7,7] + epsilon[1,7] <= 0
 v[1,7,8]: 2173.067638545975 beta[1,7] + gamma[1,7,8] + delta[1,7,8]
   + epsilon[1,7] <= 0
 v[1,7,9]: 2744.910425267672 beta[1,7] + gamma[1,7,9] + delta[1,7,9]
   + epsilon[1,7] <= 0
 v[1,8,0]: 1923.536721993297 beta[1,8] + gamma[1,8,0] + delta[1,8,0]
   + epsilon[1,8] <= 0
 v[1,8,1]: 848.7199266587268 beta[1,8] + gamma[1,8,1] + delta[1,8,1]
   + epsilon[1,8] <= 0
 v[1,8,2]: 477.8316446538374 beta[1,8] + gamma[1,8,2] + delta[1,8,2]
   + epsilon[1,8] <= 0
 v[1,8,3]: 2277.515084508714 beta[1,8] + gamma[1,8,3] + delta[1,8,3]
   + epsilon[1,8] <= 0
 v[1,8,4]: 3354.977763772072 beta[1,8] + gamma[1,8,4] + delta[1,8,4]
   + epsilon[1,8] <= 0
 v[1,8,5]: 3541.833942476251 beta[1,8] + gamma[1,8,5] + delta[1,8,5]
   + epsilon[1,8] <= 0
 v[1,8,6]: 1454.844747752973 beta[1,8] + gamma[1,8,6] + delta[1,8,6]
   + epsilon[1,8] <= 0
 v[1,8,7]: 1543.598417131297 beta[1,8] + gamma[1,8,7] + delta[1,8,7]
   + epsilon[1,8] <= 0
 v[1,8,8]: gamma[1,8,8] + delta[1,8,8] + epsilon[1,8] <= 0
 v[1,8,9]: 4125.14515260707 beta[1,8] + gamma[1,8,9] + delta[1,8,9]
   + epsilon[1,8] <= 0
 v[1,9,0]: 2987.099307319441 beta[1,9] + gamma[1,9,0] + delta[1,9,0]
   + epsilon[1,9] <= 0
 v[1,9,1]: 677.1486287522281 beta[1,9] + gamma[1,9,1] + delta[1,9,1]
   + epsilon[1,9] <= 0
 v[1,9,2]: 1570.471507885639 beta[1,9] + gamma[1,9,2] + delta[1,9,2]
   + epsilon[1,9] <= 0
 v[1,9,3]: 2036.137353970258 beta[1,9] + gamma[1,9,3] + delta[1,9,3]
   + epsilon[1,9] <= 0
 v[1,9,4]: 1795.564479639213 beta[1,9] + gamma[1,9,4] + delta[1,9,4]
   + epsilon[1,9] <= 0
 v[1,9,5]: 2016.200771108065 beta[1,9] + gamma[1,9,5] + delta[1,9,5]
   + epsilon[1,9] <= 0
 v[1,9,6]: 7.116168750551548 beta[1,9] + gamma[1,9,6] + delta[1,9,6]
   + epsilon[1,9] <= 0
 v[1,9,7]: 783.9913519619041 beta[1,9] + gamma[1,9,7] + delta[1,9,7]
   + epsilon[1,9] <= 0
 v[1,9,8]: 1588.337928904559 beta[1,9] + gamma[1,9,8] + delta[1,9,8]
   + epsilon[1,9] <= 0
 v[1,9,9]: gamma[1,9,9] + delta[1,9,9] + epsilon[1,9] <= 0
 v[2,0,0]: gamma[2,0,0] + delta[2,0,0] + epsilon[2,0] <= 0
 v[2,0,1]: 2923.47792551402 beta[2,0] + gamma[2,0,1] + delta[2,0,1]
   + epsilon[2,0] <= 0
 v[2,0,2]: 1635.921499597475 beta[2,0] + gamma[2,0,2] + delta[2,0,2]
   + epsilon[2,0] <= 0
 v[2,0,3]: 910.7342124088085 beta[2,0] + gamma[2,0,3] + delta[2,0,3]
   + epsilon[2,0] <= 0
 v[2,0,4]: 154.3053870962829 beta[2,0] + gamma[2,0,4] + delta[2,0,4]
   + epsilon[2,0] <= 0
 v[2,0,5]: 81.21760284086732 beta[2,0] + gamma[2,0,5] + delta[2,0,5]
   + epsilon[2,0] <= 0
 v[2,0,6]: 5409.360017671171 beta[2,0] + gamma[2,0,6] + delta[2,0,6]
   + epsilon[2,0] <= 0
 v[2,0,7]: 5979.439970468467 beta[2,0] + gamma[2,0,7] + delta[2,0,7]
   + epsilon[2,0] <= 0
 v[2,0,8]: 3226.959265893815 beta[2,0] + gamma[2,0,8] + delta[2,0,8]
   + epsilon[2,0] <= 0
 v[2,0,9]: 7046.512256886525 beta[2,0] + gamma[2,0,9] + delta[2,0,9]
   + epsilon[2,0] <= 0
 v[2,1,0]: 2697.334467378802 beta[2,1] + gamma[2,1,0] + delta[2,1,0]
   + epsilon[2,1] <= 0
 v[2,1,1]: gamma[2,1,1] + delta[2,1,1] + epsilon[2,1] <= 0
 v[2,1,2]: 1656.58500227162 beta[2,1] + gamma[2,1,2] + delta[2,1,2]
   + epsilon[2,1] <= 0
 v[2,1,3]: 1802.236311876137 beta[2,1] + gamma[2,1,3] + delta[2,1,3]
   + epsilon[2,1] <= 0
 v[2,1,4]: 2714.612627676679 beta[2,1] + gamma[2,1,4] + delta[2,1,4]
   + epsilon[2,1] <= 0
 v[2,1,5]: 3060.652088749507 beta[2,1] + gamma[2,1,5] + delta[2,1,5]
   + epsilon[2,1] <= 0
 v[2,1,6]: 637.8863910044912 beta[2,1] + gamma[2,1,6] + delta[2,1,6]
   + epsilon[2,1] <= 0
 v[2,1,7]: 1731.835069506326 beta[2,1] + gamma[2,1,7] + delta[2,1,7]
   + epsilon[2,1] <= 0
 v[2,1,8]: 1066.346895254466 beta[2,1] + gamma[2,1,8] + delta[2,1,8]
   + epsilon[2,1] <= 0
 v[2,1,9]: 1631.056700159953 beta[2,1] + gamma[2,1,9] + delta[2,1,9]
   + epsilon[2,1] <= 0
 v[2,2,0]: 1130.31178914973 beta[2,2] + gamma[2,2,0] + delta[2,2,0]
   + epsilon[2,2] <= 0
 v[2,2,1]: 1188.10290924052 beta[2,2] + gamma[2,2,1] + delta[2,2,1]
   + epsilon[2,2] <= 0
 v[2,2,2]: gamma[2,2,2] + delta[2,2,2] + epsilon[2,2] <= 0
 v[2,2,3]: 695.437527845769 beta[2,2] + gamma[2,2,3] + delta[2,2,3]
   + epsilon[2,2] <= 0
 v[2,2,4]: 1685.08999537336 beta[2,2] + gamma[2,2,4] + delta[2,2,4]
   + epsilon[2,2] <= 0
 v[2,2,5]: 1512.968284450181 beta[2,2] + gamma[2,2,5] + delta[2,2,5]
   + epsilon[2,2] <= 0
 v[2,2,6]: 795.2862057597631 beta[2,2] + gamma[2,2,6] + delta[2,2,6]
   + epsilon[2,2] <= 0
 v[2,2,7]: 859.7177634419101 beta[2,2] + gamma[2,2,7] + delta[2,2,7]
   + epsilon[2,2] <= 0
 v[2,2,8]: 232.2997584173099 beta[2,2] + gamma[2,2,8] + delta[2,2,8]
   + epsilon[2,2] <= 0
 v[2,2,9]: 2048.663969907535 beta[2,2] + gamma[2,2,9] + delta[2,2,9]
   + epsilon[2,2] <= 0
 v[2,3,0]: 143.7675309280787 beta[2,3] + gamma[2,3,0] + delta[2,3,0]
   + epsilon[2,3] <= 0
 v[2,3,1]: 1749.574026185072 beta[2,3] + gamma[2,3,1] + delta[2,3,1]
   + epsilon[2,3] <= 0
 v[2,3,2]: 688.4379956214432 beta[2,3] + gamma[2,3,2] + delta[2,3,2]
   + epsilon[2,3] <= 0
 v[2,3,3]: gamma[2,3,3] + delta[2,3,3] + epsilon[2,3] <= 0
 v[2,3,4]: 260.2735781798842 beta[2,3] + gamma[2,3,4] + delta[2,3,4]
   + epsilon[2,3] <= 0
 v[2,3,5]: 357.1028429747111 beta[2,3] + gamma[2,3,5] + delta[2,3,5]
   + epsilon[2,3] <= 0
 v[2,3,6]: 1920.955027430695 beta[2,3] + gamma[2,3,6] + delta[2,3,6]
   + epsilon[2,3] <= 0
 v[2,3,7]: 1699.699777859638 beta[2,3] + gamma[2,3,7] + delta[2,3,7]
   + epsilon[2,3] <= 0
 v[2,3,8]: 553.2279977228129 beta[2,3] + gamma[2,3,8] + delta[2,3,8]
   + epsilon[2,3] <= 0
 v[2,3,9]: 1949.257273431724 beta[2,3] + gamma[2,3,9] + delta[2,3,9]
   + epsilon[2,3] <= 0
 v[2,4,0]: 154.1343125838028 beta[2,4] + gamma[2,4,0] + delta[2,4,0]
   + epsilon[2,4] <= 0
 v[2,4,1]: 3142.745180522283 beta[2,4] + gamma[2,4,1] + delta[2,4,1]
   + epsilon[2,4] <= 0
 v[2,4,2]: 2072.766232554187 beta[2,4] + gamma[2,4,2] + delta[2,4,2]
   + epsilon[2,4] <= 0
 v[2,4,3]: 461.2394603100202 beta[2,4] + gamma[2,4,3] + delta[2,4,3]
   + epsilon[2,4] <= 0
 v[2,4,4]: gamma[2,4,4] + delta[2,4,4] + epsilon[2,4] <= 0
 v[2,4,5]: 70.31367714799787 beta[2,4] + gamma[2,4,5] + delta[2,4,5]
   + epsilon[2,4] <= 0
 v[2,4,6]: 3869.458219880374 beta[2,4] + gamma[2,4,6] + delta[2,4,6]
   + epsilon[2,4] <= 0
 v[2,4,7]: 2687.033844126893 beta[2,4] + gamma[2,4,7] + delta[2,4,7]
   + epsilon[2,4] <= 0
 v[2,4,8]: 1351.24740928853 beta[2,4] + gamma[2,4,8] + delta[2,4,8]
   + epsilon[2,4] <= 0
 v[2,4,9]: 3019.293958742109 beta[2,4] + gamma[2,4,9] + delta[2,4,9]
   + epsilon[2,4] <= 0
 v[2,5,0]: 37.20401399649348 beta[2,5] + gamma[2,5,0] + delta[2,5,0]
   + epsilon[2,5] <= 0
 v[2,5,1]: 2082.586272964985 beta[2,5] + gamma[2,5,1] + delta[2,5,1]
   + epsilon[2,5] <= 0
 v[2,5,2]: 2088.71027861214 beta[2,5] + gamma[2,5,2] + delta[2,5,2]
   + epsilon[2,5] <= 0
 v[2,5,3]: 442.8363321474174 beta[2,5] + gamma[2,5,3] + delta[2,5,3]
   + epsilon[2,5] <= 0
 v[2,5,4]: 136.5539123679538 beta[2,5] + gamma[2,5,4] + delta[2,5,4]
   + epsilon[2,5] <= 0
 v[2,5,5]: gamma[2,5,5] + delta[2,5,5] + epsilon[2,5] <= 0
 v[2,5,6]: 2536.996650279392 beta[2,5] + gamma[2,5,6] + delta[2,5,6]
   + epsilon[2,5] <= 0
 v[2,5,7]: 2907.309731505086 beta[2,5] + gamma[2,5,7] + delta[2,5,7]
   + epsilon[2,5] <= 0
 v[2,5,8]: 1453.543729321354 beta[2,5] + gamma[2,5,8] + delta[2,5,8]
   + epsilon[2,5] <= 0
 v[2,5,9]: 4726.516344443321 beta[2,5] + gamma[2,5,9] + delta[2,5,9]
   + epsilon[2,5] <= 0
 v[2,6,0]: 3328.496477077857 beta[2,6] + gamma[2,6,0] + delta[2,6,0]
   + epsilon[2,6] <= 0
 v[2,6,1]: 1047.820863333645 beta[2,6] + gamma[2,6,1] + delta[2,6,1]
   + epsilon[2,6] <= 0
 v[2,6,2]: 1506.057526974149 beta[2,6] + gamma[2,6,2] + delta[2,6,2]
   + epsilon[2,6] <= 0
 v[2,6,3]: 1607.820537333284 beta[2,6] + gamma[2,6,3] + delta[2,6,3]
   + epsilon[2,6] <= 0
 v[2,6,4]: 1777.197878206649 beta[2,6] + gamma[2,6,4] + delta[2,6,4]
   + epsilon[2,6] <= 0
 v[2,6,5]: 3244.102677488952 beta[2,6] + gamma[2,6,5] + delta[2,6,5]
   + epsilon[2,6] <= 0
 v[2,6,6]: gamma[2,6,6] + delta[2,6,6] + epsilon[2,6] <= 0
 v[2,6,7]: 1134.373726667051 beta[2,6] + gamma[2,6,7] + delta[2,6,7]
   + epsilon[2,6] <= 0
 v[2,6,8]: 1394.741162663459 beta[2,6] + gamma[2,6,8] + delta[2,6,8]
   + epsilon[2,6] <= 0
 v[2,6,9]: 9.331955556029492 beta[2,6] + gamma[2,6,9] + delta[2,6,9]
   + epsilon[2,6] <= 0
 v[2,7,0]: 3627.451705147775 beta[2,7] + gamma[2,7,0] + delta[2,7,0]
   + epsilon[2,7] <= 0
 v[2,7,1]: 2065.855383631183 beta[2,7] + gamma[2,7,1] + delta[2,7,1]
   + epsilon[2,7] <= 0
 v[2,7,2]: 1532.833752393671 beta[2,7] + gamma[2,7,2] + delta[2,7,2]
   + epsilon[2,7] <= 0
 v[2,7,3]: 2891.354043678962 beta[2,7] + gamma[2,7,3] + delta[2,7,3]
   + epsilon[2,7] <= 0
 v[2,7,4]: 2895.801598006956 beta[2,7] + gamma[2,7,4] + delta[2,7,4]
   + epsilon[2,7] <= 0
 v[2,7,5]: 1999.653443624281 beta[2,7] + gamma[2,7,5] + delta[2,7,5]
   + epsilon[2,7] <= 0
 v[2,7,6]: 1301.263285154166 beta[2,7] + gamma[2,7,6] + delta[2,7,6]
   + epsilon[2,7] <= 0
 v[2,7,7]: gamma[2,7,7] + delta[2,7,7] + epsilon[2,7] <= 0
 v[2,7,8]: 2436.323210603106 beta[2,7] + gamma[2,7,8] + delta[2,7,8]
   + epsilon[2,7] <= 0
 v[2,7,9]: 1017.344008787003 beta[2,7] + gamma[2,7,9] + delta[2,7,9]
   + epsilon[2,7] <= 0
 v[2,8,0]: 2600.699126643872 beta[2,8] + gamma[2,8,0] + delta[2,8,0]
   + epsilon[2,8] <= 0
 v[2,8,1]: 1622.050620526517 beta[2,8] + gamma[2,8,1] + delta[2,8,1]
   + epsilon[2,8] <= 0
 v[2,8,2]: 361.6997731859004 beta[2,8] + gamma[2,8,2] + delta[2,8,2]
   + epsilon[2,8] <= 0
 v[2,8,3]: 1506.218186130015 beta[2,8] + gamma[2,8,3] + delta[2,8,3]
   + epsilon[2,8] <= 0
 v[2,8,4]: 1239.146010761215 beta[2,8] + gamma[2,8,4] + delta[2,8,4]
   + epsilon[2,8] <= 0
 v[2,8,5]: 1881.318090061797 beta[2,8] + gamma[2,8,5] + delta[2,8,5]
   + epsilon[2,8] <= 0
 v[2,8,6]: 1882.069140556595 beta[2,8] + gamma[2,8,6] + delta[2,8,6]
   + epsilon[2,8] <= 0
 v[2,8,7]: 1198.897391486791 beta[2,8] + gamma[2,8,7] + delta[2,8,7]
   + epsilon[2,8] <= 0
 v[2,8,8]: gamma[2,8,8] + delta[2,8,8] + epsilon[2,8] <= 0
 v[2,8,9]: 1819.123902424711 beta[2,8] + gamma[2,8,9] + delta[2,8,9]
   + epsilon[2,8] <= 0
 v[2,9,0]: 6537.626867278094 beta[2,9] + gamma[2,9,0] + delta[2,9,0]
   + epsilon[2,9] <= 0
 v[2,9,1]: 2024.38276562609 beta[2,9] + gamma[2,9,1] + delta[2,9,1]
   + epsilon[2,9] <= 0
 v[2,9,2]: 1882.285770881166 beta[2,9] + gamma[2,9,2] + delta[2,9,2]
   + epsilon[2,9] <= 0
 v[2,9,3]: 2441.766700612025 beta[2,9] + gamma[2,9,3] + delta[2,9,3]
   + epsilon[2,9] <= 0
 v[2,9,4]: 3132.934915623242 beta[2,9] + gamma[2,9,4] + delta[2,9,4]
   + epsilon[2,9] <= 0
 v[2,9,5]: 4580.518083966218 beta[2,9] + gamma[2,9,5] + delta[2,9,5]
   + epsilon[2,9] <= 0
 v[2,9,6]: 17.74168024829745 beta[2,9] + gamma[2,9,6] + delta[2,9,6]
   + epsilon[2,9] <= 0
 v[2,9,7]: 1185.48420154402 beta[2,9] + gamma[2,9,7] + delta[2,9,7]
   + epsilon[2,9] <= 0
 v[2,9,8]: 1866.910386318885 beta[2,9] + gamma[2,9,8] + delta[2,9,8]
   + epsilon[2,9] <= 0
 v[2,9,9]: gamma[2,9,9] + delta[2,9,9] + epsilon[2,9] <= 0
 v[3,0,0]: gamma[3,0,0] + delta[3,0,0] + epsilon[3,0] <= 0
 v[3,0,1]: 4572.794255122402 beta[3,0] + gamma[3,0,1] + delta[3,0,1]
   + epsilon[3,0] <= 0
 v[3,0,2]: 1946.703784145207 beta[3,0] + gamma[3,0,2] + delta[3,0,2]
   + epsilon[3,0] <= 0
 v[3,0,3]: 737.2438967893161 beta[3,0] + gamma[3,0,3] + delta[3,0,3]
   + epsilon[3,0] <= 0
 v[3,0,4]: 244.1522670561566 beta[3,0] + gamma[3,0,4] + delta[3,0,4]
   + epsilon[3,0] <= 0
 v[3,0,5]: 30.63362784833252 beta[3,0] + gamma[3,0,5] + delta[3,0,5]
   + epsilon[3,0] <= 0
 v[3,0,6]: 4285.771244179145 beta[3,0] + gamma[3,0,6] + delta[3,0,6]
   + epsilon[3,0] <= 0
 v[3,0,7]: 6191.458842212887 beta[3,0] + gamma[3,0,7] + delta[3,0,7]
   + epsilon[3,0] <= 0
 v[3,0,8]: 3355.701386956742 beta[3,0] + gamma[3,0,8] + delta[3,0,8]
   + epsilon[3,0] <= 0
 v[3,0,9]: 5414.47232519466 beta[3,0] + gamma[3,0,9] + delta[3,0,9]
   + epsilon[3,0] <= 0
 v[3,1,0]: 3267.340049538999 beta[3,1] + gamma[3,1,0] + delta[3,1,0]
   + epsilon[3,1] <= 0
 v[3,1,1]: gamma[3,1,1] + delta[3,1,1] + epsilon[3,1] <= 0
 v[3,1,2]: 1439.615392489828 beta[3,1] + gamma[3,1,2] + delta[3,1,2]
   + epsilon[3,1] <= 0
 v[3,1,3]: 2642.608363757437 beta[3,1] + gamma[3,1,3] + delta[3,1,3]
   + epsilon[3,1] <= 0
 v[3,1,4]: 1396.480130970392 beta[3,1] + gamma[3,1,4] + delta[3,1,4]
   + epsilon[3,1] <= 0
 v[3,1,5]: 2727.692281369059 beta[3,1] + gamma[3,1,5] + delta[3,1,5]
   + epsilon[3,1] <= 0
 v[3,1,6]: 1266.984817576494 beta[3,1] + gamma[3,1,6] + delta[3,1,6]
   + epsilon[3,1] <= 0
 v[3,1,7]: 1496.045719700173 beta[3,1] + gamma[3,1,7] + delta[3,1,7]
   + epsilon[3,1] <= 0
 v[3,1,8]: 1262.843418550899 beta[3,1] + gamma[3,1,8] + delta[3,1,8]
   + epsilon[3,1] <= 0
 v[3,1,9]: 567.4293749871292 beta[3,1] + gamma[3,1,9] + delta[3,1,9]
   + epsilon[3,1] <= 0
 v[3,2,0]: 2347.83809944685 beta[3,2] + gamma[3,2,0] + delta[3,2,0]
   + epsilon[3,2] <= 0
 v[3,2,1]: 1400.635865576153 beta[3,2] + gamma[3,2,1] + delta[3,2,1]
   + epsilon[3,2] <= 0
 v[3,2,2]: gamma[3,2,2] + delta[3,2,2] + epsilon[3,2] <= 0
 v[3,2,3]: 1026.133235152438 beta[3,2] + gamma[3,2,3] + delta[3,2,3]
   + epsilon[3,2] <= 0
 v[3,2,4]: 1974.997064105249 beta[3,2] + gamma[3,2,4] + delta[3,2,4]
   + epsilon[3,2] <= 0
 v[3,2,5]: 1533.186551957861 beta[3,2] + gamma[3,2,5] + delta[3,2,5]
   + epsilon[3,2] <= 0
 v[3,2,6]: 2481.75905551514 beta[3,2] + gamma[3,2,6] + delta[3,2,6]
   + epsilon[3,2] <= 0
 v[3,2,7]: 1178.440049811752 beta[3,2] + gamma[3,2,7] + delta[3,2,7]
   + epsilon[3,2] <= 0
 v[3,2,8]: 336.9178324846126 beta[3,2] + gamma[3,2,8] + delta[3,2,8]
   + epsilon[3,2] <= 0
 v[3,2,9]: 2869.502213018426 beta[3,2] + gamma[3,2,9] + delta[3,2,9]
   + epsilon[3,2] <= 0
 v[3,3,0]: 372.0097943505882 beta[3,3] + gamma[3,3,0] + delta[3,3,0]
   + epsilon[3,3] <= 0
 v[3,3,1]: 1232.609182263739 beta[3,3] + gamma[3,3,1] + delta[3,3,1]
   + epsilon[3,3] <= 0
 v[3,3,2]: 1780.716124665969 beta[3,3] + gamma[3,3,2] + delta[3,3,2]
   + epsilon[3,3] <= 0
 v[3,3,3]: gamma[3,3,3] + delta[3,3,3] + epsilon[3,3] <= 0
 v[3,3,4]: 562.1895968135387 beta[3,3] + gamma[3,3,4] + delta[3,3,4]
   + epsilon[3,3] <= 0
 v[3,3,5]: 494.2149591513437 beta[3,3] + gamma[3,3,5] + delta[3,3,5]
   + epsilon[3,3] <= 0
 v[3,3,6]: 4320.544174830559 beta[3,3] + gamma[3,3,6] + delta[3,3,6]
   + epsilon[3,3] <= 0
 v[3,3,7]: 2550.24966692996 beta[3,3] + gamma[3,3,7] + delta[3,3,7]
   + epsilon[3,3] <= 0
 v[3,3,8]: 1169.583522913847 beta[3,3] + gamma[3,3,8] + delta[3,3,8]
   + epsilon[3,3] <= 0
 v[3,3,9]: 1843.997505321707 beta[3,3] + gamma[3,3,9] + delta[3,3,9]
   + epsilon[3,3] <= 0
 v[3,4,0]: 135.3543388516725 beta[3,4] + gamma[3,4,0] + delta[3,4,0]
   + epsilon[3,4] <= 0
 v[3,4,1]: 1561.857769164875 beta[3,4] + gamma[3,4,1] + delta[3,4,1]
   + epsilon[3,4] <= 0
 v[3,4,2]: 1399.098603608482 beta[3,4] + gamma[3,4,2] + delta[3,4,2]
   + epsilon[3,4] <= 0
 v[3,4,3]: 615.1863935307424 beta[3,4] + gamma[3,4,3] + delta[3,4,3]
   + epsilon[3,4] <= 0
 v[3,4,4]: gamma[3,4,4] + delta[3,4,4] + epsilon[3,4] <= 0
 v[3,4,5]: 95.00900764974996 beta[3,4] + gamma[3,4,5] + delta[3,4,5]
   + epsilon[3,4] <= 0
 v[3,4,6]: 5775.975966486836 beta[3,4] + gamma[3,4,6] + delta[3,4,6]
   + epsilon[3,4] <= 0
 v[3,4,7]: 4323.54086907186 beta[3,4] + gamma[3,4,7] + delta[3,4,7]
   + epsilon[3,4] <= 0
 v[3,4,8]: 2340.305916370557 beta[3,4] + gamma[3,4,8] + delta[3,4,8]
   + epsilon[3,4] <= 0
 v[3,4,9]: 2928.420658919549 beta[3,4] + gamma[3,4,9] + delta[3,4,9]
   + epsilon[3,4] <= 0
 v[3,5,0]: 35.89871725990424 beta[3,5] + gamma[3,5,0] + delta[3,5,0]
   + epsilon[3,5] <= 0
 v[3,5,1]: 2544.911949467479 beta[3,5] + gamma[3,5,1] + delta[3,5,1]
   + epsilon[3,5] <= 0
 v[3,5,2]: 2561.348697561447 beta[3,5] + gamma[3,5,2] + delta[3,5,2]
   + epsilon[3,5] <= 0
 v[3,5,3]: 309.4552839372835 beta[3,5] + gamma[3,5,3] + delta[3,5,3]
   + epsilon[3,5] <= 0
 v[3,5,4]: 63.98080942306058 beta[3,5] + gamma[3,5,4] + delta[3,5,4]
   + epsilon[3,5] <= 0
 v[3,5,5]: gamma[3,5,5] + delta[3,5,5] + epsilon[3,5] <= 0
 v[3,5,6]: 4192.935674092563 beta[3,5] + gamma[3,5,6] + delta[3,5,6]
   + epsilon[3,5] <= 0
 v[3,5,7]: 2772.621988404359 beta[3,5] + gamma[3,5,7] + delta[3,5,7]
   + epsilon[3,5] <= 0
 v[3,5,8]: 886.4630131418273 beta[3,5] + gamma[3,5,8] + delta[3,5,8]
   + epsilon[3,5] <= 0
 v[3,5,9]: 3392.411231952558 beta[3,5] + gamma[3,5,9] + delta[3,5,9]
   + epsilon[3,5] <= 0
 v[3,6,0]: 3695.811456736657 beta[3,6] + gamma[3,6,0] + delta[3,6,0]
   + epsilon[3,6] <= 0
 v[3,6,1]: 2086.076355466993 beta[3,6] + gamma[3,6,1] + delta[3,6,1]
   + epsilon[3,6] <= 0
 v[3,6,2]: 2222.746343172719 beta[3,6] + gamma[3,6,2] + delta[3,6,2]
   + epsilon[3,6] <= 0
 v[3,6,3]: 2370.411214549593 beta[3,6] + gamma[3,6,3] + delta[3,6,3]
   + epsilon[3,6] <= 0
 v[3,6,4]: 5558.82160591535 beta[3,6] + gamma[3,6,4] + delta[3,6,4]
   + epsilon[3,6] <= 0
 v[3,6,5]: 6925.007839256799 beta[3,6] + gamma[3,6,5] + delta[3,6,5]
   + epsilon[3,6] <= 0
 v[3,6,6]: gamma[3,6,6] + delta[3,6,6] + epsilon[3,6] <= 0
 v[3,6,7]: 2343.347998971427 beta[3,6] + gamma[3,6,7] + delta[3,6,7]
   + epsilon[3,6] <= 0
 v[3,6,8]: 2830.071788585302 beta[3,6] + gamma[3,6,8] + delta[3,6,8]
   + epsilon[3,6] <= 0
 v[3,6,9]: 16.95409328193277 beta[3,6] + gamma[3,6,9] + delta[3,6,9]
   + epsilon[3,6] <= 0
 v[3,7,0]: 2487.212326766596 beta[3,7] + gamma[3,7,0] + delta[3,7,0]
   + epsilon[3,7] <= 0
 v[3,7,1]: 1014.762495143273 beta[3,7] + gamma[3,7,1] + delta[3,7,1]
   + epsilon[3,7] <= 0
 v[3,7,2]: 365.4737072833987 beta[3,7] + gamma[3,7,2] + delta[3,7,2]
   + epsilon[3,7] <= 0
 v[3,7,3]: 1733.507770752834 beta[3,7] + gamma[3,7,3] + delta[3,7,3]
   + epsilon[3,7] <= 0
 v[3,7,4]: 1587.578880190179 beta[3,7] + gamma[3,7,4] + delta[3,7,4]
   + epsilon[3,7] <= 0
 v[3,7,5]: 2438.85187235103 beta[3,7] + gamma[3,7,5] + delta[3,7,5]
   + epsilon[3,7] <= 0
 v[3,7,6]: 511.9290997356711 beta[3,7] + gamma[3,7,6] + delta[3,7,6]
   + epsilon[3,7] <= 0
 v[3,7,7]: gamma[3,7,7] + delta[3,7,7] + epsilon[3,7] <= 0
 v[3,7,8]: 457.386554656099 beta[3,7] + gamma[3,7,8] + delta[3,7,8]
   + epsilon[3,7] <= 0
 v[3,7,9]: 663.5495417234112 beta[3,7] + gamma[3,7,9] + delta[3,7,9]
   + epsilon[3,7] <= 0
 v[3,8,0]: 1473.179303492142 beta[3,8] + gamma[3,8,0] + delta[3,8,0]
   + epsilon[3,8] <= 0
 v[3,8,1]: 1632.026039568535 beta[3,8] + gamma[3,8,1] + delta[3,8,1]
   + epsilon[3,8] <= 0
 v[3,8,2]: 483.8767948552897 beta[3,8] + gamma[3,8,2] + delta[3,8,2]
   + epsilon[3,8] <= 0
 v[3,8,3]: 2056.615051977302 beta[3,8] + gamma[3,8,3] + delta[3,8,3]
   + epsilon[3,8] <= 0
 v[3,8,4]: 2896.46233947692 beta[3,8] + gamma[3,8,4] + delta[3,8,4]
   + epsilon[3,8] <= 0
 v[3,8,5]: 2000.127980524926 beta[3,8] + gamma[3,8,5] + delta[3,8,5]
   + epsilon[3,8] <= 0
 v[3,8,6]: 3656.869537795579 beta[3,8] + gamma[3,8,6] + delta[3,8,6]
   + epsilon[3,8] <= 0
 v[3,8,7]: 1873.913595973519 beta[3,8] + gamma[3,8,7] + delta[3,8,7]
   + epsilon[3,8] <= 0
 v[3,8,8]: gamma[3,8,8] + delta[3,8,8] + epsilon[3,8] <= 0
 v[3,8,9]: 3078.31849212894 beta[3,8] + gamma[3,8,9] + delta[3,8,9]
   + epsilon[3,8] <= 0
 v[3,9,0]: 1610.006973847658 beta[3,9] + gamma[3,9,0] + delta[3,9,0]
   + epsilon[3,9] <= 0
 v[3,9,1]: 450.2691728477851 beta[3,9] + gamma[3,9,1] + delta[3,9,1]
   + epsilon[3,9] <= 0
 v[3,9,2]: 698.551373446965 beta[3,9] + gamma[3,9,2] + delta[3,9,2]
   + epsilon[3,9] <= 0
 v[3,9,3]: 1227.949670236101 beta[3,9] + gamma[3,9,3] + delta[3,9,3]
   + epsilon[3,9] <= 0
 v[3,9,4]: 2349.305833181948 beta[3,9] + gamma[3,9,4] + delta[3,9,4]
   + epsilon[3,9] <= 0
 v[3,9,5]: 1523.582072286289 beta[3,9] + gamma[3,9,5] + delta[3,9,5]
   + epsilon[3,9] <= 0
 v[3,9,6]: 4.131549486155818 beta[3,9] + gamma[3,9,6] + delta[3,9,6]
   + epsilon[3,9] <= 0
 v[3,9,7]: 1138.693964869174 beta[3,9] + gamma[3,9,7] + delta[3,9,7]
   + epsilon[3,9] <= 0
 v[3,9,8]: 1210.656218726193 beta[3,9] + gamma[3,9,8] + delta[3,9,8]
   + epsilon[3,9] <= 0
 v[3,9,9]: gamma[3,9,9] + delta[3,9,9] + epsilon[3,9] <= 0
 v[4,0,0]: gamma[4,0,0] + delta[4,0,0] + epsilon[4,0] <= 0
 v[4,0,1]: 1942.419043769738 beta[4,0] + gamma[4,0,1] + delta[4,0,1]
   + epsilon[4,0] <= 0
 v[4,0,2]: 2916.851789115131 beta[4,0] + gamma[4,0,2] + delta[4,0,2]
   + epsilon[4,0] <= 0
 v[4,0,3]: 378.9588797397554 beta[4,0] + gamma[4,0,3] + delta[4,0,3]
   + epsilon[4,0] <= 0
 v[4,0,4]: 287.1471587661869 beta[4,0] + gamma[4,0,4] + delta[4,0,4]
   + epsilon[4,0] <= 0
 v[4,0,5]: 85.85689957038079 beta[4,0] + gamma[4,0,5] + delta[4,0,5]
   + epsilon[4,0] <= 0
 v[4,0,6]: 7018.207775271378 beta[4,0] + gamma[4,0,6] + delta[4,0,6]
   + epsilon[4,0] <= 0
 v[4,0,7]: 6648.277691402679 beta[4,0] + gamma[4,0,7] + delta[4,0,7]
   + epsilon[4,0] <= 0
 v[4,0,8]: 1678.076441221633 beta[4,0] + gamma[4,0,8] + delta[4,0,8]
   + epsilon[4,0] <= 0
 v[4,0,9]: 3572.265233603591 beta[4,0] + gamma[4,0,9] + delta[4,0,9]
   + epsilon[4,0] <= 0
 v[4,1,0]: 2128.53216941664 beta[4,1] + gamma[4,1,0] + delta[4,1,0]
   + epsilon[4,1] <= 0
 v[4,1,1]: gamma[4,1,1] + delta[4,1,1] + epsilon[4,1] <= 0
 v[4,1,2]: 1975.468937456705 beta[4,1] + gamma[4,1,2] + delta[4,1,2]
   + epsilon[4,1] <= 0
 v[4,1,3]: 3764.647982613453 beta[4,1] + gamma[4,1,3] + delta[4,1,3]
   + epsilon[4,1] <= 0
 v[4,1,4]: 3657.262505363413 beta[4,1] + gamma[4,1,4] + delta[4,1,4]
   + epsilon[4,1] <= 0
 v[4,1,5]: 3778.501816337739 beta[4,1] + gamma[4,1,5] + delta[4,1,5]
   + epsilon[4,1] <= 0
 v[4,1,6]: 1920.855073701297 beta[4,1] + gamma[4,1,6] + delta[4,1,6]
   + epsilon[4,1] <= 0
 v[4,1,7]: 1571.685539211679 beta[4,1] + gamma[4,1,7] + delta[4,1,7]
   + epsilon[4,1] <= 0
 v[4,1,8]: 943.3732290082844 beta[4,1] + gamma[4,1,8] + delta[4,1,8]
   + epsilon[4,1] <= 0
 v[4,1,9]: 925.4924203401357 beta[4,1] + gamma[4,1,9] + delta[4,1,9]
   + epsilon[4,1] <= 0
 v[4,2,0]: 1133.913995214109 beta[4,2] + gamma[4,2,0] + delta[4,2,0]
   + epsilon[4,2] <= 0
 v[4,2,1]: 1478.19674522727 beta[4,2] + gamma[4,2,1] + delta[4,2,1]
   + epsilon[4,2] <= 0
 v[4,2,2]: gamma[4,2,2] + delta[4,2,2] + epsilon[4,2] <= 0
 v[4,2,3]: 1025.927497057943 beta[4,2] + gamma[4,2,3] + delta[4,2,3]
   + epsilon[4,2] <= 0
 v[4,2,4]: 2278.578282078485 beta[4,2] + gamma[4,2,4] + delta[4,2,4]
   + epsilon[4,2] <= 0
 v[4,2,5]: 1456.033365763043 beta[4,2] + gamma[4,2,5] + delta[4,2,5]
   + epsilon[4,2] <= 0
 v[4,2,6]: 2810.049887158864 beta[4,2] + gamma[4,2,6] + delta[4,2,6]
   + epsilon[4,2] <= 0
 v[4,2,7]: 945.6405589067092 beta[4,2] + gamma[4,2,7] + delta[4,2,7]
   + epsilon[4,2] <= 0
 v[4,2,8]: 268.870533356824 beta[4,2] + gamma[4,2,8] + delta[4,2,8]
   + epsilon[4,2] <= 0
 v[4,2,9]: 1552.896192967944 beta[4,2] + gamma[4,2,9] + delta[4,2,9]
   + epsilon[4,2] <= 0
 v[4,3,0]: 644.6464945662184 beta[4,3] + gamma[4,3,0] + delta[4,3,0]
   + epsilon[4,3] <= 0
 v[4,3,1]: 1632.835338742609 beta[4,3] + gamma[4,3,1] + delta[4,3,1]
   + epsilon[4,3] <= 0
 v[4,3,2]: 1166.146663384412 beta[4,3] + gamma[4,3,2] + delta[4,3,2]
   + epsilon[4,3] <= 0
 v[4,3,3]: gamma[4,3,3] + delta[4,3,3] + epsilon[4,3] <= 0
 v[4,3,4]: 539.1737153718011 beta[4,3] + gamma[4,3,4] + delta[4,3,4]
   + epsilon[4,3] <= 0
 v[4,3,5]: 226.9312532470406 beta[4,3] + gamma[4,3,5] + delta[4,3,5]
   + epsilon[4,3] <= 0
 v[4,3,6]: 3459.356803186715 beta[4,3] + gamma[4,3,6] + delta[4,3,6]
   + epsilon[4,3] <= 0
 v[4,3,7]: 3333.416737401949 beta[4,3] + gamma[4,3,7] + delta[4,3,7]
   + epsilon[4,3] <= 0
 v[4,3,8]: 590.5896320785747 beta[4,3] + gamma[4,3,8] + delta[4,3,8]
   + epsilon[4,3] <= 0
 v[4,3,9]: 2915.423602922974 beta[4,3] + gamma[4,3,9] + delta[4,3,9]
   + epsilon[4,3] <= 0
 v[4,4,0]: 63.93150125661513 beta[4,4] + gamma[4,4,0] + delta[4,4,0]
   + epsilon[4,4] <= 0
 v[4,4,1]: 1974.437094425586 beta[4,4] + gamma[4,4,1] + delta[4,4,1]
   + epsilon[4,4] <= 0
 v[4,4,2]: 1097.322460050981 beta[4,4] + gamma[4,4,2] + delta[4,4,2]
   + epsilon[4,4] <= 0
 v[4,4,3]: 287.6321686632712 beta[4,4] + gamma[4,4,3] + delta[4,4,3]
   + epsilon[4,4] <= 0
 v[4,4,4]: gamma[4,4,4] + delta[4,4,4] + epsilon[4,4] <= 0
 v[4,4,5]: 60.53337222469835 beta[4,4] + gamma[4,4,5] + delta[4,4,5]
   + epsilon[4,4] <= 0
 v[4,4,6]: 1639.118658562001 beta[4,4] + gamma[4,4,6] + delta[4,4,6]
   + epsilon[4,4] <= 0
 v[4,4,7]: 903.4908275932709 beta[4,4] + gamma[4,4,7] + delta[4,4,7]
   + epsilon[4,4] <= 0
 v[4,4,8]: 1035.360109310845 beta[4,4] + gamma[4,4,8] + delta[4,4,8]
   + epsilon[4,4] <= 0
 v[4,4,9]: 1688.472877900076 beta[4,4] + gamma[4,4,9] + delta[4,4,9]
   + epsilon[4,4] <= 0
 v[4,5,0]: 21.95018630102432 beta[4,5] + gamma[4,5,0] + delta[4,5,0]
   + epsilon[4,5] <= 0
 v[4,5,1]: 1886.834053431536 beta[4,5] + gamma[4,5,1] + delta[4,5,1]
   + epsilon[4,5] <= 0
 v[4,5,2]: 901.8187472921205 beta[4,5] + gamma[4,5,2] + delta[4,5,2]
   + epsilon[4,5] <= 0
 v[4,5,3]: 186.0221784916917 beta[4,5] + gamma[4,5,3] + delta[4,5,3]
   + epsilon[4,5] <= 0
 v[4,5,4]: 42.20837987138769 beta[4,5] + gamma[4,5,4] + delta[4,5,4]
   + epsilon[4,5] <= 0
 v[4,5,5]: gamma[4,5,5] + delta[4,5,5] + epsilon[4,5] <= 0
 v[4,5,6]: 2728.325080986924 beta[4,5] + gamma[4,5,6] + delta[4,5,6]
   + epsilon[4,5] <= 0
 v[4,5,7]: 2205.091600058674 beta[4,5] + gamma[4,5,7] + delta[4,5,7]
   + epsilon[4,5] <= 0
 v[4,5,8]: 1279.801999935741 beta[4,5] + gamma[4,5,8] + delta[4,5,8]
   + epsilon[4,5] <= 0
 v[4,5,9]: 1386.207765884536 beta[4,5] + gamma[4,5,9] + delta[4,5,9]
   + epsilon[4,5] <= 0
 v[4,6,0]: 2731.783150019083 beta[4,6] + gamma[4,6,0] + delta[4,6,0]
   + epsilon[4,6] <= 0
 v[4,6,1]: 1046.574559607173 beta[4,6] + gamma[4,6,1] + delta[4,6,1]
   + epsilon[4,6] <= 0
 v[4,6,2]: 2033.63838952564 beta[4,6] + gamma[4,6,2] + delta[4,6,2]
   + epsilon[4,6] <= 0
 v[4,6,3]: 2715.591614861054 beta[4,6] + gamma[4,6,3] + delta[4,6,3]
   + epsilon[4,6] <= 0
 v[4,6,4]: 5041.650418980679 beta[4,6] + gamma[4,6,4] + delta[4,6,4]
   + epsilon[4,6] <= 0
 v[4,6,5]: 4217.155251587935 beta[4,6] + gamma[4,6,5] + delta[4,6,5]
   + epsilon[4,6] <= 0
 v[4,6,6]: gamma[4,6,6] + delta[4,6,6] + epsilon[4,6] <= 0
 v[4,6,7]: 1150.658984055146 beta[4,6] + gamma[4,6,7] + delta[4,6,7]
   + epsilon[4,6] <= 0
 v[4,6,8]: 3013.302765313596 beta[4,6] + gamma[4,6,8] + delta[4,6,8]
   + epsilon[4,6] <= 0
 v[4,6,9]: 10.80144144752907 beta[4,6] + gamma[4,6,9] + delta[4,6,9]
   + epsilon[4,6] <= 0
 v[4,7,0]: 2795.029462821647 beta[4,7] + gamma[4,7,0] + delta[4,7,0]
   + epsilon[4,7] <= 0
 v[4,7,1]: 1008.710184326442 beta[4,7] + gamma[4,7,1] + delta[4,7,1]
   + epsilon[4,7] <= 0
 v[4,7,2]: 583.812659592401 beta[4,7] + gamma[4,7,2] + delta[4,7,2]
   + epsilon[4,7] <= 0
 v[4,7,3]: 1594.931867573106 beta[4,7] + gamma[4,7,3] + delta[4,7,3]
   + epsilon[4,7] <= 0
 v[4,7,4]: 2923.908607119819 beta[4,7] + gamma[4,7,4] + delta[4,7,4]
   + epsilon[4,7] <= 0
 v[4,7,5]: 2565.201419522921 beta[4,7] + gamma[4,7,5] + delta[4,7,5]
   + epsilon[4,7] <= 0
 v[4,7,6]: 1103.977412173892 beta[4,7] + gamma[4,7,6] + delta[4,7,6]
   + epsilon[4,7] <= 0
 v[4,7,7]: gamma[4,7,7] + delta[4,7,7] + epsilon[4,7] <= 0
 v[4,7,8]: 1285.856695083922 beta[4,7] + gamma[4,7,8] + delta[4,7,8]
   + epsilon[4,7] <= 0
 v[4,7,9]: 1632.869013178633 beta[4,7] + gamma[4,7,9] + delta[4,7,9]
   + epsilon[4,7] <= 0
 v[4,8,0]: 3711.725455859405 beta[4,8] + gamma[4,8,0] + delta[4,8,0]
   + epsilon[4,8] <= 0
 v[4,8,1]: 1835.610088938056 beta[4,8] + gamma[4,8,1] + delta[4,8,1]
   + epsilon[4,8] <= 0
 v[4,8,2]: 638.7123106839757 beta[4,8] + gamma[4,8,2] + delta[4,8,2]
   + epsilon[4,8] <= 0
 v[4,8,3]: 2806.488441112138 beta[4,8] + gamma[4,8,3] + delta[4,8,3]
   + epsilon[4,8] <= 0
 v[4,8,4]: 2920.359754924157 beta[4,8] + gamma[4,8,4] + delta[4,8,4]
   + epsilon[4,8] <= 0
 v[4,8,5]: 1472.958668631701 beta[4,8] + gamma[4,8,5] + delta[4,8,5]
   + epsilon[4,8] <= 0
 v[4,8,6]: 1879.355663565471 beta[4,8] + gamma[4,8,6] + delta[4,8,6]
   + epsilon[4,8] <= 0
 v[4,8,7]: 1329.215958057111 beta[4,8] + gamma[4,8,7] + delta[4,8,7]
   + epsilon[4,8] <= 0
 v[4,8,8]: gamma[4,8,8] + delta[4,8,8] + epsilon[4,8] <= 0
 v[4,8,9]: 2787.57892372973 beta[4,8] + gamma[4,8,9] + delta[4,8,9]
   + epsilon[4,8] <= 0
 v[4,9,0]: 4249.858826170344 beta[4,9] + gamma[4,9,0] + delta[4,9,0]
   + epsilon[4,9] <= 0
 v[4,9,1]: 598.6426527560576 beta[4,9] + gamma[4,9,1] + delta[4,9,1]
   + epsilon[4,9] <= 0
 v[4,9,2]: 2284.957869294031 beta[4,9] + gamma[4,9,2] + delta[4,9,2]
   + epsilon[4,9] <= 0
 v[4,9,3]: 3161.519348201367 beta[4,9] + gamma[4,9,3] + delta[4,9,3]
   + epsilon[4,9] <= 0
 v[4,9,4]: 3736.145221799648 beta[4,9] + gamma[4,9,4] + delta[4,9,4]
   + epsilon[4,9] <= 0
 v[4,9,5]: 3238.87446134205 beta[4,9] + gamma[4,9,5] + delta[4,9,5]
   + epsilon[4,9] <= 0
 v[4,9,6]: 10.53224093438459 beta[4,9] + gamma[4,9,6] + delta[4,9,6]
   + epsilon[4,9] <= 0
 v[4,9,7]: 1020.473603658251 beta[4,9] + gamma[4,9,7] + delta[4,9,7]
   + epsilon[4,9] <= 0
 v[4,9,8]: 799.1925699262664 beta[4,9] + gamma[4,9,8] + delta[4,9,8]
   + epsilon[4,9] <= 0
 v[4,9,9]: gamma[4,9,9] + delta[4,9,9] + epsilon[4,9] <= 0
 v[5,0,0]: gamma[5,0,0] + delta[5,0,0] + epsilon[5,0] <= 0
 v[5,0,1]: 3366.568942109262 beta[5,0] + gamma[5,0,1] + delta[5,0,1]
   + epsilon[5,0] <= 0
 v[5,0,2]: 3044.596438350439 beta[5,0] + gamma[5,0,2] + delta[5,0,2]
   + epsilon[5,0] <= 0
 v[5,0,3]: 627.3219655023812 beta[5,0] + gamma[5,0,3] + delta[5,0,3]
   + epsilon[5,0] <= 0
 v[5,0,4]: 106.8082353804894 beta[5,0] + gamma[5,0,4] + delta[5,0,4]
   + epsilon[5,0] <= 0
 v[5,0,5]: 42.88127269300274 beta[5,0] + gamma[5,0,5] + delta[5,0,5]
   + epsilon[5,0] <= 0
 v[5,0,6]: 4699.584505235594 beta[5,0] + gamma[5,0,6] + delta[5,0,6]
   + epsilon[5,0] <= 0
 v[5,0,7]: 4763.758078641688 beta[5,0] + gamma[5,0,7] + delta[5,0,7]
   + epsilon[5,0] <= 0
 v[5,0,8]: 2852.333162725611 beta[5,0] + gamma[5,0,8] + delta[5,0,8]
   + epsilon[5,0] <= 0
 v[5,0,9]: 2957.582695163331 beta[5,0] + gamma[5,0,9] + delta[5,0,9]
   + epsilon[5,0] <= 0
 v[5,1,0]: 800.144603101903 beta[5,1] + gamma[5,1,0] + delta[5,1,0]
   + epsilon[5,1] <= 0
 v[5,1,1]: gamma[5,1,1] + delta[5,1,1] + epsilon[5,1] <= 0
 v[5,1,2]: 909.9640000268997 beta[5,1] + gamma[5,1,2] + delta[5,1,2]
   + epsilon[5,1] <= 0
 v[5,1,3]: 640.2445470946896 beta[5,1] + gamma[5,1,3] + delta[5,1,3]
   + epsilon[5,1] <= 0
 v[5,1,4]: 1750.179864946978 beta[5,1] + gamma[5,1,4] + delta[5,1,4]
   + epsilon[5,1] <= 0
 v[5,1,5]: 1799.586659107858 beta[5,1] + gamma[5,1,5] + delta[5,1,5]
   + epsilon[5,1] <= 0
 v[5,1,6]: 639.565881033687 beta[5,1] + gamma[5,1,6] + delta[5,1,6]
   + epsilon[5,1] <= 0
 v[5,1,7]: 1220.130370205479 beta[5,1] + gamma[5,1,7] + delta[5,1,7]
   + epsilon[5,1] <= 0
 v[5,1,8]: 711.680004877402 beta[5,1] + gamma[5,1,8] + delta[5,1,8]
   + epsilon[5,1] <= 0
 v[5,1,9]: 343.2344011712821 beta[5,1] + gamma[5,1,9] + delta[5,1,9]
   + epsilon[5,1] <= 0
 v[5,2,0]: 2395.303895286227 beta[5,2] + gamma[5,2,0] + delta[5,2,0]
   + epsilon[5,2] <= 0
 v[5,2,1]: 2034.814836264249 beta[5,2] + gamma[5,2,1] + delta[5,2,1]
   + epsilon[5,2] <= 0
 v[5,2,2]: gamma[5,2,2] + delta[5,2,2] + epsilon[5,2] <= 0
 v[5,2,3]: 2050.407111401628 beta[5,2] + gamma[5,2,3] + delta[5,2,3]
   + epsilon[5,2] <= 0
 v[5,2,4]: 2694.969587332299 beta[5,2] + gamma[5,2,4] + delta[5,2,4]
   + epsilon[5,2] <= 0
 v[5,2,5]: 2575.129589844264 beta[5,2] + gamma[5,2,5] + delta[5,2,5]
   + epsilon[5,2] <= 0
 v[5,2,6]: 1950.559598660871 beta[5,2] + gamma[5,2,6] + delta[5,2,6]
   + epsilon[5,2] <= 0
 v[5,2,7]: 1738.864961723707 beta[5,2] + gamma[5,2,7] + delta[5,2,7]
   + epsilon[5,2] <= 0
 v[5,2,8]: 405.1383126519848 beta[5,2] + gamma[5,2,8] + delta[5,2,8]
   + epsilon[5,2] <= 0
 v[5,2,9]: 2160.404210567409 beta[5,2] + gamma[5,2,9] + delta[5,2,9]
   + epsilon[5,2] <= 0
 v[5,3,0]: 392.5224161595936 beta[5,3] + gamma[5,3,0] + delta[5,3,0]
   + epsilon[5,3] <= 0
 v[5,3,1]: 829.7530784666234 beta[5,3] + gamma[5,3,1] + delta[5,3,1]
   + epsilon[5,3] <= 0
 v[5,3,2]: 987.1184690046235 beta[5,3] + gamma[5,3,2] + delta[5,3,2]
   + epsilon[5,3] <= 0
 v[5,3,3]: gamma[5,3,3] + delta[5,3,3] + epsilon[5,3] <= 0
 v[5,3,4]: 302.4115623895616 beta[5,3] + gamma[5,3,4] + delta[5,3,4]
   + epsilon[5,3] <= 0
 v[5,3,5]: 251.2243551613103 beta[5,3] + gamma[5,3,5] + delta[5,3,5]
   + epsilon[5,3] <= 0
 v[5,3,6]: 2355.325241633947 beta[5,3] + gamma[5,3,6] + delta[5,3,6]
   + epsilon[5,3] <= 0
 v[5,3,7]: 1295.904618799056 beta[5,3] + gamma[5,3,7] + delta[5,3,7]
   + epsilon[5,3] <= 0
 v[5,3,8]: 506.4602456960559 beta[5,3] + gamma[5,3,8] + delta[5,3,8]
   + epsilon[5,3] <= 0
 v[5,3,9]: 2290.331404905034 beta[5,3] + gamma[5,3,9] + delta[5,3,9]
   + epsilon[5,3] <= 0
 v[5,4,0]: 145.8274023015942 beta[5,4] + gamma[5,4,0] + delta[5,4,0]
   + epsilon[5,4] <= 0
 v[5,4,1]: 4173.107280177796 beta[5,4] + gamma[5,4,1] + delta[5,4,1]
   + epsilon[5,4] <= 0
 v[5,4,2]: 3005.832726032942 beta[5,4] + gamma[5,4,2] + delta[5,4,2]
   + epsilon[5,4] <= 0
 v[5,4,3]: 406.1852550242829 beta[5,4] + gamma[5,4,3] + delta[5,4,3]
   + epsilon[5,4] <= 0
 v[5,4,4]: gamma[5,4,4] + delta[5,4,4] + epsilon[5,4] <= 0
 v[5,4,5]: 87.36837973605078 beta[5,4] + gamma[5,4,5] + delta[5,4,5]
   + epsilon[5,4] <= 0
 v[5,4,6]: 3105.691053175436 beta[5,4] + gamma[5,4,6] + delta[5,4,6]
   + epsilon[5,4] <= 0
 v[5,4,7]: 3624.864741237798 beta[5,4] + gamma[5,4,7] + delta[5,4,7]
   + epsilon[5,4] <= 0
 v[5,4,8]: 1179.833344904962 beta[5,4] + gamma[5,4,8] + delta[5,4,8]
   + epsilon[5,4] <= 0
 v[5,4,9]: 2961.847457061712 beta[5,4] + gamma[5,4,9] + delta[5,4,9]
   + epsilon[5,4] <= 0
 v[5,5,0]: 33.59273862449271 beta[5,5] + gamma[5,5,0] + delta[5,5,0]
   + epsilon[5,5] <= 0
 v[5,5,1]: 2294.987917628485 beta[5,5] + gamma[5,5,1] + delta[5,5,1]
   + epsilon[5,5] <= 0
 v[5,5,2]: 748.9894234173127 beta[5,5] + gamma[5,5,2] + delta[5,5,2]
   + epsilon[5,5] <= 0
 v[5,5,3]: 241.3871123680144 beta[5,5] + gamma[5,5,3] + delta[5,5,3]
   + epsilon[5,5] <= 0
 v[5,5,4]: 45.43634302848009 beta[5,5] + gamma[5,5,4] + delta[5,5,4]
   + epsilon[5,5] <= 0
 v[5,5,5]: gamma[5,5,5] + delta[5,5,5] + epsilon[5,5] <= 0
 v[5,5,6]: 3048.609279614253 beta[5,5] + gamma[5,5,6] + delta[5,5,6]
   + epsilon[5,5] <= 0
 v[5,5,7]: 2291.287265885019 beta[5,5] + gamma[5,5,7] + delta[5,5,7]
   + epsilon[5,5] <= 0
 v[5,5,8]: 960.9322307882425 beta[5,5] + gamma[5,5,8] + delta[5,5,8]
   + epsilon[5,5] <= 0
 v[5,5,9]: 3369.426097816719 beta[5,5] + gamma[5,5,9] + delta[5,5,9]
   + epsilon[5,5] <= 0
 v[5,6,0]: 3706.699551005164 beta[5,6] + gamma[5,6,0] + delta[5,6,0]
   + epsilon[5,6] <= 0
 v[5,6,1]: 1745.729969414586 beta[5,6] + gamma[5,6,1] + delta[5,6,1]
   + epsilon[5,6] <= 0
 v[5,6,2]: 2566.826899153472 beta[5,6] + gamma[5,6,2] + delta[5,6,2]
   + epsilon[5,6] <= 0
 v[5,6,3]: 3354.64282635005 beta[5,6] + gamma[5,6,3] + delta[5,6,3]
   + epsilon[5,6] <= 0
 v[5,6,4]: 2430.470876307733 beta[5,6] + gamma[5,6,4] + delta[5,6,4]
   + epsilon[5,6] <= 0
 v[5,6,5]: 5806.993106197813 beta[5,6] + gamma[5,6,5] + delta[5,6,5]
   + epsilon[5,6] <= 0
 v[5,6,6]: gamma[5,6,6] + delta[5,6,6] + epsilon[5,6] <= 0
 v[5,6,7]: 1005.692096762047 beta[5,6] + gamma[5,6,7] + delta[5,6,7]
   + epsilon[5,6] <= 0
 v[5,6,8]: 1989.754778629184 beta[5,6] + gamma[5,6,8] + delta[5,6,8]
   + epsilon[5,6] <= 0
 v[5,6,9]: 11.24317410567547 beta[5,6] + gamma[5,6,9] + delta[5,6,9]
   + epsilon[5,6] <= 0
 v[5,7,0]: 1118.189267332832 beta[5,7] + gamma[5,7,0] + delta[5,7,0]
   + epsilon[5,7] <= 0
 v[5,7,1]: 970.5412275802296 beta[5,7] + gamma[5,7,1] + delta[5,7,1]
   + epsilon[5,7] <= 0
 v[5,7,2]: 479.7212271686007 beta[5,7] + gamma[5,7,2] + delta[5,7,2]
   + epsilon[5,7] <= 0
 v[5,7,3]: 1663.398840640628 beta[5,7] + gamma[5,7,3] + delta[5,7,3]
   + epsilon[5,7] <= 0
 v[5,7,4]: 2093.439415321746 beta[5,7] + gamma[5,7,4] + delta[5,7,4]
   + epsilon[5,7] <= 0
 v[5,7,5]: 1214.808601015171 beta[5,7] + gamma[5,7,5] + delta[5,7,5]
   + epsilon[5,7] <= 0
 v[5,7,6]: 792.561233558415 beta[5,7] + gamma[5,7,6] + delta[5,7,6]
   + epsilon[5,7] <= 0
 v[5,7,7]: gamma[5,7,7] + delta[5,7,7] + epsilon[5,7] <= 0
 v[5,7,8]: 908.2981334073946 beta[5,7] + gamma[5,7,8] + delta[5,7,8]
   + epsilon[5,7] <= 0
 v[5,7,9]: 795.3210421156497 beta[5,7] + gamma[5,7,9] + delta[5,7,9]
   + epsilon[5,7] <= 0
 v[5,8,0]: 1944.22554268966 beta[5,8] + gamma[5,8,0] + delta[5,8,0]
   + epsilon[5,8] <= 0
 v[5,8,1]: 1791.260960541458 beta[5,8] + gamma[5,8,1] + delta[5,8,1]
   + epsilon[5,8] <= 0
 v[5,8,2]: 542.3184352763733 beta[5,8] + gamma[5,8,2] + delta[5,8,2]
   + epsilon[5,8] <= 0
 v[5,8,3]: 1607.043658656979 beta[5,8] + gamma[5,8,3] + delta[5,8,3]
   + epsilon[5,8] <= 0
 v[5,8,4]: 2963.483293069975 beta[5,8] + gamma[5,8,4] + delta[5,8,4]
   + epsilon[5,8] <= 0
 v[5,8,5]: 1334.097798755169 beta[5,8] + gamma[5,8,5] + delta[5,8,5]
   + epsilon[5,8] <= 0
 v[5,8,6]: 3665.836899749547 beta[5,8] + gamma[5,8,6] + delta[5,8,6]
   + epsilon[5,8] <= 0
 v[5,8,7]: 1411.22456414484 beta[5,8] + gamma[5,8,7] + delta[5,8,7]
   + epsilon[5,8] <= 0
 v[5,8,8]: gamma[5,8,8] + delta[5,8,8] + epsilon[5,8] <= 0
 v[5,8,9]: 2209.067642659423 beta[5,8] + gamma[5,8,9] + delta[5,8,9]
   + epsilon[5,8] <= 0
 v[5,9,0]: 1599.971582161084 beta[5,9] + gamma[5,9,0] + delta[5,9,0]
   + epsilon[5,9] <= 0
 v[5,9,1]: 1205.289988172213 beta[5,9] + gamma[5,9,1] + delta[5,9,1]
   + epsilon[5,9] <= 0
 v[5,9,2]: 1292.592205599043 beta[5,9] + gamma[5,9,2] + delta[5,9,2]
   + epsilon[5,9] <= 0
 v[5,9,3]: 2461.585509607941 beta[5,9] + gamma[5,9,3] + delta[5,9,3]
   + epsilon[5,9] <= 0
 v[5,9,4]: 2685.481580963258 beta[5,9] + gamma[5,9,4] + delta[5,9,4]
   + epsilon[5,9] <= 0
 v[5,9,5]: 1730.87829832792 beta[5,9] + gamma[5,9,5] + delta[5,9,5]
   + epsilon[5,9] <= 0
 v[5,9,6]: 6.095603808027843 beta[5,9] + gamma[5,9,6] + delta[5,9,6]
   + epsilon[5,9] <= 0
 v[5,9,7]: 576.3039307707015 beta[5,9] + gamma[5,9,7] + delta[5,9,7]
   + epsilon[5,9] <= 0
 v[5,9,8]: 1503.047444355329 beta[5,9] + gamma[5,9,8] + delta[5,9,8]
   + epsilon[5,9] <= 0
 v[5,9,9]: gamma[5,9,9] + delta[5,9,9] + epsilon[5,9] <= 0
 v[6,0,0]: gamma[6,0,0] + delta[6,0,0] + epsilon[6,0] <= 0
 v[6,0,1]: 3280.463764599395 beta[6,0] + gamma[6,0,1] + delta[6,0,1]
   + epsilon[6,0] <= 0
 v[6,0,2]: 1268.662029794709 beta[6,0] + gamma[6,0,2] + delta[6,0,2]
   + epsilon[6,0] <= 0
 v[6,0,3]: 550.3906730108777 beta[6,0] + gamma[6,0,3] + delta[6,0,3]
   + epsilon[6,0] <= 0
 v[6,0,4]: 98.55759018708729 beta[6,0] + gamma[6,0,4] + delta[6,0,4]
   + epsilon[6,0] <= 0
 v[6,0,5]: 53.46772503225299 beta[6,0] + gamma[6,0,5] + delta[6,0,5]
   + epsilon[6,0] <= 0
 v[6,0,6]: 4603.120759598152 beta[6,0] + gamma[6,0,6] + delta[6,0,6]
   + epsilon[6,0] <= 0
 v[6,0,7]: 2208.978723996565 beta[6,0] + gamma[6,0,7] + delta[6,0,7]
   + epsilon[6,0] <= 0
 v[6,0,8]: 2460.744318263748 beta[6,0] + gamma[6,0,8] + delta[6,0,8]
   + epsilon[6,0] <= 0
 v[6,0,9]: 1962.163531537228 beta[6,0] + gamma[6,0,9] + delta[6,0,9]
   + epsilon[6,0] <= 0
 v[6,1,0]: 5014.596055325719 beta[6,1] + gamma[6,1,0] + delta[6,1,0]
   + epsilon[6,1] <= 0
 v[6,1,1]: gamma[6,1,1] + delta[6,1,1] + epsilon[6,1] <= 0
 v[6,1,2]: 1511.591812669658 beta[6,1] + gamma[6,1,2] + delta[6,1,2]
   + epsilon[6,1] <= 0
 v[6,1,3]: 3787.796337679728 beta[6,1] + gamma[6,1,3] + delta[6,1,3]
   + epsilon[6,1] <= 0
 v[6,1,4]: 4376.534365131678 beta[6,1] + gamma[6,1,4] + delta[6,1,4]
   + epsilon[6,1] <= 0
 v[6,1,5]: 5433.535864978067 beta[6,1] + gamma[6,1,5] + delta[6,1,5]
   + epsilon[6,1] <= 0
 v[6,1,6]: 1263.040676074112 beta[6,1] + gamma[6,1,6] + delta[6,1,6]
   + epsilon[6,1] <= 0
 v[6,1,7]: 1046.770482482137 beta[6,1] + gamma[6,1,7] + delta[6,1,7]
   + epsilon[6,1] <= 0
 v[6,1,8]: 2148.450927290571 beta[6,1] + gamma[6,1,8] + delta[6,1,8]
   + epsilon[6,1] <= 0
 v[6,1,9]: 2222.512359162442 beta[6,1] + gamma[6,1,9] + delta[6,1,9]
   + epsilon[6,1] <= 0
 v[6,2,0]: 1844.375684889682 beta[6,2] + gamma[6,2,0] + delta[6,2,0]
   + epsilon[6,2] <= 0
 v[6,2,1]: 938.7839573293779 beta[6,2] + gamma[6,2,1] + delta[6,2,1]
   + epsilon[6,2] <= 0
 v[6,2,2]: gamma[6,2,2] + delta[6,2,2] + epsilon[6,2] <= 0
 v[6,2,3]: 1409.370369953523 beta[6,2] + gamma[6,2,3] + delta[6,2,3]
   + epsilon[6,2] <= 0
 v[6,2,4]: 1239.131963650771 beta[6,2] + gamma[6,2,4] + delta[6,2,4]
   + epsilon[6,2] <= 0
 v[6,2,5]: 1870.970094466238 beta[6,2] + gamma[6,2,5] + delta[6,2,5]
   + epsilon[6,2] <= 0
 v[6,2,6]: 1943.800413119796 beta[6,2] + gamma[6,2,6] + delta[6,2,6]
   + epsilon[6,2] <= 0
 v[6,2,7]: 1125.885462489996 beta[6,2] + gamma[6,2,7] + delta[6,2,7]
   + epsilon[6,2] <= 0
 v[6,2,8]: 192.3497211409067 beta[6,2] + gamma[6,2,8] + delta[6,2,8]
   + epsilon[6,2] <= 0
 v[6,2,9]: 1871.78562162801 beta[6,2] + gamma[6,2,9] + delta[6,2,9]
   + epsilon[6,2] <= 0
 v[6,3,0]: 144.9046994270823 beta[6,3] + gamma[6,3,0] + delta[6,3,0]
   + epsilon[6,3] <= 0
 v[6,3,1]: 1514.164045177136 beta[6,3] + gamma[6,3,1] + delta[6,3,1]
   + epsilon[6,3] <= 0
 v[6,3,2]: 1099.419459563269 beta[6,3] + gamma[6,3,2] + delta[6,3,2]
   + epsilon[6,3] <= 0
 v[6,3,3]: gamma[6,3,3] + delta[6,3,3] + epsilon[6,3] <= 0
 v[6,3,4]: 146.6167024708515 beta[6,3] + gamma[6,3,4] + delta[6,3,4]
   + epsilon[6,3] <= 0
 v[6,3,5]: 194.0429084022978 beta[6,3] + gamma[6,3,5] + delta[6,3,5]
   + epsilon[6,3] <= 0
 v[6,3,6]: 1437.233541517561 beta[6,3] + gamma[6,3,6] + delta[6,3,6]
   + epsilon[6,3] <= 0
 v[6,3,7]: 2089.008524045307 beta[6,3] + gamma[6,3,7] + delta[6,3,7]
   + epsilon[6,3] <= 0
 v[6,3,8]: 736.4524119547174 beta[6,3] + gamma[6,3,8] + delta[6,3,8]
   + epsilon[6,3] <= 0
 v[6,3,9]: 1414.480908067322 beta[6,3] + gamma[6,3,9] + delta[6,3,9]
   + epsilon[6,3] <= 0
 v[6,4,0]: 49.92404153558132 beta[6,4] + gamma[6,4,0] + delta[6,4,0]
   + epsilon[6,4] <= 0
 v[6,4,1]: 1159.528620143912 beta[6,4] + gamma[6,4,1] + delta[6,4,1]
   + epsilon[6,4] <= 0
 v[6,4,2]: 1166.753245844569 beta[6,4] + gamma[6,4,2] + delta[6,4,2]
   + epsilon[6,4] <= 0
 v[6,4,3]: 201.2403342160736 beta[6,4] + gamma[6,4,3] + delta[6,4,3]
   + epsilon[6,4] <= 0
 v[6,4,4]: gamma[6,4,4] + delta[6,4,4] + epsilon[6,4] <= 0
 v[6,4,5]: 51.17998318869012 beta[6,4] + gamma[6,4,5] + delta[6,4,5]
   + epsilon[6,4] <= 0
 v[6,4,6]: 2784.068218396464 beta[6,4] + gamma[6,4,6] + delta[6,4,6]
   + epsilon[6,4] <= 0
 v[6,4,7]: 1660.314499542357 beta[6,4] + gamma[6,4,7] + delta[6,4,7]
   + epsilon[6,4] <= 0
 v[6,4,8]: 821.66984640082 beta[6,4] + gamma[6,4,8] + delta[6,4,8]
   + epsilon[6,4] <= 0
 v[6,4,9]: 2470.416972931881 beta[6,4] + gamma[6,4,9] + delta[6,4,9]
   + epsilon[6,4] <= 0
 v[6,5,0]: 15.66546539077708 beta[6,5] + gamma[6,5,0] + delta[6,5,0]
   + epsilon[6,5] <= 0
 v[6,5,1]: 1296.574944580181 beta[6,5] + gamma[6,5,1] + delta[6,5,1]
   + epsilon[6,5] <= 0
 v[6,5,2]: 798.0556070114657 beta[6,5] + gamma[6,5,2] + delta[6,5,2]
   + epsilon[6,5] <= 0
 v[6,5,3]: 282.2259686283834 beta[6,5] + gamma[6,5,3] + delta[6,5,3]
   + epsilon[6,5] <= 0
 v[6,5,4]: 67.14233843728211 beta[6,5] + gamma[6,5,4] + delta[6,5,4]
   + epsilon[6,5] <= 0
 v[6,5,5]: gamma[6,5,5] + delta[6,5,5] + epsilon[6,5] <= 0
 v[6,5,6]: 2454.220657771137 beta[6,5] + gamma[6,5,6] + delta[6,5,6]
   + epsilon[6,5] <= 0
 v[6,5,7]: 1579.544631261194 beta[6,5] + gamma[6,5,7] + delta[6,5,7]
   + epsilon[6,5] <= 0
 v[6,5,8]: 1048.455579577194 beta[6,5] + gamma[6,5,8] + delta[6,5,8]
   + epsilon[6,5] <= 0
 v[6,5,9]: 1210.818339433429 beta[6,5] + gamma[6,5,9] + delta[6,5,9]
   + epsilon[6,5] <= 0
 v[6,6,0]: 4192.393370810643 beta[6,6] + gamma[6,6,0] + delta[6,6,0]
   + epsilon[6,6] <= 0
 v[6,6,1]: 1959.726233995054 beta[6,6] + gamma[6,6,1] + delta[6,6,1]
   + epsilon[6,6] <= 0
 v[6,6,2]: 2667.696517673884 beta[6,6] + gamma[6,6,2] + delta[6,6,2]
   + epsilon[6,6] <= 0
 v[6,6,3]: 4926.609925383583 beta[6,6] + gamma[6,6,3] + delta[6,6,3]
   + epsilon[6,6] <= 0
 v[6,6,4]: 4421.425598233825 beta[6,6] + gamma[6,6,4] + delta[6,6,4]
   + epsilon[6,6] <= 0
 v[6,6,5]: 4452.072119719795 beta[6,6] + gamma[6,6,5] + delta[6,6,5]
   + epsilon[6,6] <= 0
 v[6,6,6]: gamma[6,6,6] + delta[6,6,6] + epsilon[6,6] <= 0
 v[6,6,7]: 921.5344506484532 beta[6,6] + gamma[6,6,7] + delta[6,6,7]
   + epsilon[6,6] <= 0
 v[6,6,8]: 2600.341055466342 beta[6,6] + gamma[6,6,8] + delta[6,6,8]
   + epsilon[6,6] <= 0
 v[6,6,9]: 20.34106519031348 beta[6,6] + gamma[6,6,9] + delta[6,6,9]
   + epsilon[6,6] <= 0
 v[6,7,0]: 3438.251966369617 beta[6,7] + gamma[6,7,0] + delta[6,7,0]
   + epsilon[6,7] <= 0
 v[6,7,1]: 1071.203630877096 beta[6,7] + gamma[6,7,1] + delta[6,7,1]
   + epsilon[6,7] <= 0
 v[6,7,2]: 1834.309199399208 beta[6,7] + gamma[6,7,2] + delta[6,7,2]
   + epsilon[6,7] <= 0
 v[6,7,3]: 2135.458424354362 beta[6,7] + gamma[6,7,3] + delta[6,7,3]
   + epsilon[6,7] <= 0
 v[6,7,4]: 3915.367991922476 beta[6,7] + gamma[6,7,4] + delta[6,7,4]
   + epsilon[6,7] <= 0
 v[6,7,5]: 3473.887779895068 beta[6,7] + gamma[6,7,5] + delta[6,7,5]
   + epsilon[6,7] <= 0
 v[6,7,6]: 2037.767415092303 beta[6,7] + gamma[6,7,6] + delta[6,7,6]
   + epsilon[6,7] <= 0
 v[6,7,7]: gamma[6,7,7] + delta[6,7,7] + epsilon[6,7] <= 0
 v[6,7,8]: 2086.043715794987 beta[6,7] + gamma[6,7,8] + delta[6,7,8]
   + epsilon[6,7] <= 0
 v[6,7,9]: 719.3885440320067 beta[6,7] + gamma[6,7,9] + delta[6,7,9]
   + epsilon[6,7] <= 0
 v[6,8,0]: 3661.175564309758 beta[6,8] + gamma[6,8,0] + delta[6,8,0]
   + epsilon[6,8] <= 0
 v[6,8,1]: 845.1051980476722 beta[6,8] + gamma[6,8,1] + delta[6,8,1]
   + epsilon[6,8] <= 0
 v[6,8,2]: 310.2734549898541 beta[6,8] + gamma[6,8,2] + delta[6,8,2]
   + epsilon[6,8] <= 0
 v[6,8,3]: 2525.498183092289 beta[6,8] + gamma[6,8,3] + delta[6,8,3]
   + epsilon[6,8] <= 0
 v[6,8,4]: 1183.12121642995 beta[6,8] + gamma[6,8,4] + delta[6,8,4]
   + epsilon[6,8] <= 0
 v[6,8,5]: 3191.42381860685 beta[6,8] + gamma[6,8,5] + delta[6,8,5]
   + epsilon[6,8] <= 0
 v[6,8,6]: 3758.717904200735 beta[6,8] + gamma[6,8,6] + delta[6,8,6]
   + epsilon[6,8] <= 0
 v[6,8,7]: 1927.051786870497 beta[6,8] + gamma[6,8,7] + delta[6,8,7]
   + epsilon[6,8] <= 0
 v[6,8,8]: gamma[6,8,8] + delta[6,8,8] + epsilon[6,8] <= 0
 v[6,8,9]: 3617.649645914808 beta[6,8] + gamma[6,8,9] + delta[6,8,9]
   + epsilon[6,8] <= 0
 v[6,9,0]: 3938.648158290816 beta[6,9] + gamma[6,9,0] + delta[6,9,0]
   + epsilon[6,9] <= 0
 v[6,9,1]: 794.9016235253828 beta[6,9] + gamma[6,9,1] + delta[6,9,1]
   + epsilon[6,9] <= 0
 v[6,9,2]: 3097.494488591146 beta[6,9] + gamma[6,9,2] + delta[6,9,2]
   + epsilon[6,9] <= 0
 v[6,9,3]: 4868.809074900084 beta[6,9] + gamma[6,9,3] + delta[6,9,3]
   + epsilon[6,9] <= 0
 v[6,9,4]: 3911.653994578981 beta[6,9] + gamma[6,9,4] + delta[6,9,4]
   + epsilon[6,9] <= 0
 v[6,9,5]: 5424.832080634838 beta[6,9] + gamma[6,9,5] + delta[6,9,5]
   + epsilon[6,9] <= 0
 v[6,9,6]: 11.57475477019147 beta[6,9] + gamma[6,9,6] + delta[6,9,6]
   + epsilon[6,9] <= 0
 v[6,9,7]: 1728.642093930039 beta[6,9] + gamma[6,9,7] + delta[6,9,7]
   + epsilon[6,9] <= 0
 v[6,9,8]: 1955.602767890864 beta[6,9] + gamma[6,9,8] + delta[6,9,8]
   + epsilon[6,9] <= 0
 v[6,9,9]: gamma[6,9,9] + delta[6,9,9] + epsilon[6,9] <= 0
 v[7,0,0]: gamma[7,0,0] + delta[7,0,0] + epsilon[7,0] <= 0
 v[7,0,1]: 2536.798454202475 beta[7,0] + gamma[7,0,1] + delta[7,0,1]
   + epsilon[7,0] <= 0
 v[7,0,2]: 2000.376461196805 beta[7,0] + gamma[7,0,2] + delta[7,0,2]
   + epsilon[7,0] <= 0
 v[7,0,3]: 380.6676229486025 beta[7,0] + gamma[7,0,3] + delta[7,0,3]
   + epsilon[7,0] <= 0
 v[7,0,4]: 126.5936617582461 beta[7,0] + gamma[7,0,4] + delta[7,0,4]
   + epsilon[7,0] <= 0
 v[7,0,5]: 23.14887452509511 beta[7,0] + gamma[7,0,5] + delta[7,0,5]
   + epsilon[7,0] <= 0
 v[7,0,6]: 2480.771275593424 beta[7,0] + gamma[7,0,6] + delta[7,0,6]
   + epsilon[7,0] <= 0
 v[7,0,7]: 1761.189593984904 beta[7,0] + gamma[7,0,7] + delta[7,0,7]
   + epsilon[7,0] <= 0
 v[7,0,8]: 1316.571902164783 beta[7,0] + gamma[7,0,8] + delta[7,0,8]
   + epsilon[7,0] <= 0
 v[7,0,9]: 2244.591444744732 beta[7,0] + gamma[7,0,9] + delta[7,0,9]
   + epsilon[7,0] <= 0
 v[7,1,0]: 1596.400731572047 beta[7,1] + gamma[7,1,0] + delta[7,1,0]
   + epsilon[7,1] <= 0
 v[7,1,1]: gamma[7,1,1] + delta[7,1,1] + epsilon[7,1] <= 0
 v[7,1,2]: 655.5024418258821 beta[7,1] + gamma[7,1,2] + delta[7,1,2]
   + epsilon[7,1] <= 0
 v[7,1,3]: 2130.055192871159 beta[7,1] + gamma[7,1,3] + delta[7,1,3]
   + epsilon[7,1] <= 0
 v[7,1,4]: 1951.241669673628 beta[7,1] + gamma[7,1,4] + delta[7,1,4]
   + epsilon[7,1] <= 0
 v[7,1,5]: 1526.441468180149 beta[7,1] + gamma[7,1,5] + delta[7,1,5]
   + epsilon[7,1] <= 0
 v[7,1,6]: 1052.1453944205 beta[7,1] + gamma[7,1,6] + delta[7,1,6]
   + epsilon[7,1] <= 0
 v[7,1,7]: 1110.726458739631 beta[7,1] + gamma[7,1,7] + delta[7,1,7]
   + epsilon[7,1] <= 0
 v[7,1,8]: 538.0137072814733 beta[7,1] + gamma[7,1,8] + delta[7,1,8]
   + epsilon[7,1] <= 0
 v[7,1,9]: 1010.36427515229 beta[7,1] + gamma[7,1,9] + delta[7,1,9]
   + epsilon[7,1] <= 0
 v[7,2,0]: 1626.859579761239 beta[7,2] + gamma[7,2,0] + delta[7,2,0]
   + epsilon[7,2] <= 0
 v[7,2,1]: 919.9500705154006 beta[7,2] + gamma[7,2,1] + delta[7,2,1]
   + epsilon[7,2] <= 0
 v[7,2,2]: gamma[7,2,2] + delta[7,2,2] + epsilon[7,2] <= 0
 v[7,2,3]: 1531.123914027067 beta[7,2] + gamma[7,2,3] + delta[7,2,3]
   + epsilon[7,2] <= 0
 v[7,2,4]: 1264.87215343127 beta[7,2] + gamma[7,2,4] + delta[7,2,4]
   + epsilon[7,2] <= 0
 v[7,2,5]: 2238.7197712198 beta[7,2] + gamma[7,2,5] + delta[7,2,5]
   + epsilon[7,2] <= 0
 v[7,2,6]: 1436.13910844541 beta[7,2] + gamma[7,2,6] + delta[7,2,6]
   + epsilon[7,2] <= 0
 v[7,2,7]: 1226.91499748658 beta[7,2] + gamma[7,2,7] + delta[7,2,7]
   + epsilon[7,2] <= 0
 v[7,2,8]: 244.8041983702453 beta[7,2] + gamma[7,2,8] + delta[7,2,8]
   + epsilon[7,2] <= 0
 v[7,2,9]: 2461.787575718464 beta[7,2] + gamma[7,2,9] + delta[7,2,9]
   + epsilon[7,2] <= 0
 v[7,3,0]: 402.5611630802731 beta[7,3] + gamma[7,3,0] + delta[7,3,0]
   + epsilon[7,3] <= 0
 v[7,3,1]: 818.2811328301515 beta[7,3] + gamma[7,3,1] + delta[7,3,1]
   + epsilon[7,3] <= 0
 v[7,3,2]: 860.2448249430311 beta[7,3] + gamma[7,3,2] + delta[7,3,2]
   + epsilon[7,3] <= 0
 v[7,3,3]: gamma[7,3,3] + delta[7,3,3] + epsilon[7,3] <= 0
 v[7,3,4]: 224.9593778503113 beta[7,3] + gamma[7,3,4] + delta[7,3,4]
   + epsilon[7,3] <= 0
 v[7,3,5]: 385.6184879322661 beta[7,3] + gamma[7,3,5] + delta[7,3,5]
   + epsilon[7,3] <= 0
 v[7,3,6]: 1405.197716241144 beta[7,3] + gamma[7,3,6] + delta[7,3,6]
   + epsilon[7,3] <= 0
 v[7,3,7]: 1686.083062959439 beta[7,3] + gamma[7,3,7] + delta[7,3,7]
   + epsilon[7,3] <= 0
 v[7,3,8]: 840.3877108467427 beta[7,3] + gamma[7,3,8] + delta[7,3,8]
   + epsilon[7,3] <= 0
 v[7,3,9]: 1318.074813470211 beta[7,3] + gamma[7,3,9] + delta[7,3,9]
   + epsilon[7,3] <= 0
 v[7,4,0]: 70.07365237352232 beta[7,4] + gamma[7,4,0] + delta[7,4,0]
   + epsilon[7,4] <= 0
 v[7,4,1]: 2346.140346718066 beta[7,4] + gamma[7,4,1] + delta[7,4,1]
   + epsilon[7,4] <= 0
 v[7,4,2]: 2627.779486229107 beta[7,4] + gamma[7,4,2] + delta[7,4,2]
   + epsilon[7,4] <= 0
 v[7,4,3]: 359.0023573637434 beta[7,4] + gamma[7,4,3] + delta[7,4,3]
   + epsilon[7,4] <= 0
 v[7,4,4]: gamma[7,4,4] + delta[7,4,4] + epsilon[7,4] <= 0
 v[7,4,5]: 63.78943986733662 beta[7,4] + gamma[7,4,5] + delta[7,4,5]
   + epsilon[7,4] <= 0
 v[7,4,6]: 3077.979552443946 beta[7,4] + gamma[7,4,6] + delta[7,4,6]
   + epsilon[7,4] <= 0
 v[7,4,7]: 2014.862058582201 beta[7,4] + gamma[7,4,7] + delta[7,4,7]
   + epsilon[7,4] <= 0
 v[7,4,8]: 2061.725495539773 beta[7,4] + gamma[7,4,8] + delta[7,4,8]
   + epsilon[7,4] <= 0
 v[7,4,9]: 4397.58172053641 beta[7,4] + gamma[7,4,9] + delta[7,4,9]
   + epsilon[7,4] <= 0
 v[7,5,0]: 27.50740011325044 beta[7,5] + gamma[7,5,0] + delta[7,5,0]
   + epsilon[7,5] <= 0
 v[7,5,1]: 1076.508994242035 beta[7,5] + gamma[7,5,1] + delta[7,5,1]
   + epsilon[7,5] <= 0
 v[7,5,2]: 1244.91997700965 beta[7,5] + gamma[7,5,2] + delta[7,5,2]
   + epsilon[7,5] <= 0
 v[7,5,3]: 194.1665412599118 beta[7,5] + gamma[7,5,3] + delta[7,5,3]
   + epsilon[7,5] <= 0
 v[7,5,4]: 48.04623671234744 beta[7,5] + gamma[7,5,4] + delta[7,5,4]
   + epsilon[7,5] <= 0
 v[7,5,5]: gamma[7,5,5] + delta[7,5,5] + epsilon[7,5] <= 0
 v[7,5,6]: 1742.007242674401 beta[7,5] + gamma[7,5,6] + delta[7,5,6]
   + epsilon[7,5] <= 0
 v[7,5,7]: 981.1480638344783 beta[7,5] + gamma[7,5,7] + delta[7,5,7]
   + epsilon[7,5] <= 0
 v[7,5,8]: 1093.831578439132 beta[7,5] + gamma[7,5,8] + delta[7,5,8]
   + epsilon[7,5] <= 0
 v[7,5,9]: 2466.18668053514 beta[7,5] + gamma[7,5,9] + delta[7,5,9]
   + epsilon[7,5] <= 0
 v[7,6,0]: 3228.948694875151 beta[7,6] + gamma[7,6,0] + delta[7,6,0]
   + epsilon[7,6] <= 0
 v[7,6,1]: 1030.376133195051 beta[7,6] + gamma[7,6,1] + delta[7,6,1]
   + epsilon[7,6] <= 0
 v[7,6,2]: 1043.204946351889 beta[7,6] + gamma[7,6,2] + delta[7,6,2]
   + epsilon[7,6] <= 0
 v[7,6,3]: 1629.820427723519 beta[7,6] + gamma[7,6,3] + delta[7,6,3]
   + epsilon[7,6] <= 0
 v[7,6,4]: 2879.956425255231 beta[7,6] + gamma[7,6,4] + delta[7,6,4]
   + epsilon[7,6] <= 0
 v[7,6,5]: 3958.790895627637 beta[7,6] + gamma[7,6,5] + delta[7,6,5]
   + epsilon[7,6] <= 0
 v[7,6,6]: gamma[7,6,6] + delta[7,6,6] + epsilon[7,6] <= 0
 v[7,6,7]: 1008.351506329531 beta[7,6] + gamma[7,6,7] + delta[7,6,7]
   + epsilon[7,6] <= 0
 v[7,6,8]: 1452.28490969301 beta[7,6] + gamma[7,6,8] + delta[7,6,8]
   + epsilon[7,6] <= 0
 v[7,6,9]: 9.072924027096441 beta[7,6] + gamma[7,6,9] + delta[7,6,9]
   + epsilon[7,6] <= 0
 v[7,7,0]: 1177.694990361826 beta[7,7] + gamma[7,7,0] + delta[7,7,0]
   + epsilon[7,7] <= 0
 v[7,7,1]: 684.8678087158104 beta[7,7] + gamma[7,7,1] + delta[7,7,1]
   + epsilon[7,7] <= 0
 v[7,7,2]: 546.4309104711999 beta[7,7] + gamma[7,7,2] + delta[7,7,2]
   + epsilon[7,7] <= 0
 v[7,7,3]: 2013.249120384493 beta[7,7] + gamma[7,7,3] + delta[7,7,3]
   + epsilon[7,7] <= 0
 v[7,7,4]: 2234.751651355337 beta[7,7] + gamma[7,7,4] + delta[7,7,4]
   + epsilon[7,7] <= 0
 v[7,7,5]: 1630.911586606244 beta[7,7] + gamma[7,7,5] + delta[7,7,5]
   + epsilon[7,7] <= 0
 v[7,7,6]: 820.6070826708792 beta[7,7] + gamma[7,7,6] + delta[7,7,6]
   + epsilon[7,7] <= 0
 v[7,7,7]: gamma[7,7,7] + delta[7,7,7] + epsilon[7,7] <= 0
 v[7,7,8]: 829.6134456235094 beta[7,7] + gamma[7,7,8] + delta[7,7,8]
   + epsilon[7,7] <= 0
 v[7,7,9]: 594.884916117662 beta[7,7] + gamma[7,7,9] + delta[7,7,9]
   + epsilon[7,7] <= 0
 v[7,8,0]: 3107.268960686792 beta[7,8] + gamma[7,8,0] + delta[7,8,0]
   + epsilon[7,8] <= 0
 v[7,8,1]: 1593.312433241548 beta[7,8] + gamma[7,8,1] + delta[7,8,1]
   + epsilon[7,8] <= 0
 v[7,8,2]: 397.6863457849407 beta[7,8] + gamma[7,8,2] + delta[7,8,2]
   + epsilon[7,8] <= 0
 v[7,8,3]: 1014.29643217323 beta[7,8] + gamma[7,8,3] + delta[7,8,3]
   + epsilon[7,8] <= 0
 v[7,8,4]: 2097.601653952857 beta[7,8] + gamma[7,8,4] + delta[7,8,4]
   + epsilon[7,8] <= 0
 v[7,8,5]: 3334.449796237298 beta[7,8] + gamma[7,8,5] + delta[7,8,5]
   + epsilon[7,8] <= 0
 v[7,8,6]: 3717.894759107851 beta[7,8] + gamma[7,8,6] + delta[7,8,6]
   + epsilon[7,8] <= 0
 v[7,8,7]: 2182.762100402356 beta[7,8] + gamma[7,8,7] + delta[7,8,7]
   + epsilon[7,8] <= 0
 v[7,8,8]: gamma[7,8,8] + delta[7,8,8] + epsilon[7,8] <= 0
 v[7,8,9]: 2940.914537032579 beta[7,8] + gamma[7,8,9] + delta[7,8,9]
   + epsilon[7,8] <= 0
 v[7,9,0]: 4980.095005564373 beta[7,9] + gamma[7,9,0] + delta[7,9,0]
   + epsilon[7,9] <= 0
 v[7,9,1]: 1408.363332431962 beta[7,9] + gamma[7,9,1] + delta[7,9,1]
   + epsilon[7,9] <= 0
 v[7,9,2]: 2335.835535311257 beta[7,9] + gamma[7,9,2] + delta[7,9,2]
   + epsilon[7,9] <= 0
 v[7,9,3]: 4448.354972397377 beta[7,9] + gamma[7,9,3] + delta[7,9,3]
   + epsilon[7,9] <= 0
 v[7,9,4]: 6179.032882683153 beta[7,9] + gamma[7,9,4] + delta[7,9,4]
   + epsilon[7,9] <= 0
 v[7,9,5]: 5804.295178069965 beta[7,9] + gamma[7,9,5] + delta[7,9,5]
   + epsilon[7,9] <= 0
 v[7,9,6]: 15.78013704603913 beta[7,9] + gamma[7,9,6] + delta[7,9,6]
   + epsilon[7,9] <= 0
 v[7,9,7]: 2414.475975505661 beta[7,9] + gamma[7,9,7] + delta[7,9,7]
   + epsilon[7,9] <= 0
 v[7,9,8]: 1426.94047735203 beta[7,9] + gamma[7,9,8] + delta[7,9,8]
   + epsilon[7,9] <= 0
 v[7,9,9]: gamma[7,9,9] + delta[7,9,9] + epsilon[7,9] <= 0
 v[8,0,0]: gamma[8,0,0] + delta[8,0,0] + epsilon[8,0] <= 0
 v[8,0,1]: 1697.078543428207 beta[8,0] + gamma[8,0,1] + delta[8,0,1]
   + epsilon[8,0] <= 0
 v[8,0,2]: 2269.086204722981 beta[8,0] + gamma[8,0,2] + delta[8,0,2]
   + epsilon[8,0] <= 0
 v[8,0,3]: 410.2155934652194 beta[8,0] + gamma[8,0,3] + delta[8,0,3]
   + epsilon[8,0] <= 0
 v[8,0,4]: 166.0649942141103 beta[8,0] + gamma[8,0,4] + delta[8,0,4]
   + epsilon[8,0] <= 0
 v[8,0,5]: 39.10355538691039 beta[8,0] + gamma[8,0,5] + delta[8,0,5]
   + epsilon[8,0] <= 0
 v[8,0,6]: 3598.884756032989 beta[8,0] + gamma[8,0,6] + delta[8,0,6]
   + epsilon[8,0] <= 0
 v[8,0,7]: 1608.54872614617 beta[8,0] + gamma[8,0,7] + delta[8,0,7]
   + epsilon[8,0] <= 0
 v[8,0,8]: 1927.16384058249 beta[8,0] + gamma[8,0,8] + delta[8,0,8]
   + epsilon[8,0] <= 0
 v[8,0,9]: 2993.494400008991 beta[8,0] + gamma[8,0,9] + delta[8,0,9]
   + epsilon[8,0] <= 0
 v[8,1,0]: 1090.17554951146 beta[8,1] + gamma[8,1,0] + delta[8,1,0]
   + epsilon[8,1] <= 0
 v[8,1,1]: gamma[8,1,1] + delta[8,1,1] + epsilon[8,1] <= 0
 v[8,1,2]: 598.9145569887431 beta[8,1] + gamma[8,1,2] + delta[8,1,2]
   + epsilon[8,1] <= 0
 v[8,1,3]: 1183.45262848807 beta[8,1] + gamma[8,1,3] + delta[8,1,3]
   + epsilon[8,1] <= 0
 v[8,1,4]: 2157.747373711534 beta[8,1] + gamma[8,1,4] + delta[8,1,4]
   + epsilon[8,1] <= 0
 v[8,1,5]: 1280.576326591262 beta[8,1] + gamma[8,1,5] + delta[8,1,5]
   + epsilon[8,1] <= 0
 v[8,1,6]: 838.9283964458628 beta[8,1] + gamma[8,1,6] + delta[8,1,6]
   + epsilon[8,1] <= 0
 v[8,1,7]: 887.0675245922241 beta[8,1] + gamma[8,1,7] + delta[8,1,7]
   + epsilon[8,1] <= 0
 v[8,1,8]: 844.0621513899163 beta[8,1] + gamma[8,1,8] + delta[8,1,8]
   + epsilon[8,1] <= 0
 v[8,1,9]: 567.4726031397504 beta[8,1] + gamma[8,1,9] + delta[8,1,9]
   + epsilon[8,1] <= 0
 v[8,2,0]: 1930.419160956706 beta[8,2] + gamma[8,2,0] + delta[8,2,0]
   + epsilon[8,2] <= 0
 v[8,2,1]: 2247.778166884799 beta[8,2] + gamma[8,2,1] + delta[8,2,1]
   + epsilon[8,2] <= 0
 v[8,2,2]: gamma[8,2,2] + delta[8,2,2] + epsilon[8,2] <= 0
 v[8,2,3]: 2754.685485176561 beta[8,2] + gamma[8,2,3] + delta[8,2,3]
   + epsilon[8,2] <= 0
 v[8,2,4]: 2581.339969861103 beta[8,2] + gamma[8,2,4] + delta[8,2,4]
   + epsilon[8,2] <= 0
 v[8,2,5]: 2226.793131714933 beta[8,2] + gamma[8,2,5] + delta[8,2,5]
   + epsilon[8,2] <= 0
 v[8,2,6]: 4008.889430030518 beta[8,2] + gamma[8,2,6] + delta[8,2,6]
   + epsilon[8,2] <= 0
 v[8,2,7]: 974.2622659497213 beta[8,2] + gamma[8,2,7] + delta[8,2,7]
   + epsilon[8,2] <= 0
 v[8,2,8]: 374.1002590416991 beta[8,2] + gamma[8,2,8] + delta[8,2,8]
   + epsilon[8,2] <= 0
 v[8,2,9]: 1432.408901074906 beta[8,2] + gamma[8,2,9] + delta[8,2,9]
   + epsilon[8,2] <= 0
 v[8,3,0]: 376.7209080060174 beta[8,3] + gamma[8,3,0] + delta[8,3,0]
   + epsilon[8,3] <= 0
 v[8,3,1]: 986.8426811967194 beta[8,3] + gamma[8,3,1] + delta[8,3,1]
   + epsilon[8,3] <= 0
 v[8,3,2]: 1087.269727834125 beta[8,3] + gamma[8,3,2] + delta[8,3,2]
   + epsilon[8,3] <= 0
 v[8,3,3]: gamma[8,3,3] + delta[8,3,3] + epsilon[8,3] <= 0
 v[8,3,4]: 396.5231588445743 beta[8,3] + gamma[8,3,4] + delta[8,3,4]
   + epsilon[8,3] <= 0
 v[8,3,5]: 313.1109081798591 beta[8,3] + gamma[8,3,5] + delta[8,3,5]
   + epsilon[8,3] <= 0
 v[8,3,6]: 2644.083296094359 beta[8,3] + gamma[8,3,6] + delta[8,3,6]
   + epsilon[8,3] <= 0
 v[8,3,7]: 2356.1411660003 beta[8,3] + gamma[8,3,7] + delta[8,3,7]
   + epsilon[8,3] <= 0
 v[8,3,8]: 857.6056111026179 beta[8,3] + gamma[8,3,8] + delta[8,3,8]
   + epsilon[8,3] <= 0
 v[8,3,9]: 2258.093333331894 beta[8,3] + gamma[8,3,9] + delta[8,3,9]
   + epsilon[8,3] <= 0
 v[8,4,0]: 69.98774348824675 beta[8,4] + gamma[8,4,0] + delta[8,4,0]
   + epsilon[8,4] <= 0
 v[8,4,1]: 1754.412634547963 beta[8,4] + gamma[8,4,1] + delta[8,4,1]
   + epsilon[8,4] <= 0
 v[8,4,2]: 925.0515504829087 beta[8,4] + gamma[8,4,2] + delta[8,4,2]
   + epsilon[8,4] <= 0
 v[8,4,3]: 411.9636930387138 beta[8,4] + gamma[8,4,3] + delta[8,4,3]
   + epsilon[8,4] <= 0
 v[8,4,4]: gamma[8,4,4] + delta[8,4,4] + epsilon[8,4] <= 0
 v[8,4,5]: 70.73302122967914 beta[8,4] + gamma[8,4,5] + delta[8,4,5]
   + epsilon[8,4] <= 0
 v[8,4,6]: 3995.188764663305 beta[8,4] + gamma[8,4,6] + delta[8,4,6]
   + epsilon[8,4] <= 0
 v[8,4,7]: 2532.709185851092 beta[8,4] + gamma[8,4,7] + delta[8,4,7]
   + epsilon[8,4] <= 0
 v[8,4,8]: 1214.579868396084 beta[8,4] + gamma[8,4,8] + delta[8,4,8]
   + epsilon[8,4] <= 0
 v[8,4,9]: 2223.030504615768 beta[8,4] + gamma[8,4,9] + delta[8,4,9]
   + epsilon[8,4] <= 0
 v[8,5,0]: 25.8875574387193 beta[8,5] + gamma[8,5,0] + delta[8,5,0]
   + epsilon[8,5] <= 0
 v[8,5,1]: 2565.287132932669 beta[8,5] + gamma[8,5,1] + delta[8,5,1]
   + epsilon[8,5] <= 0
 v[8,5,2]: 1073.373136205287 beta[8,5] + gamma[8,5,2] + delta[8,5,2]
   + epsilon[8,5] <= 0
 v[8,5,3]: 413.0074077606531 beta[8,5] + gamma[8,5,3] + delta[8,5,3]
   + epsilon[8,5] <= 0
 v[8,5,4]: 101.9620317278275 beta[8,5] + gamma[8,5,4] + delta[8,5,4]
   + epsilon[8,5] <= 0
 v[8,5,5]: gamma[8,5,5] + delta[8,5,5] + epsilon[8,5] <= 0
 v[8,5,6]: 2426.902398280819 beta[8,5] + gamma[8,5,6] + delta[8,5,6]
   + epsilon[8,5] <= 0
 v[8,5,7]: 4035.68149199774 beta[8,5] + gamma[8,5,7] + delta[8,5,7]
   + epsilon[8,5] <= 0
 v[8,5,8]: 1889.926098944639 beta[8,5] + gamma[8,5,8] + delta[8,5,8]
   + epsilon[8,5] <= 0
 v[8,5,9]: 5101.097776833688 beta[8,5] + gamma[8,5,9] + delta[8,5,9]
   + epsilon[8,5] <= 0
 v[8,6,0]: 2439.827349543278 beta[8,6] + gamma[8,6,0] + delta[8,6,0]
   + epsilon[8,6] <= 0
 v[8,6,1]: 684.1590659587209 beta[8,6] + gamma[8,6,1] + delta[8,6,1]
   + epsilon[8,6] <= 0
 v[8,6,2]: 1830.629920128909 beta[8,6] + gamma[8,6,2] + delta[8,6,2]
   + epsilon[8,6] <= 0
 v[8,6,3]: 4255.254296548795 beta[8,6] + gamma[8,6,3] + delta[8,6,3]
   + epsilon[8,6] <= 0
 v[8,6,4]: 2567.896684660011 beta[8,6] + gamma[8,6,4] + delta[8,6,4]
   + epsilon[8,6] <= 0
 v[8,6,5]: 4770.89354848405 beta[8,6] + gamma[8,6,5] + delta[8,6,5]
   + epsilon[8,6] <= 0
 v[8,6,6]: gamma[8,6,6] + delta[8,6,6] + epsilon[8,6] <= 0
 v[8,6,7]: 981.5527739900878 beta[8,6] + gamma[8,6,7] + delta[8,6,7]
   + epsilon[8,6] <= 0
 v[8,6,8]: 1056.139572404526 beta[8,6] + gamma[8,6,8] + delta[8,6,8]
   + epsilon[8,6] <= 0
 v[8,6,9]: 10.68865668633413 beta[8,6] + gamma[8,6,9] + delta[8,6,9]
   + epsilon[8,6] <= 0
 v[8,7,0]: 1610.501809959235 beta[8,7] + gamma[8,7,0] + delta[8,7,0]
   + epsilon[8,7] <= 0
 v[8,7,1]: 967.6746046420474 beta[8,7] + gamma[8,7,1] + delta[8,7,1]
   + epsilon[8,7] <= 0
 v[8,7,2]: 420.4087305729183 beta[8,7] + gamma[8,7,2] + delta[8,7,2]
   + epsilon[8,7] <= 0
 v[8,7,3]: 1003.291105706737 beta[8,7] + gamma[8,7,3] + delta[8,7,3]
   + epsilon[8,7] <= 0
 v[8,7,4]: 911.5033982521924 beta[8,7] + gamma[8,7,4] + delta[8,7,4]
   + epsilon[8,7] <= 0
 v[8,7,5]: 1109.157720045925 beta[8,7] + gamma[8,7,5] + delta[8,7,5]
   + epsilon[8,7] <= 0
 v[8,7,6]: 876.6047836606422 beta[8,7] + gamma[8,7,6] + delta[8,7,6]
   + epsilon[8,7] <= 0
 v[8,7,7]: gamma[8,7,7] + delta[8,7,7] + epsilon[8,7] <= 0
 v[8,7,8]: 692.6710377666409 beta[8,7] + gamma[8,7,8] + delta[8,7,8]
   + epsilon[8,7] <= 0
 v[8,7,9]: 810.7785346866702 beta[8,7] + gamma[8,7,9] + delta[8,7,9]
   + epsilon[8,7] <= 0
 v[8,8,0]: 3896.894167765935 beta[8,8] + gamma[8,8,0] + delta[8,8,0]
   + epsilon[8,8] <= 0
 v[8,8,1]: 2072.095835437371 beta[8,8] + gamma[8,8,1] + delta[8,8,1]
   + epsilon[8,8] <= 0
 v[8,8,2]: 480.5597895995417 beta[8,8] + gamma[8,8,2] + delta[8,8,2]
   + epsilon[8,8] <= 0
 v[8,8,3]: 2628.975979114331 beta[8,8] + gamma[8,8,3] + delta[8,8,3]
   + epsilon[8,8] <= 0
 v[8,8,4]: 1331.794254795681 beta[8,8] + gamma[8,8,4] + delta[8,8,4]
   + epsilon[8,8] <= 0
 v[8,8,5]: 1460.630316024805 beta[8,8] + gamma[8,8,5] + delta[8,8,5]
   + epsilon[8,8] <= 0
 v[8,8,6]: 2854.789146792125 beta[8,8] + gamma[8,8,6] + delta[8,8,6]
   + epsilon[8,8] <= 0
 v[8,8,7]: 1848.141511382149 beta[8,8] + gamma[8,8,7] + delta[8,8,7]
   + epsilon[8,8] <= 0
 v[8,8,8]: gamma[8,8,8] + delta[8,8,8] + epsilon[8,8] <= 0
 v[8,8,9]: 2916.943999962828 beta[8,8] + gamma[8,8,9] + delta[8,8,9]
   + epsilon[8,8] <= 0
 v[8,9,0]: 6033.932746226605 beta[8,9] + gamma[8,9,0] + delta[8,9,0]
   + epsilon[8,9] <= 0
 v[8,9,1]: 2519.367718708036 beta[8,9] + gamma[8,9,1] + delta[8,9,1]
   + epsilon[8,9] <= 0
 v[8,9,2]: 2959.067726837245 beta[8,9] + gamma[8,9,2] + delta[8,9,2]
   + epsilon[8,9] <= 0
 v[8,9,3]: 6920.468874421754 beta[8,9] + gamma[8,9,3] + delta[8,9,3]
   + epsilon[8,9] <= 0
 v[8,9,4]: 4283.732799156965 beta[8,9] + gamma[8,9,4] + delta[8,9,4]
   + epsilon[8,9] <= 0
 v[8,9,5]: 4458.409046053743 beta[8,9] + gamma[8,9,5] + delta[8,9,5]
   + epsilon[8,9] <= 0
 v[8,9,6]: 15.04040942814681 beta[8,9] + gamma[8,9,6] + delta[8,9,6]
   + epsilon[8,9] <= 0
 v[8,9,7]: 2204.826619823533 beta[8,9] + gamma[8,9,7] + delta[8,9,7]
   + epsilon[8,9] <= 0
 v[8,9,8]: 2259.763562284803 beta[8,9] + gamma[8,9,8] + delta[8,9,8]
   + epsilon[8,9] <= 0
 v[8,9,9]: gamma[8,9,9] + delta[8,9,9] + epsilon[8,9] <= 0
 v[9,0,0]: gamma[9,0,0] + delta[9,0,0] + epsilon[9,0] <= 0
 v[9,0,1]: 3314.208360865884 beta[9,0] + gamma[9,0,1] + delta[9,0,1]
   + epsilon[9,0] <= 0
 v[9,0,2]: 3206.962508628642 beta[9,0] + gamma[9,0,2] + delta[9,0,2]
   + epsilon[9,0] <= 0
 v[9,0,3]: 702.5820939783952 beta[9,0] + gamma[9,0,3] + delta[9,0,3]
   + epsilon[9,0] <= 0
 v[9,0,4]: 231.8491239452035 beta[9,0] + gamma[9,0,4] + delta[9,0,4]
   + epsilon[9,0] <= 0
 v[9,0,5]: 100.0695226416946 beta[9,0] + gamma[9,0,5] + delta[9,0,5]
   + epsilon[9,0] <= 0
 v[9,0,6]: 4883.602027541699 beta[9,0] + gamma[9,0,6] + delta[9,0,6]
   + epsilon[9,0] <= 0
 v[9,0,7]: 6099.838360987877 beta[9,0] + gamma[9,0,7] + delta[9,0,7]
   + epsilon[9,0] <= 0
 v[9,0,8]: 2385.338308364218 beta[9,0] + gamma[9,0,8] + delta[9,0,8]
   + epsilon[9,0] <= 0
 v[9,0,9]: 7690.22334433391 beta[9,0] + gamma[9,0,9] + delta[9,0,9]
   + epsilon[9,0] <= 0
 v[9,1,0]: 3041.443721786888 beta[9,1] + gamma[9,1,0] + delta[9,1,0]
   + epsilon[9,1] <= 0
 v[9,1,1]: gamma[9,1,1] + delta[9,1,1] + epsilon[9,1] <= 0
 v[9,1,2]: 2092.772320776616 beta[9,1] + gamma[9,1,2] + delta[9,1,2]
   + epsilon[9,1] <= 0
 v[9,1,3]: 3795.799925881321 beta[9,1] + gamma[9,1,3] + delta[9,1,3]
   + epsilon[9,1] <= 0
 v[9,1,4]: 3483.22728696865 beta[9,1] + gamma[9,1,4] + delta[9,1,4]
   + epsilon[9,1] <= 0
 v[9,1,5]: 2594.061929922373 beta[9,1] + gamma[9,1,5] + delta[9,1,5]
   + epsilon[9,1] <= 0
 v[9,1,6]: 878.1061807158932 beta[9,1] + gamma[9,1,6] + delta[9,1,6]
   + epsilon[9,1] <= 0
 v[9,1,7]: 1705.919010065868 beta[9,1] + gamma[9,1,7] + delta[9,1,7]
   + epsilon[9,1] <= 0
 v[9,1,8]: 803.9878201595994 beta[9,1] + gamma[9,1,8] + delta[9,1,8]
   + epsilon[9,1] <= 0
 v[9,1,9]: 1952.123762250751 beta[9,1] + gamma[9,1,9] + delta[9,1,9]
   + epsilon[9,1] <= 0
 v[9,2,0]: 1870.394297856396 beta[9,2] + gamma[9,2,0] + delta[9,2,0]
   + epsilon[9,2] <= 0
 v[9,2,1]: 914.8460748569266 beta[9,2] + gamma[9,2,1] + delta[9,2,1]
   + epsilon[9,2] <= 0
 v[9,2,2]: gamma[9,2,2] + delta[9,2,2] + epsilon[9,2] <= 0
 v[9,2,3]: 1214.132376606873 beta[9,2] + gamma[9,2,3] + delta[9,2,3]
   + epsilon[9,2] <= 0
 v[9,2,4]: 814.178014094368 beta[9,2] + gamma[9,2,4] + delta[9,2,4]
   + epsilon[9,2] <= 0
 v[9,2,5]: 995.5251837419116 beta[9,2] + gamma[9,2,5] + delta[9,2,5]
   + epsilon[9,2] <= 0
 v[9,2,6]: 1147.588277147037 beta[9,2] + gamma[9,2,6] + delta[9,2,6]
   + epsilon[9,2] <= 0
 v[9,2,7]: 899.9438724438244 beta[9,2] + gamma[9,2,7] + delta[9,2,7]
   + epsilon[9,2] <= 0
 v[9,2,8]: 125.4301670993102 beta[9,2] + gamma[9,2,8] + delta[9,2,8]
   + epsilon[9,2] <= 0
 v[9,2,9]: 1839.609244414822 beta[9,2] + gamma[9,2,9] + delta[9,2,9]
   + epsilon[9,2] <= 0
 v[9,3,0]: 529.2215964634719 beta[9,3] + gamma[9,3,0] + delta[9,3,0]
   + epsilon[9,3] <= 0
 v[9,3,1]: 2409.595444843859 beta[9,3] + gamma[9,3,1] + delta[9,3,1]
   + epsilon[9,3] <= 0
 v[9,3,2]: 1332.734991946144 beta[9,3] + gamma[9,3,2] + delta[9,3,2]
   + epsilon[9,3] <= 0
 v[9,3,3]: gamma[9,3,3] + delta[9,3,3] + epsilon[9,3] <= 0
 v[9,3,4]: 409.1134509435045 beta[9,3] + gamma[9,3,4] + delta[9,3,4]
   + epsilon[9,3] <= 0
 v[9,3,5]: 252.0897563640121 beta[9,3] + gamma[9,3,5] + delta[9,3,5]
   + epsilon[9,3] <= 0
 v[9,3,6]: 4466.9027110576 beta[9,3] + gamma[9,3,6] + delta[9,3,6]
   + epsilon[9,3] <= 0
 v[9,3,7]: 3170.709322140762 beta[9,3] + gamma[9,3,7] + delta[9,3,7]
   + epsilon[9,3] <= 0
 v[9,3,8]: 1695.76525177106 beta[9,3] + gamma[9,3,8] + delta[9,3,8]
   + epsilon[9,3] <= 0
 v[9,3,9]: 4739.088651091251 beta[9,3] + gamma[9,3,9] + delta[9,3,9]
   + epsilon[9,3] <= 0
 v[9,4,0]: 73.48921092731011 beta[9,4] + gamma[9,4,0] + delta[9,4,0]
   + epsilon[9,4] <= 0
 v[9,4,1]: 1839.916603154463 beta[9,4] + gamma[9,4,1] + delta[9,4,1]
   + epsilon[9,4] <= 0
 v[9,4,2]: 1006.492856579538 beta[9,4] + gamma[9,4,2] + delta[9,4,2]
   + epsilon[9,4] <= 0
 v[9,4,3]: 383.8644749646768 beta[9,4] + gamma[9,4,3] + delta[9,4,3]
   + epsilon[9,4] <= 0
 v[9,4,4]: gamma[9,4,4] + delta[9,4,4] + epsilon[9,4] <= 0
 v[9,4,5]: 53.77865607649959 beta[9,4] + gamma[9,4,5] + delta[9,4,5]
   + epsilon[9,4] <= 0
 v[9,4,6]: 3479.314662132228 beta[9,4] + gamma[9,4,6] + delta[9,4,6]
   + epsilon[9,4] <= 0
 v[9,4,7]: 1912.347146160295 beta[9,4] + gamma[9,4,7] + delta[9,4,7]
   + epsilon[9,4] <= 0
 v[9,4,8]: 2109.426344587537 beta[9,4] + gamma[9,4,8] + delta[9,4,8]
   + epsilon[9,4] <= 0
 v[9,4,9]: 4810.627546875166 beta[9,4] + gamma[9,4,9] + delta[9,4,9]
   + epsilon[9,4] <= 0
 v[9,5,0]: 25.60107703255975 beta[9,5] + gamma[9,5,0] + delta[9,5,0]
   + epsilon[9,5] <= 0
 v[9,5,1]: 2869.721628834454 beta[9,5] + gamma[9,5,1] + delta[9,5,1]
   + epsilon[9,5] <= 0
 v[9,5,2]: 1571.789426644968 beta[9,5] + gamma[9,5,2] + delta[9,5,2]
   + epsilon[9,5] <= 0
 v[9,5,3]: 789.4870553537405 beta[9,5] + gamma[9,5,3] + delta[9,5,3]
   + epsilon[9,5] <= 0
 v[9,5,4]: 152.7537052156991 beta[9,5] + gamma[9,5,4] + delta[9,5,4]
   + epsilon[9,5] <= 0
 v[9,5,5]: gamma[9,5,5] + delta[9,5,5] + epsilon[9,5] <= 0
 v[9,5,6]: 6110.100156679437 beta[9,5] + gamma[9,5,6] + delta[9,5,6]
   + epsilon[9,5] <= 0
 v[9,5,7]: 3026.418457078506 beta[9,5] + gamma[9,5,7] + delta[9,5,7]
   + epsilon[9,5] <= 0
 v[9,5,8]: 2750.900321877041 beta[9,5] + gamma[9,5,8] + delta[9,5,8]
   + epsilon[9,5] <= 0
 v[9,5,9]: 3886.893599227221 beta[9,5] + gamma[9,5,9] + delta[9,5,9]
   + epsilon[9,5] <= 0
 v[9,6,0]: 1520.521944867549 beta[9,6] + gamma[9,6,0] + delta[9,6,0]
   + epsilon[9,6] <= 0
 v[9,6,1]: 921.3733960309381 beta[9,6] + gamma[9,6,1] + delta[9,6,1]
   + epsilon[9,6] <= 0
 v[9,6,2]: 1314.174030453407 beta[9,6] + gamma[9,6,2] + delta[9,6,2]
   + epsilon[9,6] <= 0
 v[9,6,3]: 1010.333134334402 beta[9,6] + gamma[9,6,3] + delta[9,6,3]
   + epsilon[9,6] <= 0
 v[9,6,4]: 2773.423320808924 beta[9,6] + gamma[9,6,4] + delta[9,6,4]
   + epsilon[9,6] <= 0
 v[9,6,5]: 2464.468254145219 beta[9,6] + gamma[9,6,5] + delta[9,6,5]
   + epsilon[9,6] <= 0
 v[9,6,6]: gamma[9,6,6] + delta[9,6,6] + epsilon[9,6] <= 0
 v[9,6,7]: 919.7906853307129 beta[9,6] + gamma[9,6,7] + delta[9,6,7]
   + epsilon[9,6] <= 0
 v[9,6,8]: 839.4673401048683 beta[9,6] + gamma[9,6,8] + delta[9,6,8]
   + epsilon[9,6] <= 0
 v[9,6,9]: 7.135095513000208 beta[9,6] + gamma[9,6,9] + delta[9,6,9]
   + epsilon[9,6] <= 0
 v[9,7,0]: 5215.84667252385 beta[9,7] + gamma[9,7,0] + delta[9,7,0]
   + epsilon[9,7] <= 0
 v[9,7,1]: 1003.91107080224 beta[9,7] + gamma[9,7,1] + delta[9,7,1]
   + epsilon[9,7] <= 0
 v[9,7,2]: 1007.789524761943 beta[9,7] + gamma[9,7,2] + delta[9,7,2]
   + epsilon[9,7] <= 0
 v[9,7,3]: 4479.74935523281 beta[9,7] + gamma[9,7,3] + delta[9,7,3]
   + epsilon[9,7] <= 0
 v[9,7,4]: 4283.315018423798 beta[9,7] + gamma[9,7,4] + delta[9,7,4]
   + epsilon[9,7] <= 0
 v[9,7,5]: 4045.093805220036 beta[9,7] + gamma[9,7,5] + delta[9,7,5]
   + epsilon[9,7] <= 0
 v[9,7,6]: 1236.177636811386 beta[9,7] + gamma[9,7,6] + delta[9,7,6]
   + epsilon[9,7] <= 0
 v[9,7,7]: gamma[9,7,7] + delta[9,7,7] + epsilon[9,7] <= 0
 v[9,7,8]: 2513.743273618425 beta[9,7] + gamma[9,7,8] + delta[9,7,8]
   + epsilon[9,7] <= 0
 v[9,7,9]: 1533.225203878871 beta[9,7] + gamma[9,7,9] + delta[9,7,9]
   + epsilon[9,7] <= 0
 v[9,8,0]: 1552.172632598535 beta[9,8] + gamma[9,8,0] + delta[9,8,0]
   + epsilon[9,8] <= 0
 v[9,8,1]: 2286.084441535923 beta[9,8] + gamma[9,8,1] + delta[9,8,1]
   + epsilon[9,8] <= 0
 v[9,8,2]: 296.9600060725494 beta[9,8] + gamma[9,8,2] + delta[9,8,2]
   + epsilon[9,8] <= 0
 v[9,8,3]: 2590.920522553082 beta[9,8] + gamma[9,8,3] + delta[9,8,3]
   + epsilon[9,8] <= 0
 v[9,8,4]: 3593.374069185722 beta[9,8] + gamma[9,8,4] + delta[9,8,4]
   + epsilon[9,8] <= 0
 v[9,8,5]: 2663.322700826857 beta[9,8] + gamma[9,8,5] + delta[9,8,5]
   + epsilon[9,8] <= 0
 v[9,8,6]: 3393.151358775756 beta[9,8] + gamma[9,8,6] + delta[9,8,6]
   + epsilon[9,8] <= 0
 v[9,8,7]: 2356.298393240122 beta[9,8] + gamma[9,8,7] + delta[9,8,7]
   + epsilon[9,8] <= 0
 v[9,8,8]: gamma[9,8,8] + delta[9,8,8] + epsilon[9,8] <= 0
 v[9,8,9]: 4209.604590571375 beta[9,8] + gamma[9,8,9] + delta[9,8,9]
   + epsilon[9,8] <= 0
 v[9,9,0]: 3236.648416829565 beta[9,9] + gamma[9,9,0] + delta[9,9,0]
   + epsilon[9,9] <= 0
 v[9,9,1]: 1781.530745968883 beta[9,9] + gamma[9,9,1] + delta[9,9,1]
   + epsilon[9,9] <= 0
 v[9,9,2]: 2379.662093967128 beta[9,9] + gamma[9,9,2] + delta[9,9,2]
   + epsilon[9,9] <= 0
 v[9,9,3]: 2204.12522157349 beta[9,9] + gamma[9,9,3] + delta[9,9,3]
   + epsilon[9,9] <= 0
 v[9,9,4]: 2946.624294007933 beta[9,9] + gamma[9,9,4] + delta[9,9,4]
   + epsilon[9,9] <= 0
 v[9,9,5]: 4502.430008642989 beta[9,9] + gamma[9,9,5] + delta[9,9,5]
   + epsilon[9,9] <= 0
 v[9,9,6]: 9.57291207874666 beta[9,9] + gamma[9,9,6] + delta[9,9,6]
   + epsilon[9,9] <= 0
 v[9,9,7]: 947.560115105373 beta[9,9] + gamma[9,9,7] + delta[9,9,7]
   + epsilon[9,9] <= 0
 v[9,9,8]: 1480.841672403755 beta[9,9] + gamma[9,9,8] + delta[9,9,8]
   + epsilon[9,9] <= 0
 v[9,9,9]: gamma[9,9,9] + delta[9,9,9] + epsilon[9,9] <= 0
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
 L3[5]: - beta[5,0] - beta[5,1] - beta[5,2] - beta[5,3] - beta[5,4]
   - beta[5,5] - beta[5,6] - beta[5,7] - beta[5,8] - beta[5,9] <= 1
 L3[6]: - beta[6,0] - beta[6,1] - beta[6,2] - beta[6,3] - beta[6,4]
   - beta[6,5] - beta[6,6] - beta[6,7] - beta[6,8] - beta[6,9] <= 1
 L3[7]: - beta[7,0] - beta[7,1] - beta[7,2] - beta[7,3] - beta[7,4]
   - beta[7,5] - beta[7,6] - beta[7,7] - beta[7,8] - beta[7,9] <= 1
 L3[8]: - beta[8,0] - beta[8,1] - beta[8,2] - beta[8,3] - beta[8,4]
   - beta[8,5] - beta[8,6] - beta[8,7] - beta[8,8] - beta[8,9] <= 1
 L3[9]: - beta[9,0] - beta[9,1] - beta[9,2] - beta[9,3] - beta[9,4]
   - beta[9,5] - beta[9,6] - beta[9,7] - beta[9,8] - beta[9,9] <= 1
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
 -infinity <= beta[5,0] <= 0
 -infinity <= beta[5,1] <= 0
 -infinity <= beta[5,2] <= 0
 -infinity <= beta[5,3] <= 0
 -infinity <= beta[5,4] <= 0
 -infinity <= beta[5,5] <= 0
 -infinity <= beta[5,6] <= 0
 -infinity <= beta[5,7] <= 0
 -infinity <= beta[5,8] <= 0
 -infinity <= beta[5,9] <= 0
 -infinity <= beta[6,0] <= 0
 -infinity <= beta[6,1] <= 0
 -infinity <= beta[6,2] <= 0
 -infinity <= beta[6,3] <= 0
 -infinity <= beta[6,4] <= 0
 -infinity <= beta[6,5] <= 0
 -infinity <= beta[6,6] <= 0
 -infinity <= beta[6,7] <= 0
 -infinity <= beta[6,8] <= 0
 -infinity <= beta[6,9] <= 0
 -infinity <= beta[7,0] <= 0
 -infinity <= beta[7,1] <= 0
 -infinity <= beta[7,2] <= 0
 -infinity <= beta[7,3] <= 0
 -infinity <= beta[7,4] <= 0
 -infinity <= beta[7,5] <= 0
 -infinity <= beta[7,6] <= 0
 -infinity <= beta[7,7] <= 0
 -infinity <= beta[7,8] <= 0
 -infinity <= beta[7,9] <= 0
 -infinity <= beta[8,0] <= 0
 -infinity <= beta[8,1] <= 0
 -infinity <= beta[8,2] <= 0
 -infinity <= beta[8,3] <= 0
 -infinity <= beta[8,4] <= 0
 -infinity <= beta[8,5] <= 0
 -infinity <= beta[8,6] <= 0
 -infinity <= beta[8,7] <= 0
 -infinity <= beta[8,8] <= 0
 -infinity <= beta[8,9] <= 0
 -infinity <= beta[9,0] <= 0
 -infinity <= beta[9,1] <= 0
 -infinity <= beta[9,2] <= 0
 -infinity <= beta[9,3] <= 0
 -infinity <= beta[9,4] <= 0
 -infinity <= beta[9,5] <= 0
 -infinity <= beta[9,6] <= 0
 -infinity <= beta[9,7] <= 0
 -infinity <= beta[9,8] <= 0
 -infinity <= beta[9,9] <= 0
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
 -infinity <= gamma[5,0,0] <= 0
 -infinity <= gamma[5,0,1] <= 0
 -infinity <= gamma[5,0,2] <= 0
 -infinity <= gamma[5,0,3] <= 0
 -infinity <= gamma[5,0,4] <= 0
 -infinity <= gamma[5,0,5] <= 0
 -infinity <= gamma[5,0,6] <= 0
 -infinity <= gamma[5,0,7] <= 0
 -infinity <= gamma[5,0,8] <= 0
 -infinity <= gamma[5,0,9] <= 0
 -infinity <= gamma[5,1,0] <= 0
 -infinity <= gamma[5,1,1] <= 0
 -infinity <= gamma[5,1,2] <= 0
 -infinity <= gamma[5,1,3] <= 0
 -infinity <= gamma[5,1,4] <= 0
 -infinity <= gamma[5,1,5] <= 0
 -infinity <= gamma[5,1,6] <= 0
 -infinity <= gamma[5,1,7] <= 0
 -infinity <= gamma[5,1,8] <= 0
 -infinity <= gamma[5,1,9] <= 0
 -infinity <= gamma[5,2,0] <= 0
 -infinity <= gamma[5,2,1] <= 0
 -infinity <= gamma[5,2,2] <= 0
 -infinity <= gamma[5,2,3] <= 0
 -infinity <= gamma[5,2,4] <= 0
 -infinity <= gamma[5,2,5] <= 0
 -infinity <= gamma[5,2,6] <= 0
 -infinity <= gamma[5,2,7] <= 0
 -infinity <= gamma[5,2,8] <= 0
 -infinity <= gamma[5,2,9] <= 0
 -infinity <= gamma[5,3,0] <= 0
 -infinity <= gamma[5,3,1] <= 0
 -infinity <= gamma[5,3,2] <= 0
 -infinity <= gamma[5,3,3] <= 0
 -infinity <= gamma[5,3,4] <= 0
 -infinity <= gamma[5,3,5] <= 0
 -infinity <= gamma[5,3,6] <= 0
 -infinity <= gamma[5,3,7] <= 0
 -infinity <= gamma[5,3,8] <= 0
 -infinity <= gamma[5,3,9] <= 0
 -infinity <= gamma[5,4,0] <= 0
 -infinity <= gamma[5,4,1] <= 0
 -infinity <= gamma[5,4,2] <= 0
 -infinity <= gamma[5,4,3] <= 0
 -infinity <= gamma[5,4,4] <= 0
 -infinity <= gamma[5,4,5] <= 0
 -infinity <= gamma[5,4,6] <= 0
 -infinity <= gamma[5,4,7] <= 0
 -infinity <= gamma[5,4,8] <= 0
 -infinity <= gamma[5,4,9] <= 0
 -infinity <= gamma[5,5,0] <= 0
 -infinity <= gamma[5,5,1] <= 0
 -infinity <= gamma[5,5,2] <= 0
 -infinity <= gamma[5,5,3] <= 0
 -infinity <= gamma[5,5,4] <= 0
 -infinity <= gamma[5,5,5] <= 0
 -infinity <= gamma[5,5,6] <= 0
 -infinity <= gamma[5,5,7] <= 0
 -infinity <= gamma[5,5,8] <= 0
 -infinity <= gamma[5,5,9] <= 0
 -infinity <= gamma[5,6,0] <= 0
 -infinity <= gamma[5,6,1] <= 0
 -infinity <= gamma[5,6,2] <= 0
 -infinity <= gamma[5,6,3] <= 0
 -infinity <= gamma[5,6,4] <= 0
 -infinity <= gamma[5,6,5] <= 0
 -infinity <= gamma[5,6,6] <= 0
 -infinity <= gamma[5,6,7] <= 0
 -infinity <= gamma[5,6,8] <= 0
 -infinity <= gamma[5,6,9] <= 0
 -infinity <= gamma[5,7,0] <= 0
 -infinity <= gamma[5,7,1] <= 0
 -infinity <= gamma[5,7,2] <= 0
 -infinity <= gamma[5,7,3] <= 0
 -infinity <= gamma[5,7,4] <= 0
 -infinity <= gamma[5,7,5] <= 0
 -infinity <= gamma[5,7,6] <= 0
 -infinity <= gamma[5,7,7] <= 0
 -infinity <= gamma[5,7,8] <= 0
 -infinity <= gamma[5,7,9] <= 0
 -infinity <= gamma[5,8,0] <= 0
 -infinity <= gamma[5,8,1] <= 0
 -infinity <= gamma[5,8,2] <= 0
 -infinity <= gamma[5,8,3] <= 0
 -infinity <= gamma[5,8,4] <= 0
 -infinity <= gamma[5,8,5] <= 0
 -infinity <= gamma[5,8,6] <= 0
 -infinity <= gamma[5,8,7] <= 0
 -infinity <= gamma[5,8,8] <= 0
 -infinity <= gamma[5,8,9] <= 0
 -infinity <= gamma[5,9,0] <= 0
 -infinity <= gamma[5,9,1] <= 0
 -infinity <= gamma[5,9,2] <= 0
 -infinity <= gamma[5,9,3] <= 0
 -infinity <= gamma[5,9,4] <= 0
 -infinity <= gamma[5,9,5] <= 0
 -infinity <= gamma[5,9,6] <= 0
 -infinity <= gamma[5,9,7] <= 0
 -infinity <= gamma[5,9,8] <= 0
 -infinity <= gamma[5,9,9] <= 0
 -infinity <= gamma[6,0,0] <= 0
 -infinity <= gamma[6,0,1] <= 0
 -infinity <= gamma[6,0,2] <= 0
 -infinity <= gamma[6,0,3] <= 0
 -infinity <= gamma[6,0,4] <= 0
 -infinity <= gamma[6,0,5] <= 0
 -infinity <= gamma[6,0,6] <= 0
 -infinity <= gamma[6,0,7] <= 0
 -infinity <= gamma[6,0,8] <= 0
 -infinity <= gamma[6,0,9] <= 0
 -infinity <= gamma[6,1,0] <= 0
 -infinity <= gamma[6,1,1] <= 0
 -infinity <= gamma[6,1,2] <= 0
 -infinity <= gamma[6,1,3] <= 0
 -infinity <= gamma[6,1,4] <= 0
 -infinity <= gamma[6,1,5] <= 0
 -infinity <= gamma[6,1,6] <= 0
 -infinity <= gamma[6,1,7] <= 0
 -infinity <= gamma[6,1,8] <= 0
 -infinity <= gamma[6,1,9] <= 0
 -infinity <= gamma[6,2,0] <= 0
 -infinity <= gamma[6,2,1] <= 0
 -infinity <= gamma[6,2,2] <= 0
 -infinity <= gamma[6,2,3] <= 0
 -infinity <= gamma[6,2,4] <= 0
 -infinity <= gamma[6,2,5] <= 0
 -infinity <= gamma[6,2,6] <= 0
 -infinity <= gamma[6,2,7] <= 0
 -infinity <= gamma[6,2,8] <= 0
 -infinity <= gamma[6,2,9] <= 0
 -infinity <= gamma[6,3,0] <= 0
 -infinity <= gamma[6,3,1] <= 0
 -infinity <= gamma[6,3,2] <= 0
 -infinity <= gamma[6,3,3] <= 0
 -infinity <= gamma[6,3,4] <= 0
 -infinity <= gamma[6,3,5] <= 0
 -infinity <= gamma[6,3,6] <= 0
 -infinity <= gamma[6,3,7] <= 0
 -infinity <= gamma[6,3,8] <= 0
 -infinity <= gamma[6,3,9] <= 0
 -infinity <= gamma[6,4,0] <= 0
 -infinity <= gamma[6,4,1] <= 0
 -infinity <= gamma[6,4,2] <= 0
 -infinity <= gamma[6,4,3] <= 0
 -infinity <= gamma[6,4,4] <= 0
 -infinity <= gamma[6,4,5] <= 0
 -infinity <= gamma[6,4,6] <= 0
 -infinity <= gamma[6,4,7] <= 0
 -infinity <= gamma[6,4,8] <= 0
 -infinity <= gamma[6,4,9] <= 0
 -infinity <= gamma[6,5,0] <= 0
 -infinity <= gamma[6,5,1] <= 0
 -infinity <= gamma[6,5,2] <= 0
 -infinity <= gamma[6,5,3] <= 0
 -infinity <= gamma[6,5,4] <= 0
 -infinity <= gamma[6,5,5] <= 0
 -infinity <= gamma[6,5,6] <= 0
 -infinity <= gamma[6,5,7] <= 0
 -infinity <= gamma[6,5,8] <= 0
 -infinity <= gamma[6,5,9] <= 0
 -infinity <= gamma[6,6,0] <= 0
 -infinity <= gamma[6,6,1] <= 0
 -infinity <= gamma[6,6,2] <= 0
 -infinity <= gamma[6,6,3] <= 0
 -infinity <= gamma[6,6,4] <= 0
 -infinity <= gamma[6,6,5] <= 0
 -infinity <= gamma[6,6,6] <= 0
 -infinity <= gamma[6,6,7] <= 0
 -infinity <= gamma[6,6,8] <= 0
 -infinity <= gamma[6,6,9] <= 0
 -infinity <= gamma[6,7,0] <= 0
 -infinity <= gamma[6,7,1] <= 0
 -infinity <= gamma[6,7,2] <= 0
 -infinity <= gamma[6,7,3] <= 0
 -infinity <= gamma[6,7,4] <= 0
 -infinity <= gamma[6,7,5] <= 0
 -infinity <= gamma[6,7,6] <= 0
 -infinity <= gamma[6,7,7] <= 0
 -infinity <= gamma[6,7,8] <= 0
 -infinity <= gamma[6,7,9] <= 0
 -infinity <= gamma[6,8,0] <= 0
 -infinity <= gamma[6,8,1] <= 0
 -infinity <= gamma[6,8,2] <= 0
 -infinity <= gamma[6,8,3] <= 0
 -infinity <= gamma[6,8,4] <= 0
 -infinity <= gamma[6,8,5] <= 0
 -infinity <= gamma[6,8,6] <= 0
 -infinity <= gamma[6,8,7] <= 0
 -infinity <= gamma[6,8,8] <= 0
 -infinity <= gamma[6,8,9] <= 0
 -infinity <= gamma[6,9,0] <= 0
 -infinity <= gamma[6,9,1] <= 0
 -infinity <= gamma[6,9,2] <= 0
 -infinity <= gamma[6,9,3] <= 0
 -infinity <= gamma[6,9,4] <= 0
 -infinity <= gamma[6,9,5] <= 0
 -infinity <= gamma[6,9,6] <= 0
 -infinity <= gamma[6,9,7] <= 0
 -infinity <= gamma[6,9,8] <= 0
 -infinity <= gamma[6,9,9] <= 0
 -infinity <= gamma[7,0,0] <= 0
 -infinity <= gamma[7,0,1] <= 0
 -infinity <= gamma[7,0,2] <= 0
 -infinity <= gamma[7,0,3] <= 0
 -infinity <= gamma[7,0,4] <= 0
 -infinity <= gamma[7,0,5] <= 0
 -infinity <= gamma[7,0,6] <= 0
 -infinity <= gamma[7,0,7] <= 0
 -infinity <= gamma[7,0,8] <= 0
 -infinity <= gamma[7,0,9] <= 0
 -infinity <= gamma[7,1,0] <= 0
 -infinity <= gamma[7,1,1] <= 0
 -infinity <= gamma[7,1,2] <= 0
 -infinity <= gamma[7,1,3] <= 0
 -infinity <= gamma[7,1,4] <= 0
 -infinity <= gamma[7,1,5] <= 0
 -infinity <= gamma[7,1,6] <= 0
 -infinity <= gamma[7,1,7] <= 0
 -infinity <= gamma[7,1,8] <= 0
 -infinity <= gamma[7,1,9] <= 0
 -infinity <= gamma[7,2,0] <= 0
 -infinity <= gamma[7,2,1] <= 0
 -infinity <= gamma[7,2,2] <= 0
 -infinity <= gamma[7,2,3] <= 0
 -infinity <= gamma[7,2,4] <= 0
 -infinity <= gamma[7,2,5] <= 0
 -infinity <= gamma[7,2,6] <= 0
 -infinity <= gamma[7,2,7] <= 0
 -infinity <= gamma[7,2,8] <= 0
 -infinity <= gamma[7,2,9] <= 0
 -infinity <= gamma[7,3,0] <= 0
 -infinity <= gamma[7,3,1] <= 0
 -infinity <= gamma[7,3,2] <= 0
 -infinity <= gamma[7,3,3] <= 0
 -infinity <= gamma[7,3,4] <= 0
 -infinity <= gamma[7,3,5] <= 0
 -infinity <= gamma[7,3,6] <= 0
 -infinity <= gamma[7,3,7] <= 0
 -infinity <= gamma[7,3,8] <= 0
 -infinity <= gamma[7,3,9] <= 0
 -infinity <= gamma[7,4,0] <= 0
 -infinity <= gamma[7,4,1] <= 0
 -infinity <= gamma[7,4,2] <= 0
 -infinity <= gamma[7,4,3] <= 0
 -infinity <= gamma[7,4,4] <= 0
 -infinity <= gamma[7,4,5] <= 0
 -infinity <= gamma[7,4,6] <= 0
 -infinity <= gamma[7,4,7] <= 0
 -infinity <= gamma[7,4,8] <= 0
 -infinity <= gamma[7,4,9] <= 0
 -infinity <= gamma[7,5,0] <= 0
 -infinity <= gamma[7,5,1] <= 0
 -infinity <= gamma[7,5,2] <= 0
 -infinity <= gamma[7,5,3] <= 0
 -infinity <= gamma[7,5,4] <= 0
 -infinity <= gamma[7,5,5] <= 0
 -infinity <= gamma[7,5,6] <= 0
 -infinity <= gamma[7,5,7] <= 0
 -infinity <= gamma[7,5,8] <= 0
 -infinity <= gamma[7,5,9] <= 0
 -infinity <= gamma[7,6,0] <= 0
 -infinity <= gamma[7,6,1] <= 0
 -infinity <= gamma[7,6,2] <= 0
 -infinity <= gamma[7,6,3] <= 0
 -infinity <= gamma[7,6,4] <= 0
 -infinity <= gamma[7,6,5] <= 0
 -infinity <= gamma[7,6,6] <= 0
 -infinity <= gamma[7,6,7] <= 0
 -infinity <= gamma[7,6,8] <= 0
 -infinity <= gamma[7,6,9] <= 0
 -infinity <= gamma[7,7,0] <= 0
 -infinity <= gamma[7,7,1] <= 0
 -infinity <= gamma[7,7,2] <= 0
 -infinity <= gamma[7,7,3] <= 0
 -infinity <= gamma[7,7,4] <= 0
 -infinity <= gamma[7,7,5] <= 0
 -infinity <= gamma[7,7,6] <= 0
 -infinity <= gamma[7,7,7] <= 0
 -infinity <= gamma[7,7,8] <= 0
 -infinity <= gamma[7,7,9] <= 0
 -infinity <= gamma[7,8,0] <= 0
 -infinity <= gamma[7,8,1] <= 0
 -infinity <= gamma[7,8,2] <= 0
 -infinity <= gamma[7,8,3] <= 0
 -infinity <= gamma[7,8,4] <= 0
 -infinity <= gamma[7,8,5] <= 0
 -infinity <= gamma[7,8,6] <= 0
 -infinity <= gamma[7,8,7] <= 0
 -infinity <= gamma[7,8,8] <= 0
 -infinity <= gamma[7,8,9] <= 0
 -infinity <= gamma[7,9,0] <= 0
 -infinity <= gamma[7,9,1] <= 0
 -infinity <= gamma[7,9,2] <= 0
 -infinity <= gamma[7,9,3] <= 0
 -infinity <= gamma[7,9,4] <= 0
 -infinity <= gamma[7,9,5] <= 0
 -infinity <= gamma[7,9,6] <= 0
 -infinity <= gamma[7,9,7] <= 0
 -infinity <= gamma[7,9,8] <= 0
 -infinity <= gamma[7,9,9] <= 0
 -infinity <= gamma[8,0,0] <= 0
 -infinity <= gamma[8,0,1] <= 0
 -infinity <= gamma[8,0,2] <= 0
 -infinity <= gamma[8,0,3] <= 0
 -infinity <= gamma[8,0,4] <= 0
 -infinity <= gamma[8,0,5] <= 0
 -infinity <= gamma[8,0,6] <= 0
 -infinity <= gamma[8,0,7] <= 0
 -infinity <= gamma[8,0,8] <= 0
 -infinity <= gamma[8,0,9] <= 0
 -infinity <= gamma[8,1,0] <= 0
 -infinity <= gamma[8,1,1] <= 0
 -infinity <= gamma[8,1,2] <= 0
 -infinity <= gamma[8,1,3] <= 0
 -infinity <= gamma[8,1,4] <= 0
 -infinity <= gamma[8,1,5] <= 0
 -infinity <= gamma[8,1,6] <= 0
 -infinity <= gamma[8,1,7] <= 0
 -infinity <= gamma[8,1,8] <= 0
 -infinity <= gamma[8,1,9] <= 0
 -infinity <= gamma[8,2,0] <= 0
 -infinity <= gamma[8,2,1] <= 0
 -infinity <= gamma[8,2,2] <= 0
 -infinity <= gamma[8,2,3] <= 0
 -infinity <= gamma[8,2,4] <= 0
 -infinity <= gamma[8,2,5] <= 0
 -infinity <= gamma[8,2,6] <= 0
 -infinity <= gamma[8,2,7] <= 0
 -infinity <= gamma[8,2,8] <= 0
 -infinity <= gamma[8,2,9] <= 0
 -infinity <= gamma[8,3,0] <= 0
 -infinity <= gamma[8,3,1] <= 0
 -infinity <= gamma[8,3,2] <= 0
 -infinity <= gamma[8,3,3] <= 0
 -infinity <= gamma[8,3,4] <= 0
 -infinity <= gamma[8,3,5] <= 0
 -infinity <= gamma[8,3,6] <= 0
 -infinity <= gamma[8,3,7] <= 0
 -infinity <= gamma[8,3,8] <= 0
 -infinity <= gamma[8,3,9] <= 0
 -infinity <= gamma[8,4,0] <= 0
 -infinity <= gamma[8,4,1] <= 0
 -infinity <= gamma[8,4,2] <= 0
 -infinity <= gamma[8,4,3] <= 0
 -infinity <= gamma[8,4,4] <= 0
 -infinity <= gamma[8,4,5] <= 0
 -infinity <= gamma[8,4,6] <= 0
 -infinity <= gamma[8,4,7] <= 0
 -infinity <= gamma[8,4,8] <= 0
 -infinity <= gamma[8,4,9] <= 0
 -infinity <= gamma[8,5,0] <= 0
 -infinity <= gamma[8,5,1] <= 0
 -infinity <= gamma[8,5,2] <= 0
 -infinity <= gamma[8,5,3] <= 0
 -infinity <= gamma[8,5,4] <= 0
 -infinity <= gamma[8,5,5] <= 0
 -infinity <= gamma[8,5,6] <= 0
 -infinity <= gamma[8,5,7] <= 0
 -infinity <= gamma[8,5,8] <= 0
 -infinity <= gamma[8,5,9] <= 0
 -infinity <= gamma[8,6,0] <= 0
 -infinity <= gamma[8,6,1] <= 0
 -infinity <= gamma[8,6,2] <= 0
 -infinity <= gamma[8,6,3] <= 0
 -infinity <= gamma[8,6,4] <= 0
 -infinity <= gamma[8,6,5] <= 0
 -infinity <= gamma[8,6,6] <= 0
 -infinity <= gamma[8,6,7] <= 0
 -infinity <= gamma[8,6,8] <= 0
 -infinity <= gamma[8,6,9] <= 0
 -infinity <= gamma[8,7,0] <= 0
 -infinity <= gamma[8,7,1] <= 0
 -infinity <= gamma[8,7,2] <= 0
 -infinity <= gamma[8,7,3] <= 0
 -infinity <= gamma[8,7,4] <= 0
 -infinity <= gamma[8,7,5] <= 0
 -infinity <= gamma[8,7,6] <= 0
 -infinity <= gamma[8,7,7] <= 0
 -infinity <= gamma[8,7,8] <= 0
 -infinity <= gamma[8,7,9] <= 0
 -infinity <= gamma[8,8,0] <= 0
 -infinity <= gamma[8,8,1] <= 0
 -infinity <= gamma[8,8,2] <= 0
 -infinity <= gamma[8,8,3] <= 0
 -infinity <= gamma[8,8,4] <= 0
 -infinity <= gamma[8,8,5] <= 0
 -infinity <= gamma[8,8,6] <= 0
 -infinity <= gamma[8,8,7] <= 0
 -infinity <= gamma[8,8,8] <= 0
 -infinity <= gamma[8,8,9] <= 0
 -infinity <= gamma[8,9,0] <= 0
 -infinity <= gamma[8,9,1] <= 0
 -infinity <= gamma[8,9,2] <= 0
 -infinity <= gamma[8,9,3] <= 0
 -infinity <= gamma[8,9,4] <= 0
 -infinity <= gamma[8,9,5] <= 0
 -infinity <= gamma[8,9,6] <= 0
 -infinity <= gamma[8,9,7] <= 0
 -infinity <= gamma[8,9,8] <= 0
 -infinity <= gamma[8,9,9] <= 0
 -infinity <= gamma[9,0,0] <= 0
 -infinity <= gamma[9,0,1] <= 0
 -infinity <= gamma[9,0,2] <= 0
 -infinity <= gamma[9,0,3] <= 0
 -infinity <= gamma[9,0,4] <= 0
 -infinity <= gamma[9,0,5] <= 0
 -infinity <= gamma[9,0,6] <= 0
 -infinity <= gamma[9,0,7] <= 0
 -infinity <= gamma[9,0,8] <= 0
 -infinity <= gamma[9,0,9] <= 0
 -infinity <= gamma[9,1,0] <= 0
 -infinity <= gamma[9,1,1] <= 0
 -infinity <= gamma[9,1,2] <= 0
 -infinity <= gamma[9,1,3] <= 0
 -infinity <= gamma[9,1,4] <= 0
 -infinity <= gamma[9,1,5] <= 0
 -infinity <= gamma[9,1,6] <= 0
 -infinity <= gamma[9,1,7] <= 0
 -infinity <= gamma[9,1,8] <= 0
 -infinity <= gamma[9,1,9] <= 0
 -infinity <= gamma[9,2,0] <= 0
 -infinity <= gamma[9,2,1] <= 0
 -infinity <= gamma[9,2,2] <= 0
 -infinity <= gamma[9,2,3] <= 0
 -infinity <= gamma[9,2,4] <= 0
 -infinity <= gamma[9,2,5] <= 0
 -infinity <= gamma[9,2,6] <= 0
 -infinity <= gamma[9,2,7] <= 0
 -infinity <= gamma[9,2,8] <= 0
 -infinity <= gamma[9,2,9] <= 0
 -infinity <= gamma[9,3,0] <= 0
 -infinity <= gamma[9,3,1] <= 0
 -infinity <= gamma[9,3,2] <= 0
 -infinity <= gamma[9,3,3] <= 0
 -infinity <= gamma[9,3,4] <= 0
 -infinity <= gamma[9,3,5] <= 0
 -infinity <= gamma[9,3,6] <= 0
 -infinity <= gamma[9,3,7] <= 0
 -infinity <= gamma[9,3,8] <= 0
 -infinity <= gamma[9,3,9] <= 0
 -infinity <= gamma[9,4,0] <= 0
 -infinity <= gamma[9,4,1] <= 0
 -infinity <= gamma[9,4,2] <= 0
 -infinity <= gamma[9,4,3] <= 0
 -infinity <= gamma[9,4,4] <= 0
 -infinity <= gamma[9,4,5] <= 0
 -infinity <= gamma[9,4,6] <= 0
 -infinity <= gamma[9,4,7] <= 0
 -infinity <= gamma[9,4,8] <= 0
 -infinity <= gamma[9,4,9] <= 0
 -infinity <= gamma[9,5,0] <= 0
 -infinity <= gamma[9,5,1] <= 0
 -infinity <= gamma[9,5,2] <= 0
 -infinity <= gamma[9,5,3] <= 0
 -infinity <= gamma[9,5,4] <= 0
 -infinity <= gamma[9,5,5] <= 0
 -infinity <= gamma[9,5,6] <= 0
 -infinity <= gamma[9,5,7] <= 0
 -infinity <= gamma[9,5,8] <= 0
 -infinity <= gamma[9,5,9] <= 0
 -infinity <= gamma[9,6,0] <= 0
 -infinity <= gamma[9,6,1] <= 0
 -infinity <= gamma[9,6,2] <= 0
 -infinity <= gamma[9,6,3] <= 0
 -infinity <= gamma[9,6,4] <= 0
 -infinity <= gamma[9,6,5] <= 0
 -infinity <= gamma[9,6,6] <= 0
 -infinity <= gamma[9,6,7] <= 0
 -infinity <= gamma[9,6,8] <= 0
 -infinity <= gamma[9,6,9] <= 0
 -infinity <= gamma[9,7,0] <= 0
 -infinity <= gamma[9,7,1] <= 0
 -infinity <= gamma[9,7,2] <= 0
 -infinity <= gamma[9,7,3] <= 0
 -infinity <= gamma[9,7,4] <= 0
 -infinity <= gamma[9,7,5] <= 0
 -infinity <= gamma[9,7,6] <= 0
 -infinity <= gamma[9,7,7] <= 0
 -infinity <= gamma[9,7,8] <= 0
 -infinity <= gamma[9,7,9] <= 0
 -infinity <= gamma[9,8,0] <= 0
 -infinity <= gamma[9,8,1] <= 0
 -infinity <= gamma[9,8,2] <= 0
 -infinity <= gamma[9,8,3] <= 0
 -infinity <= gamma[9,8,4] <= 0
 -infinity <= gamma[9,8,5] <= 0
 -infinity <= gamma[9,8,6] <= 0
 -infinity <= gamma[9,8,7] <= 0
 -infinity <= gamma[9,8,8] <= 0
 -infinity <= gamma[9,8,9] <= 0
 -infinity <= gamma[9,9,0] <= 0
 -infinity <= gamma[9,9,1] <= 0
 -infinity <= gamma[9,9,2] <= 0
 -infinity <= gamma[9,9,3] <= 0
 -infinity <= gamma[9,9,4] <= 0
 -infinity <= gamma[9,9,5] <= 0
 -infinity <= gamma[9,9,6] <= 0
 -infinity <= gamma[9,9,7] <= 0
 -infinity <= gamma[9,9,8] <= 0
 -infinity <= gamma[9,9,9] <= 0
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
 -infinity <= delta[5,0,0] <= 0
 -infinity <= delta[5,0,1] <= 0
 -infinity <= delta[5,0,2] <= 0
 -infinity <= delta[5,0,3] <= 0
 -infinity <= delta[5,0,4] <= 0
 -infinity <= delta[5,0,5] <= 0
 -infinity <= delta[5,0,6] <= 0
 -infinity <= delta[5,0,7] <= 0
 -infinity <= delta[5,0,8] <= 0
 -infinity <= delta[5,0,9] <= 0
 -infinity <= delta[5,1,0] <= 0
 -infinity <= delta[5,1,1] <= 0
 -infinity <= delta[5,1,2] <= 0
 -infinity <= delta[5,1,3] <= 0
 -infinity <= delta[5,1,4] <= 0
 -infinity <= delta[5,1,5] <= 0
 -infinity <= delta[5,1,6] <= 0
 -infinity <= delta[5,1,7] <= 0
 -infinity <= delta[5,1,8] <= 0
 -infinity <= delta[5,1,9] <= 0
 -infinity <= delta[5,2,0] <= 0
 -infinity <= delta[5,2,1] <= 0
 -infinity <= delta[5,2,2] <= 0
 -infinity <= delta[5,2,3] <= 0
 -infinity <= delta[5,2,4] <= 0
 -infinity <= delta[5,2,5] <= 0
 -infinity <= delta[5,2,6] <= 0
 -infinity <= delta[5,2,7] <= 0
 -infinity <= delta[5,2,8] <= 0
 -infinity <= delta[5,2,9] <= 0
 -infinity <= delta[5,3,0] <= 0
 -infinity <= delta[5,3,1] <= 0
 -infinity <= delta[5,3,2] <= 0
 -infinity <= delta[5,3,3] <= 0
 -infinity <= delta[5,3,4] <= 0
 -infinity <= delta[5,3,5] <= 0
 -infinity <= delta[5,3,6] <= 0
 -infinity <= delta[5,3,7] <= 0
 -infinity <= delta[5,3,8] <= 0
 -infinity <= delta[5,3,9] <= 0
 -infinity <= delta[5,4,0] <= 0
 -infinity <= delta[5,4,1] <= 0
 -infinity <= delta[5,4,2] <= 0
 -infinity <= delta[5,4,3] <= 0
 -infinity <= delta[5,4,4] <= 0
 -infinity <= delta[5,4,5] <= 0
 -infinity <= delta[5,4,6] <= 0
 -infinity <= delta[5,4,7] <= 0
 -infinity <= delta[5,4,8] <= 0
 -infinity <= delta[5,4,9] <= 0
 -infinity <= delta[5,5,0] <= 0
 -infinity <= delta[5,5,1] <= 0
 -infinity <= delta[5,5,2] <= 0
 -infinity <= delta[5,5,3] <= 0
 -infinity <= delta[5,5,4] <= 0
 -infinity <= delta[5,5,5] <= 0
 -infinity <= delta[5,5,6] <= 0
 -infinity <= delta[5,5,7] <= 0
 -infinity <= delta[5,5,8] <= 0
 -infinity <= delta[5,5,9] <= 0
 -infinity <= delta[5,6,0] <= 0
 -infinity <= delta[5,6,1] <= 0
 -infinity <= delta[5,6,2] <= 0
 -infinity <= delta[5,6,3] <= 0
 -infinity <= delta[5,6,4] <= 0
 -infinity <= delta[5,6,5] <= 0
 -infinity <= delta[5,6,6] <= 0
 -infinity <= delta[5,6,7] <= 0
 -infinity <= delta[5,6,8] <= 0
 -infinity <= delta[5,6,9] <= 0
 -infinity <= delta[5,7,0] <= 0
 -infinity <= delta[5,7,1] <= 0
 -infinity <= delta[5,7,2] <= 0
 -infinity <= delta[5,7,3] <= 0
 -infinity <= delta[5,7,4] <= 0
 -infinity <= delta[5,7,5] <= 0
 -infinity <= delta[5,7,6] <= 0
 -infinity <= delta[5,7,7] <= 0
 -infinity <= delta[5,7,8] <= 0
 -infinity <= delta[5,7,9] <= 0
 -infinity <= delta[5,8,0] <= 0
 -infinity <= delta[5,8,1] <= 0
 -infinity <= delta[5,8,2] <= 0
 -infinity <= delta[5,8,3] <= 0
 -infinity <= delta[5,8,4] <= 0
 -infinity <= delta[5,8,5] <= 0
 -infinity <= delta[5,8,6] <= 0
 -infinity <= delta[5,8,7] <= 0
 -infinity <= delta[5,8,8] <= 0
 -infinity <= delta[5,8,9] <= 0
 -infinity <= delta[5,9,0] <= 0
 -infinity <= delta[5,9,1] <= 0
 -infinity <= delta[5,9,2] <= 0
 -infinity <= delta[5,9,3] <= 0
 -infinity <= delta[5,9,4] <= 0
 -infinity <= delta[5,9,5] <= 0
 -infinity <= delta[5,9,6] <= 0
 -infinity <= delta[5,9,7] <= 0
 -infinity <= delta[5,9,8] <= 0
 -infinity <= delta[5,9,9] <= 0
 -infinity <= delta[6,0,0] <= 0
 -infinity <= delta[6,0,1] <= 0
 -infinity <= delta[6,0,2] <= 0
 -infinity <= delta[6,0,3] <= 0
 -infinity <= delta[6,0,4] <= 0
 -infinity <= delta[6,0,5] <= 0
 -infinity <= delta[6,0,6] <= 0
 -infinity <= delta[6,0,7] <= 0
 -infinity <= delta[6,0,8] <= 0
 -infinity <= delta[6,0,9] <= 0
 -infinity <= delta[6,1,0] <= 0
 -infinity <= delta[6,1,1] <= 0
 -infinity <= delta[6,1,2] <= 0
 -infinity <= delta[6,1,3] <= 0
 -infinity <= delta[6,1,4] <= 0
 -infinity <= delta[6,1,5] <= 0
 -infinity <= delta[6,1,6] <= 0
 -infinity <= delta[6,1,7] <= 0
 -infinity <= delta[6,1,8] <= 0
 -infinity <= delta[6,1,9] <= 0
 -infinity <= delta[6,2,0] <= 0
 -infinity <= delta[6,2,1] <= 0
 -infinity <= delta[6,2,2] <= 0
 -infinity <= delta[6,2,3] <= 0
 -infinity <= delta[6,2,4] <= 0
 -infinity <= delta[6,2,5] <= 0
 -infinity <= delta[6,2,6] <= 0
 -infinity <= delta[6,2,7] <= 0
 -infinity <= delta[6,2,8] <= 0
 -infinity <= delta[6,2,9] <= 0
 -infinity <= delta[6,3,0] <= 0
 -infinity <= delta[6,3,1] <= 0
 -infinity <= delta[6,3,2] <= 0
 -infinity <= delta[6,3,3] <= 0
 -infinity <= delta[6,3,4] <= 0
 -infinity <= delta[6,3,5] <= 0
 -infinity <= delta[6,3,6] <= 0
 -infinity <= delta[6,3,7] <= 0
 -infinity <= delta[6,3,8] <= 0
 -infinity <= delta[6,3,9] <= 0
 -infinity <= delta[6,4,0] <= 0
 -infinity <= delta[6,4,1] <= 0
 -infinity <= delta[6,4,2] <= 0
 -infinity <= delta[6,4,3] <= 0
 -infinity <= delta[6,4,4] <= 0
 -infinity <= delta[6,4,5] <= 0
 -infinity <= delta[6,4,6] <= 0
 -infinity <= delta[6,4,7] <= 0
 -infinity <= delta[6,4,8] <= 0
 -infinity <= delta[6,4,9] <= 0
 -infinity <= delta[6,5,0] <= 0
 -infinity <= delta[6,5,1] <= 0
 -infinity <= delta[6,5,2] <= 0
 -infinity <= delta[6,5,3] <= 0
 -infinity <= delta[6,5,4] <= 0
 -infinity <= delta[6,5,5] <= 0
 -infinity <= delta[6,5,6] <= 0
 -infinity <= delta[6,5,7] <= 0
 -infinity <= delta[6,5,8] <= 0
 -infinity <= delta[6,5,9] <= 0
 -infinity <= delta[6,6,0] <= 0
 -infinity <= delta[6,6,1] <= 0
 -infinity <= delta[6,6,2] <= 0
 -infinity <= delta[6,6,3] <= 0
 -infinity <= delta[6,6,4] <= 0
 -infinity <= delta[6,6,5] <= 0
 -infinity <= delta[6,6,6] <= 0
 -infinity <= delta[6,6,7] <= 0
 -infinity <= delta[6,6,8] <= 0
 -infinity <= delta[6,6,9] <= 0
 -infinity <= delta[6,7,0] <= 0
 -infinity <= delta[6,7,1] <= 0
 -infinity <= delta[6,7,2] <= 0
 -infinity <= delta[6,7,3] <= 0
 -infinity <= delta[6,7,4] <= 0
 -infinity <= delta[6,7,5] <= 0
 -infinity <= delta[6,7,6] <= 0
 -infinity <= delta[6,7,7] <= 0
 -infinity <= delta[6,7,8] <= 0
 -infinity <= delta[6,7,9] <= 0
 -infinity <= delta[6,8,0] <= 0
 -infinity <= delta[6,8,1] <= 0
 -infinity <= delta[6,8,2] <= 0
 -infinity <= delta[6,8,3] <= 0
 -infinity <= delta[6,8,4] <= 0
 -infinity <= delta[6,8,5] <= 0
 -infinity <= delta[6,8,6] <= 0
 -infinity <= delta[6,8,7] <= 0
 -infinity <= delta[6,8,8] <= 0
 -infinity <= delta[6,8,9] <= 0
 -infinity <= delta[6,9,0] <= 0
 -infinity <= delta[6,9,1] <= 0
 -infinity <= delta[6,9,2] <= 0
 -infinity <= delta[6,9,3] <= 0
 -infinity <= delta[6,9,4] <= 0
 -infinity <= delta[6,9,5] <= 0
 -infinity <= delta[6,9,6] <= 0
 -infinity <= delta[6,9,7] <= 0
 -infinity <= delta[6,9,8] <= 0
 -infinity <= delta[6,9,9] <= 0
 -infinity <= delta[7,0,0] <= 0
 -infinity <= delta[7,0,1] <= 0
 -infinity <= delta[7,0,2] <= 0
 -infinity <= delta[7,0,3] <= 0
 -infinity <= delta[7,0,4] <= 0
 -infinity <= delta[7,0,5] <= 0
 -infinity <= delta[7,0,6] <= 0
 -infinity <= delta[7,0,7] <= 0
 -infinity <= delta[7,0,8] <= 0
 -infinity <= delta[7,0,9] <= 0
 -infinity <= delta[7,1,0] <= 0
 -infinity <= delta[7,1,1] <= 0
 -infinity <= delta[7,1,2] <= 0
 -infinity <= delta[7,1,3] <= 0
 -infinity <= delta[7,1,4] <= 0
 -infinity <= delta[7,1,5] <= 0
 -infinity <= delta[7,1,6] <= 0
 -infinity <= delta[7,1,7] <= 0
 -infinity <= delta[7,1,8] <= 0
 -infinity <= delta[7,1,9] <= 0
 -infinity <= delta[7,2,0] <= 0
 -infinity <= delta[7,2,1] <= 0
 -infinity <= delta[7,2,2] <= 0
 -infinity <= delta[7,2,3] <= 0
 -infinity <= delta[7,2,4] <= 0
 -infinity <= delta[7,2,5] <= 0
 -infinity <= delta[7,2,6] <= 0
 -infinity <= delta[7,2,7] <= 0
 -infinity <= delta[7,2,8] <= 0
 -infinity <= delta[7,2,9] <= 0
 -infinity <= delta[7,3,0] <= 0
 -infinity <= delta[7,3,1] <= 0
 -infinity <= delta[7,3,2] <= 0
 -infinity <= delta[7,3,3] <= 0
 -infinity <= delta[7,3,4] <= 0
 -infinity <= delta[7,3,5] <= 0
 -infinity <= delta[7,3,6] <= 0
 -infinity <= delta[7,3,7] <= 0
 -infinity <= delta[7,3,8] <= 0
 -infinity <= delta[7,3,9] <= 0
 -infinity <= delta[7,4,0] <= 0
 -infinity <= delta[7,4,1] <= 0
 -infinity <= delta[7,4,2] <= 0
 -infinity <= delta[7,4,3] <= 0
 -infinity <= delta[7,4,4] <= 0
 -infinity <= delta[7,4,5] <= 0
 -infinity <= delta[7,4,6] <= 0
 -infinity <= delta[7,4,7] <= 0
 -infinity <= delta[7,4,8] <= 0
 -infinity <= delta[7,4,9] <= 0
 -infinity <= delta[7,5,0] <= 0
 -infinity <= delta[7,5,1] <= 0
 -infinity <= delta[7,5,2] <= 0
 -infinity <= delta[7,5,3] <= 0
 -infinity <= delta[7,5,4] <= 0
 -infinity <= delta[7,5,5] <= 0
 -infinity <= delta[7,5,6] <= 0
 -infinity <= delta[7,5,7] <= 0
 -infinity <= delta[7,5,8] <= 0
 -infinity <= delta[7,5,9] <= 0
 -infinity <= delta[7,6,0] <= 0
 -infinity <= delta[7,6,1] <= 0
 -infinity <= delta[7,6,2] <= 0
 -infinity <= delta[7,6,3] <= 0
 -infinity <= delta[7,6,4] <= 0
 -infinity <= delta[7,6,5] <= 0
 -infinity <= delta[7,6,6] <= 0
 -infinity <= delta[7,6,7] <= 0
 -infinity <= delta[7,6,8] <= 0
 -infinity <= delta[7,6,9] <= 0
 -infinity <= delta[7,7,0] <= 0
 -infinity <= delta[7,7,1] <= 0
 -infinity <= delta[7,7,2] <= 0
 -infinity <= delta[7,7,3] <= 0
 -infinity <= delta[7,7,4] <= 0
 -infinity <= delta[7,7,5] <= 0
 -infinity <= delta[7,7,6] <= 0
 -infinity <= delta[7,7,7] <= 0
 -infinity <= delta[7,7,8] <= 0
 -infinity <= delta[7,7,9] <= 0
 -infinity <= delta[7,8,0] <= 0
 -infinity <= delta[7,8,1] <= 0
 -infinity <= delta[7,8,2] <= 0
 -infinity <= delta[7,8,3] <= 0
 -infinity <= delta[7,8,4] <= 0
 -infinity <= delta[7,8,5] <= 0
 -infinity <= delta[7,8,6] <= 0
 -infinity <= delta[7,8,7] <= 0
 -infinity <= delta[7,8,8] <= 0
 -infinity <= delta[7,8,9] <= 0
 -infinity <= delta[7,9,0] <= 0
 -infinity <= delta[7,9,1] <= 0
 -infinity <= delta[7,9,2] <= 0
 -infinity <= delta[7,9,3] <= 0
 -infinity <= delta[7,9,4] <= 0
 -infinity <= delta[7,9,5] <= 0
 -infinity <= delta[7,9,6] <= 0
 -infinity <= delta[7,9,7] <= 0
 -infinity <= delta[7,9,8] <= 0
 -infinity <= delta[7,9,9] <= 0
 -infinity <= delta[8,0,0] <= 0
 -infinity <= delta[8,0,1] <= 0
 -infinity <= delta[8,0,2] <= 0
 -infinity <= delta[8,0,3] <= 0
 -infinity <= delta[8,0,4] <= 0
 -infinity <= delta[8,0,5] <= 0
 -infinity <= delta[8,0,6] <= 0
 -infinity <= delta[8,0,7] <= 0
 -infinity <= delta[8,0,8] <= 0
 -infinity <= delta[8,0,9] <= 0
 -infinity <= delta[8,1,0] <= 0
 -infinity <= delta[8,1,1] <= 0
 -infinity <= delta[8,1,2] <= 0
 -infinity <= delta[8,1,3] <= 0
 -infinity <= delta[8,1,4] <= 0
 -infinity <= delta[8,1,5] <= 0
 -infinity <= delta[8,1,6] <= 0
 -infinity <= delta[8,1,7] <= 0
 -infinity <= delta[8,1,8] <= 0
 -infinity <= delta[8,1,9] <= 0
 -infinity <= delta[8,2,0] <= 0
 -infinity <= delta[8,2,1] <= 0
 -infinity <= delta[8,2,2] <= 0
 -infinity <= delta[8,2,3] <= 0
 -infinity <= delta[8,2,4] <= 0
 -infinity <= delta[8,2,5] <= 0
 -infinity <= delta[8,2,6] <= 0
 -infinity <= delta[8,2,7] <= 0
 -infinity <= delta[8,2,8] <= 0
 -infinity <= delta[8,2,9] <= 0
 -infinity <= delta[8,3,0] <= 0
 -infinity <= delta[8,3,1] <= 0
 -infinity <= delta[8,3,2] <= 0
 -infinity <= delta[8,3,3] <= 0
 -infinity <= delta[8,3,4] <= 0
 -infinity <= delta[8,3,5] <= 0
 -infinity <= delta[8,3,6] <= 0
 -infinity <= delta[8,3,7] <= 0
 -infinity <= delta[8,3,8] <= 0
 -infinity <= delta[8,3,9] <= 0
 -infinity <= delta[8,4,0] <= 0
 -infinity <= delta[8,4,1] <= 0
 -infinity <= delta[8,4,2] <= 0
 -infinity <= delta[8,4,3] <= 0
 -infinity <= delta[8,4,4] <= 0
 -infinity <= delta[8,4,5] <= 0
 -infinity <= delta[8,4,6] <= 0
 -infinity <= delta[8,4,7] <= 0
 -infinity <= delta[8,4,8] <= 0
 -infinity <= delta[8,4,9] <= 0
 -infinity <= delta[8,5,0] <= 0
 -infinity <= delta[8,5,1] <= 0
 -infinity <= delta[8,5,2] <= 0
 -infinity <= delta[8,5,3] <= 0
 -infinity <= delta[8,5,4] <= 0
 -infinity <= delta[8,5,5] <= 0
 -infinity <= delta[8,5,6] <= 0
 -infinity <= delta[8,5,7] <= 0
 -infinity <= delta[8,5,8] <= 0
 -infinity <= delta[8,5,9] <= 0
 -infinity <= delta[8,6,0] <= 0
 -infinity <= delta[8,6,1] <= 0
 -infinity <= delta[8,6,2] <= 0
 -infinity <= delta[8,6,3] <= 0
 -infinity <= delta[8,6,4] <= 0
 -infinity <= delta[8,6,5] <= 0
 -infinity <= delta[8,6,6] <= 0
 -infinity <= delta[8,6,7] <= 0
 -infinity <= delta[8,6,8] <= 0
 -infinity <= delta[8,6,9] <= 0
 -infinity <= delta[8,7,0] <= 0
 -infinity <= delta[8,7,1] <= 0
 -infinity <= delta[8,7,2] <= 0
 -infinity <= delta[8,7,3] <= 0
 -infinity <= delta[8,7,4] <= 0
 -infinity <= delta[8,7,5] <= 0
 -infinity <= delta[8,7,6] <= 0
 -infinity <= delta[8,7,7] <= 0
 -infinity <= delta[8,7,8] <= 0
 -infinity <= delta[8,7,9] <= 0
 -infinity <= delta[8,8,0] <= 0
 -infinity <= delta[8,8,1] <= 0
 -infinity <= delta[8,8,2] <= 0
 -infinity <= delta[8,8,3] <= 0
 -infinity <= delta[8,8,4] <= 0
 -infinity <= delta[8,8,5] <= 0
 -infinity <= delta[8,8,6] <= 0
 -infinity <= delta[8,8,7] <= 0
 -infinity <= delta[8,8,8] <= 0
 -infinity <= delta[8,8,9] <= 0
 -infinity <= delta[8,9,0] <= 0
 -infinity <= delta[8,9,1] <= 0
 -infinity <= delta[8,9,2] <= 0
 -infinity <= delta[8,9,3] <= 0
 -infinity <= delta[8,9,4] <= 0
 -infinity <= delta[8,9,5] <= 0
 -infinity <= delta[8,9,6] <= 0
 -infinity <= delta[8,9,7] <= 0
 -infinity <= delta[8,9,8] <= 0
 -infinity <= delta[8,9,9] <= 0
 -infinity <= delta[9,0,0] <= 0
 -infinity <= delta[9,0,1] <= 0
 -infinity <= delta[9,0,2] <= 0
 -infinity <= delta[9,0,3] <= 0
 -infinity <= delta[9,0,4] <= 0
 -infinity <= delta[9,0,5] <= 0
 -infinity <= delta[9,0,6] <= 0
 -infinity <= delta[9,0,7] <= 0
 -infinity <= delta[9,0,8] <= 0
 -infinity <= delta[9,0,9] <= 0
 -infinity <= delta[9,1,0] <= 0
 -infinity <= delta[9,1,1] <= 0
 -infinity <= delta[9,1,2] <= 0
 -infinity <= delta[9,1,3] <= 0
 -infinity <= delta[9,1,4] <= 0
 -infinity <= delta[9,1,5] <= 0
 -infinity <= delta[9,1,6] <= 0
 -infinity <= delta[9,1,7] <= 0
 -infinity <= delta[9,1,8] <= 0
 -infinity <= delta[9,1,9] <= 0
 -infinity <= delta[9,2,0] <= 0
 -infinity <= delta[9,2,1] <= 0
 -infinity <= delta[9,2,2] <= 0
 -infinity <= delta[9,2,3] <= 0
 -infinity <= delta[9,2,4] <= 0
 -infinity <= delta[9,2,5] <= 0
 -infinity <= delta[9,2,6] <= 0
 -infinity <= delta[9,2,7] <= 0
 -infinity <= delta[9,2,8] <= 0
 -infinity <= delta[9,2,9] <= 0
 -infinity <= delta[9,3,0] <= 0
 -infinity <= delta[9,3,1] <= 0
 -infinity <= delta[9,3,2] <= 0
 -infinity <= delta[9,3,3] <= 0
 -infinity <= delta[9,3,4] <= 0
 -infinity <= delta[9,3,5] <= 0
 -infinity <= delta[9,3,6] <= 0
 -infinity <= delta[9,3,7] <= 0
 -infinity <= delta[9,3,8] <= 0
 -infinity <= delta[9,3,9] <= 0
 -infinity <= delta[9,4,0] <= 0
 -infinity <= delta[9,4,1] <= 0
 -infinity <= delta[9,4,2] <= 0
 -infinity <= delta[9,4,3] <= 0
 -infinity <= delta[9,4,4] <= 0
 -infinity <= delta[9,4,5] <= 0
 -infinity <= delta[9,4,6] <= 0
 -infinity <= delta[9,4,7] <= 0
 -infinity <= delta[9,4,8] <= 0
 -infinity <= delta[9,4,9] <= 0
 -infinity <= delta[9,5,0] <= 0
 -infinity <= delta[9,5,1] <= 0
 -infinity <= delta[9,5,2] <= 0
 -infinity <= delta[9,5,3] <= 0
 -infinity <= delta[9,5,4] <= 0
 -infinity <= delta[9,5,5] <= 0
 -infinity <= delta[9,5,6] <= 0
 -infinity <= delta[9,5,7] <= 0
 -infinity <= delta[9,5,8] <= 0
 -infinity <= delta[9,5,9] <= 0
 -infinity <= delta[9,6,0] <= 0
 -infinity <= delta[9,6,1] <= 0
 -infinity <= delta[9,6,2] <= 0
 -infinity <= delta[9,6,3] <= 0
 -infinity <= delta[9,6,4] <= 0
 -infinity <= delta[9,6,5] <= 0
 -infinity <= delta[9,6,6] <= 0
 -infinity <= delta[9,6,7] <= 0
 -infinity <= delta[9,6,8] <= 0
 -infinity <= delta[9,6,9] <= 0
 -infinity <= delta[9,7,0] <= 0
 -infinity <= delta[9,7,1] <= 0
 -infinity <= delta[9,7,2] <= 0
 -infinity <= delta[9,7,3] <= 0
 -infinity <= delta[9,7,4] <= 0
 -infinity <= delta[9,7,5] <= 0
 -infinity <= delta[9,7,6] <= 0
 -infinity <= delta[9,7,7] <= 0
 -infinity <= delta[9,7,8] <= 0
 -infinity <= delta[9,7,9] <= 0
 -infinity <= delta[9,8,0] <= 0
 -infinity <= delta[9,8,1] <= 0
 -infinity <= delta[9,8,2] <= 0
 -infinity <= delta[9,8,3] <= 0
 -infinity <= delta[9,8,4] <= 0
 -infinity <= delta[9,8,5] <= 0
 -infinity <= delta[9,8,6] <= 0
 -infinity <= delta[9,8,7] <= 0
 -infinity <= delta[9,8,8] <= 0
 -infinity <= delta[9,8,9] <= 0
 -infinity <= delta[9,9,0] <= 0
 -infinity <= delta[9,9,1] <= 0
 -infinity <= delta[9,9,2] <= 0
 -infinity <= delta[9,9,3] <= 0
 -infinity <= delta[9,9,4] <= 0
 -infinity <= delta[9,9,5] <= 0
 -infinity <= delta[9,9,6] <= 0
 -infinity <= delta[9,9,7] <= 0
 -infinity <= delta[9,9,8] <= 0
 -infinity <= delta[9,9,9] <= 0
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
 -infinity <= lambda[5,0] <= 0
 -infinity <= lambda[5,1] <= 0
 -infinity <= lambda[5,2] <= 0
 -infinity <= lambda[5,3] <= 0
 -infinity <= lambda[5,4] <= 0
 -infinity <= lambda[5,5] <= 0
 -infinity <= lambda[5,6] <= 0
 -infinity <= lambda[5,7] <= 0
 -infinity <= lambda[5,8] <= 0
 -infinity <= lambda[5,9] <= 0
 -infinity <= lambda[6,0] <= 0
 -infinity <= lambda[6,1] <= 0
 -infinity <= lambda[6,2] <= 0
 -infinity <= lambda[6,3] <= 0
 -infinity <= lambda[6,4] <= 0
 -infinity <= lambda[6,5] <= 0
 -infinity <= lambda[6,6] <= 0
 -infinity <= lambda[6,7] <= 0
 -infinity <= lambda[6,8] <= 0
 -infinity <= lambda[6,9] <= 0
 -infinity <= lambda[7,0] <= 0
 -infinity <= lambda[7,1] <= 0
 -infinity <= lambda[7,2] <= 0
 -infinity <= lambda[7,3] <= 0
 -infinity <= lambda[7,4] <= 0
 -infinity <= lambda[7,5] <= 0
 -infinity <= lambda[7,6] <= 0
 -infinity <= lambda[7,7] <= 0
 -infinity <= lambda[7,8] <= 0
 -infinity <= lambda[7,9] <= 0
 -infinity <= lambda[8,0] <= 0
 -infinity <= lambda[8,1] <= 0
 -infinity <= lambda[8,2] <= 0
 -infinity <= lambda[8,3] <= 0
 -infinity <= lambda[8,4] <= 0
 -infinity <= lambda[8,5] <= 0
 -infinity <= lambda[8,6] <= 0
 -infinity <= lambda[8,7] <= 0
 -infinity <= lambda[8,8] <= 0
 -infinity <= lambda[8,9] <= 0
 -infinity <= lambda[9,0] <= 0
 -infinity <= lambda[9,1] <= 0
 -infinity <= lambda[9,2] <= 0
 -infinity <= lambda[9,3] <= 0
 -infinity <= lambda[9,4] <= 0
 -infinity <= lambda[9,5] <= 0
 -infinity <= lambda[9,6] <= 0
 -infinity <= lambda[9,7] <= 0
 -infinity <= lambda[9,8] <= 0
 -infinity <= lambda[9,9] <= 0
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
 -infinity <= mu[5,0] <= 0
 -infinity <= mu[5,1] <= 0
 -infinity <= mu[5,2] <= 0
 -infinity <= mu[5,3] <= 0
 -infinity <= mu[5,4] <= 0
 -infinity <= mu[5,5] <= 0
 -infinity <= mu[5,6] <= 0
 -infinity <= mu[5,7] <= 0
 -infinity <= mu[5,8] <= 0
 -infinity <= mu[5,9] <= 0
 -infinity <= mu[6,0] <= 0
 -infinity <= mu[6,1] <= 0
 -infinity <= mu[6,2] <= 0
 -infinity <= mu[6,3] <= 0
 -infinity <= mu[6,4] <= 0
 -infinity <= mu[6,5] <= 0
 -infinity <= mu[6,6] <= 0
 -infinity <= mu[6,7] <= 0
 -infinity <= mu[6,8] <= 0
 -infinity <= mu[6,9] <= 0
 -infinity <= mu[7,0] <= 0
 -infinity <= mu[7,1] <= 0
 -infinity <= mu[7,2] <= 0
 -infinity <= mu[7,3] <= 0
 -infinity <= mu[7,4] <= 0
 -infinity <= mu[7,5] <= 0
 -infinity <= mu[7,6] <= 0
 -infinity <= mu[7,7] <= 0
 -infinity <= mu[7,8] <= 0
 -infinity <= mu[7,9] <= 0
 -infinity <= mu[8,0] <= 0
 -infinity <= mu[8,1] <= 0
 -infinity <= mu[8,2] <= 0
 -infinity <= mu[8,3] <= 0
 -infinity <= mu[8,4] <= 0
 -infinity <= mu[8,5] <= 0
 -infinity <= mu[8,6] <= 0
 -infinity <= mu[8,7] <= 0
 -infinity <= mu[8,8] <= 0
 -infinity <= mu[8,9] <= 0
 -infinity <= mu[9,0] <= 0
 -infinity <= mu[9,1] <= 0
 -infinity <= mu[9,2] <= 0
 -infinity <= mu[9,3] <= 0
 -infinity <= mu[9,4] <= 0
 -infinity <= mu[9,5] <= 0
 -infinity <= mu[9,6] <= 0
 -infinity <= mu[9,7] <= 0
 -infinity <= mu[9,8] <= 0
 -infinity <= mu[9,9] <= 0
End
