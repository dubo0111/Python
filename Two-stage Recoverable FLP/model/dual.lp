\ Model sub model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Qk[0] + Qk[1] + Qk[2] + Qk[3] + Qk[4] + Qk[5] + Qk[6] + Qk[7] + Qk[8]
   + Qk[9]
Subject To
 Q(k)[0]: - gamma[0,0,0] - gamma[0,0,1] - gamma[0,0,9] - gamma[0,1,0]
   - gamma[0,1,1] - gamma[0,1,9] - gamma[0,2,0] - gamma[0,2,1]
   - gamma[0,2,9] - gamma[0,3,0] - gamma[0,3,1] - gamma[0,3,9]
   - gamma[0,4,0] - gamma[0,4,1] - gamma[0,4,9] - gamma[0,5,0]
   - gamma[0,5,1] - gamma[0,5,9] - gamma[0,6,0] - gamma[0,6,1]
   - gamma[0,6,9] - gamma[0,7,0] - gamma[0,7,1] - gamma[0,7,9]
   - gamma[0,8,0] - gamma[0,8,1] - gamma[0,8,9] - gamma[0,9,0]
   - gamma[0,9,1] - gamma[0,9,9] - delta[0,0,1] - delta[0,0,2]
   - delta[0,0,3] - delta[0,0,5] - delta[0,0,6] - delta[0,0,7]
   - delta[0,0,8] - delta[0,0,9] - delta[0,1,1] - delta[0,1,2]
   - delta[0,1,3] - delta[0,1,5] - delta[0,1,6] - delta[0,1,7]
   - delta[0,1,8] - delta[0,1,9] - delta[0,2,1] - delta[0,2,2]
   - delta[0,2,3] - delta[0,2,5] - delta[0,2,6] - delta[0,2,7]
   - delta[0,2,8] - delta[0,2,9] - delta[0,3,1] - delta[0,3,2]
   - delta[0,3,3] - delta[0,3,5] - delta[0,3,6] - delta[0,3,7]
   - delta[0,3,8] - delta[0,3,9] - delta[0,4,1] - delta[0,4,2]
   - delta[0,4,3] - delta[0,4,5] - delta[0,4,6] - delta[0,4,7]
   - delta[0,4,8] - delta[0,4,9] - delta[0,5,1] - delta[0,5,2]
   - delta[0,5,3] - delta[0,5,5] - delta[0,5,6] - delta[0,5,7]
   - delta[0,5,8] - delta[0,5,9] - delta[0,6,1] - delta[0,6,2]
   - delta[0,6,3] - delta[0,6,5] - delta[0,6,6] - delta[0,6,7]
   - delta[0,6,8] - delta[0,6,9] - delta[0,7,1] - delta[0,7,2]
   - delta[0,7,3] - delta[0,7,5] - delta[0,7,6] - delta[0,7,7]
   - delta[0,7,8] - delta[0,7,9] - delta[0,8,1] - delta[0,8,2]
   - delta[0,8,3] - delta[0,8,5] - delta[0,8,6] - delta[0,8,7]
   - delta[0,8,8] - delta[0,8,9] - delta[0,9,1] - delta[0,9,2]
   - delta[0,9,3] - delta[0,9,5] - delta[0,9,6] - delta[0,9,7]
   - delta[0,9,8] - delta[0,9,9] - epsilon[0,0] - epsilon[0,1]
   - epsilon[0,2] - epsilon[0,3] - epsilon[0,4] - epsilon[0,5]
   - epsilon[0,6] - epsilon[0,7] - epsilon[0,8] - epsilon[0,9]
   - lambda[0,2] - lambda[0,3] - lambda[0,4] - lambda[0,5] - lambda[0,6]
   - lambda[0,7] - lambda[0,8] - mu[0,1] - mu[0,2] - mu[0,3] - mu[0,5]
   - mu[0,6] - mu[0,7] - mu[0,8] - mu[0,9] - nu[0] + Qk[0] = 0
 Q(k)[1]: - gamma[1,0,0] - gamma[1,0,1] - gamma[1,0,9] - gamma[1,1,0]
   - gamma[1,1,1] - gamma[1,1,9] - gamma[1,2,0] - gamma[1,2,1]
   - gamma[1,2,9] - gamma[1,3,0] - gamma[1,3,1] - gamma[1,3,9]
   - gamma[1,4,0] - gamma[1,4,1] - gamma[1,4,9] - gamma[1,5,0]
   - gamma[1,5,1] - gamma[1,5,9] - gamma[1,6,0] - gamma[1,6,1]
   - gamma[1,6,9] - gamma[1,7,0] - gamma[1,7,1] - gamma[1,7,9]
   - gamma[1,8,0] - gamma[1,8,1] - gamma[1,8,9] - gamma[1,9,0]
   - gamma[1,9,1] - gamma[1,9,9] - delta[1,0,0] - delta[1,0,2]
   - delta[1,0,4] - delta[1,0,5] - delta[1,0,6] - delta[1,0,7]
   - delta[1,0,8] - delta[1,0,9] - delta[1,1,0] - delta[1,1,2]
   - delta[1,1,4] - delta[1,1,5] - delta[1,1,6] - delta[1,1,7]
   - delta[1,1,8] - delta[1,1,9] - delta[1,2,0] - delta[1,2,2]
   - delta[1,2,4] - delta[1,2,5] - delta[1,2,6] - delta[1,2,7]
   - delta[1,2,8] - delta[1,2,9] - delta[1,3,0] - delta[1,3,2]
   - delta[1,3,4] - delta[1,3,5] - delta[1,3,6] - delta[1,3,7]
   - delta[1,3,8] - delta[1,3,9] - delta[1,4,0] - delta[1,4,2]
   - delta[1,4,4] - delta[1,4,5] - delta[1,4,6] - delta[1,4,7]
   - delta[1,4,8] - delta[1,4,9] - delta[1,5,0] - delta[1,5,2]
   - delta[1,5,4] - delta[1,5,5] - delta[1,5,6] - delta[1,5,7]
   - delta[1,5,8] - delta[1,5,9] - delta[1,6,0] - delta[1,6,2]
   - delta[1,6,4] - delta[1,6,5] - delta[1,6,6] - delta[1,6,7]
   - delta[1,6,8] - delta[1,6,9] - delta[1,7,0] - delta[1,7,2]
   - delta[1,7,4] - delta[1,7,5] - delta[1,7,6] - delta[1,7,7]
   - delta[1,7,8] - delta[1,7,9] - delta[1,8,0] - delta[1,8,2]
   - delta[1,8,4] - delta[1,8,5] - delta[1,8,6] - delta[1,8,7]
   - delta[1,8,8] - delta[1,8,9] - delta[1,9,0] - delta[1,9,2]
   - delta[1,9,4] - delta[1,9,5] - delta[1,9,6] - delta[1,9,7]
   - delta[1,9,8] - delta[1,9,9] - epsilon[1,0] - epsilon[1,1]
   - epsilon[1,2] - epsilon[1,3] - epsilon[1,4] - epsilon[1,5]
   - epsilon[1,6] - epsilon[1,7] - epsilon[1,8] - epsilon[1,9]
   - lambda[1,2] - lambda[1,3] - lambda[1,4] - lambda[1,5] - lambda[1,6]
   - lambda[1,7] - lambda[1,8] - mu[1,0] - mu[1,2] - mu[1,4] - mu[1,5]
   - mu[1,6] - mu[1,7] - mu[1,8] - mu[1,9] - nu[1] + Qk[1] = 0
 Q(k)[2]: - gamma[2,0,0] - gamma[2,0,1] - gamma[2,0,9] - gamma[2,1,0]
   - gamma[2,1,1] - gamma[2,1,9] - gamma[2,2,0] - gamma[2,2,1]
   - gamma[2,2,9] - gamma[2,3,0] - gamma[2,3,1] - gamma[2,3,9]
   - gamma[2,4,0] - gamma[2,4,1] - gamma[2,4,9] - gamma[2,5,0]
   - gamma[2,5,1] - gamma[2,5,9] - gamma[2,6,0] - gamma[2,6,1]
   - gamma[2,6,9] - gamma[2,7,0] - gamma[2,7,1] - gamma[2,7,9]
   - gamma[2,8,0] - gamma[2,8,1] - gamma[2,8,9] - gamma[2,9,0]
   - gamma[2,9,1] - gamma[2,9,9] - delta[2,0,0] - delta[2,0,1]
   - delta[2,0,2] - delta[2,0,3] - delta[2,0,4] - delta[2,0,6]
   - delta[2,0,7] - delta[2,0,9] - delta[2,1,0] - delta[2,1,1]
   - delta[2,1,2] - delta[2,1,3] - delta[2,1,4] - delta[2,1,6]
   - delta[2,1,7] - delta[2,1,9] - delta[2,2,0] - delta[2,2,1]
   - delta[2,2,2] - delta[2,2,3] - delta[2,2,4] - delta[2,2,6]
   - delta[2,2,7] - delta[2,2,9] - delta[2,3,0] - delta[2,3,1]
   - delta[2,3,2] - delta[2,3,3] - delta[2,3,4] - delta[2,3,6]
   - delta[2,3,7] - delta[2,3,9] - delta[2,4,0] - delta[2,4,1]
   - delta[2,4,2] - delta[2,4,3] - delta[2,4,4] - delta[2,4,6]
   - delta[2,4,7] - delta[2,4,9] - delta[2,5,0] - delta[2,5,1]
   - delta[2,5,2] - delta[2,5,3] - delta[2,5,4] - delta[2,5,6]
   - delta[2,5,7] - delta[2,5,9] - delta[2,6,0] - delta[2,6,1]
   - delta[2,6,2] - delta[2,6,3] - delta[2,6,4] - delta[2,6,6]
   - delta[2,6,7] - delta[2,6,9] - delta[2,7,0] - delta[2,7,1]
   - delta[2,7,2] - delta[2,7,3] - delta[2,7,4] - delta[2,7,6]
   - delta[2,7,7] - delta[2,7,9] - delta[2,8,0] - delta[2,8,1]
   - delta[2,8,2] - delta[2,8,3] - delta[2,8,4] - delta[2,8,6]
   - delta[2,8,7] - delta[2,8,9] - delta[2,9,0] - delta[2,9,1]
   - delta[2,9,2] - delta[2,9,3] - delta[2,9,4] - delta[2,9,6]
   - delta[2,9,7] - delta[2,9,9] - epsilon[2,0] - epsilon[2,1]
   - epsilon[2,2] - epsilon[2,3] - epsilon[2,4] - epsilon[2,5]
   - epsilon[2,6] - epsilon[2,7] - epsilon[2,8] - epsilon[2,9]
   - lambda[2,2] - lambda[2,3] - lambda[2,4] - lambda[2,5] - lambda[2,6]
   - lambda[2,7] - lambda[2,8] - mu[2,0] - mu[2,1] - mu[2,2] - mu[2,3]
   - mu[2,4] - mu[2,6] - mu[2,7] - mu[2,9] + Qk[2] = 0
 Q(k)[3]: - gamma[3,0,0] - gamma[3,0,1] - gamma[3,0,9] - gamma[3,1,0]
   - gamma[3,1,1] - gamma[3,1,9] - gamma[3,2,0] - gamma[3,2,1]
   - gamma[3,2,9] - gamma[3,3,0] - gamma[3,3,1] - gamma[3,3,9]
   - gamma[3,4,0] - gamma[3,4,1] - gamma[3,4,9] - gamma[3,5,0]
   - gamma[3,5,1] - gamma[3,5,9] - gamma[3,6,0] - gamma[3,6,1]
   - gamma[3,6,9] - gamma[3,7,0] - gamma[3,7,1] - gamma[3,7,9]
   - gamma[3,8,0] - gamma[3,8,1] - gamma[3,8,9] - gamma[3,9,0]
   - gamma[3,9,1] - gamma[3,9,9] - delta[3,0,0] - delta[3,0,1]
   - delta[3,0,2] - delta[3,0,3] - delta[3,0,4] - delta[3,0,5]
   - delta[3,0,6] - delta[3,0,9] - delta[3,1,0] - delta[3,1,1]
   - delta[3,1,2] - delta[3,1,3] - delta[3,1,4] - delta[3,1,5]
   - delta[3,1,6] - delta[3,1,9] - delta[3,2,0] - delta[3,2,1]
   - delta[3,2,2] - delta[3,2,3] - delta[3,2,4] - delta[3,2,5]
   - delta[3,2,6] - delta[3,2,9] - delta[3,3,0] - delta[3,3,1]
   - delta[3,3,2] - delta[3,3,3] - delta[3,3,4] - delta[3,3,5]
   - delta[3,3,6] - delta[3,3,9] - delta[3,4,0] - delta[3,4,1]
   - delta[3,4,2] - delta[3,4,3] - delta[3,4,4] - delta[3,4,5]
   - delta[3,4,6] - delta[3,4,9] - delta[3,5,0] - delta[3,5,1]
   - delta[3,5,2] - delta[3,5,3] - delta[3,5,4] - delta[3,5,5]
   - delta[3,5,6] - delta[3,5,9] - delta[3,6,0] - delta[3,6,1]
   - delta[3,6,2] - delta[3,6,3] - delta[3,6,4] - delta[3,6,5]
   - delta[3,6,6] - delta[3,6,9] - delta[3,7,0] - delta[3,7,1]
   - delta[3,7,2] - delta[3,7,3] - delta[3,7,4] - delta[3,7,5]
   - delta[3,7,6] - delta[3,7,9] - delta[3,8,0] - delta[3,8,1]
   - delta[3,8,2] - delta[3,8,3] - delta[3,8,4] - delta[3,8,5]
   - delta[3,8,6] - delta[3,8,9] - delta[3,9,0] - delta[3,9,1]
   - delta[3,9,2] - delta[3,9,3] - delta[3,9,4] - delta[3,9,5]
   - delta[3,9,6] - delta[3,9,9] - epsilon[3,0] - epsilon[3,1]
   - epsilon[3,2] - epsilon[3,3] - epsilon[3,4] - epsilon[3,5]
   - epsilon[3,6] - epsilon[3,7] - epsilon[3,8] - epsilon[3,9]
   - lambda[3,2] - lambda[3,3] - lambda[3,4] - lambda[3,5] - lambda[3,6]
   - lambda[3,7] - lambda[3,8] - mu[3,0] - mu[3,1] - mu[3,2] - mu[3,3]
   - mu[3,4] - mu[3,5] - mu[3,6] - mu[3,9] + Qk[3] = 0
 Q(k)[4]: - gamma[4,0,0] - gamma[4,0,1] - gamma[4,0,9] - gamma[4,1,0]
   - gamma[4,1,1] - gamma[4,1,9] - gamma[4,2,0] - gamma[4,2,1]
   - gamma[4,2,9] - gamma[4,3,0] - gamma[4,3,1] - gamma[4,3,9]
   - gamma[4,4,0] - gamma[4,4,1] - gamma[4,4,9] - gamma[4,5,0]
   - gamma[4,5,1] - gamma[4,5,9] - gamma[4,6,0] - gamma[4,6,1]
   - gamma[4,6,9] - gamma[4,7,0] - gamma[4,7,1] - gamma[4,7,9]
   - gamma[4,8,0] - gamma[4,8,1] - gamma[4,8,9] - gamma[4,9,0]
   - gamma[4,9,1] - gamma[4,9,9] - delta[4,0,0] - delta[4,0,1]
   - delta[4,0,3] - delta[4,0,4] - delta[4,0,5] - delta[4,0,6]
   - delta[4,0,7] - delta[4,0,9] - delta[4,1,0] - delta[4,1,1]
   - delta[4,1,3] - delta[4,1,4] - delta[4,1,5] - delta[4,1,6]
   - delta[4,1,7] - delta[4,1,9] - delta[4,2,0] - delta[4,2,1]
   - delta[4,2,3] - delta[4,2,4] - delta[4,2,5] - delta[4,2,6]
   - delta[4,2,7] - delta[4,2,9] - delta[4,3,0] - delta[4,3,1]
   - delta[4,3,3] - delta[4,3,4] - delta[4,3,5] - delta[4,3,6]
   - delta[4,3,7] - delta[4,3,9] - delta[4,4,0] - delta[4,4,1]
   - delta[4,4,3] - delta[4,4,4] - delta[4,4,5] - delta[4,4,6]
   - delta[4,4,7] - delta[4,4,9] - delta[4,5,0] - delta[4,5,1]
   - delta[4,5,3] - delta[4,5,4] - delta[4,5,5] - delta[4,5,6]
   - delta[4,5,7] - delta[4,5,9] - delta[4,6,0] - delta[4,6,1]
   - delta[4,6,3] - delta[4,6,4] - delta[4,6,5] - delta[4,6,6]
   - delta[4,6,7] - delta[4,6,9] - delta[4,7,0] - delta[4,7,1]
   - delta[4,7,3] - delta[4,7,4] - delta[4,7,5] - delta[4,7,6]
   - delta[4,7,7] - delta[4,7,9] - delta[4,8,0] - delta[4,8,1]
   - delta[4,8,3] - delta[4,8,4] - delta[4,8,5] - delta[4,8,6]
   - delta[4,8,7] - delta[4,8,9] - delta[4,9,0] - delta[4,9,1]
   - delta[4,9,3] - delta[4,9,4] - delta[4,9,5] - delta[4,9,6]
   - delta[4,9,7] - delta[4,9,9] - epsilon[4,0] - epsilon[4,1]
   - epsilon[4,2] - epsilon[4,3] - epsilon[4,4] - epsilon[4,5]
   - epsilon[4,6] - epsilon[4,7] - epsilon[4,8] - epsilon[4,9]
   - lambda[4,2] - lambda[4,3] - lambda[4,4] - lambda[4,5] - lambda[4,6]
   - lambda[4,7] - lambda[4,8] - mu[4,0] - mu[4,1] - mu[4,3] - mu[4,4]
   - mu[4,5] - mu[4,6] - mu[4,7] - mu[4,9] + Qk[4] = 0
 Q(k)[5]: - gamma[5,0,0] - gamma[5,0,1] - gamma[5,0,9] - gamma[5,1,0]
   - gamma[5,1,1] - gamma[5,1,9] - gamma[5,2,0] - gamma[5,2,1]
   - gamma[5,2,9] - gamma[5,3,0] - gamma[5,3,1] - gamma[5,3,9]
   - gamma[5,4,0] - gamma[5,4,1] - gamma[5,4,9] - gamma[5,5,0]
   - gamma[5,5,1] - gamma[5,5,9] - gamma[5,6,0] - gamma[5,6,1]
   - gamma[5,6,9] - gamma[5,7,0] - gamma[5,7,1] - gamma[5,7,9]
   - gamma[5,8,0] - gamma[5,8,1] - gamma[5,8,9] - gamma[5,9,0]
   - gamma[5,9,1] - gamma[5,9,9] - delta[5,0,0] - delta[5,0,2]
   - delta[5,0,3] - delta[5,0,4] - delta[5,0,5] - delta[5,0,6]
   - delta[5,0,7] - delta[5,0,8] - delta[5,1,0] - delta[5,1,2]
   - delta[5,1,3] - delta[5,1,4] - delta[5,1,5] - delta[5,1,6]
   - delta[5,1,7] - delta[5,1,8] - delta[5,2,0] - delta[5,2,2]
   - delta[5,2,3] - delta[5,2,4] - delta[5,2,5] - delta[5,2,6]
   - delta[5,2,7] - delta[5,2,8] - delta[5,3,0] - delta[5,3,2]
   - delta[5,3,3] - delta[5,3,4] - delta[5,3,5] - delta[5,3,6]
   - delta[5,3,7] - delta[5,3,8] - delta[5,4,0] - delta[5,4,2]
   - delta[5,4,3] - delta[5,4,4] - delta[5,4,5] - delta[5,4,6]
   - delta[5,4,7] - delta[5,4,8] - delta[5,5,0] - delta[5,5,2]
   - delta[5,5,3] - delta[5,5,4] - delta[5,5,5] - delta[5,5,6]
   - delta[5,5,7] - delta[5,5,8] - delta[5,6,0] - delta[5,6,2]
   - delta[5,6,3] - delta[5,6,4] - delta[5,6,5] - delta[5,6,6]
   - delta[5,6,7] - delta[5,6,8] - delta[5,7,0] - delta[5,7,2]
   - delta[5,7,3] - delta[5,7,4] - delta[5,7,5] - delta[5,7,6]
   - delta[5,7,7] - delta[5,7,8] - delta[5,8,0] - delta[5,8,2]
   - delta[5,8,3] - delta[5,8,4] - delta[5,8,5] - delta[5,8,6]
   - delta[5,8,7] - delta[5,8,8] - delta[5,9,0] - delta[5,9,2]
   - delta[5,9,3] - delta[5,9,4] - delta[5,9,5] - delta[5,9,6]
   - delta[5,9,7] - delta[5,9,8] - epsilon[5,0] - epsilon[5,1]
   - epsilon[5,2] - epsilon[5,3] - epsilon[5,4] - epsilon[5,5]
   - epsilon[5,6] - epsilon[5,7] - epsilon[5,8] - epsilon[5,9]
   - lambda[5,2] - lambda[5,3] - lambda[5,4] - lambda[5,5] - lambda[5,6]
   - lambda[5,7] - lambda[5,8] - mu[5,0] - mu[5,2] - mu[5,3] - mu[5,4]
   - mu[5,5] - mu[5,6] - mu[5,7] - mu[5,8] - 2 nu[5] + Qk[5] = 0
 Q(k)[6]: - gamma[6,0,0] - gamma[6,0,1] - gamma[6,0,9] - gamma[6,1,0]
   - gamma[6,1,1] - gamma[6,1,9] - gamma[6,2,0] - gamma[6,2,1]
   - gamma[6,2,9] - gamma[6,3,0] - gamma[6,3,1] - gamma[6,3,9]
   - gamma[6,4,0] - gamma[6,4,1] - gamma[6,4,9] - gamma[6,5,0]
   - gamma[6,5,1] - gamma[6,5,9] - gamma[6,6,0] - gamma[6,6,1]
   - gamma[6,6,9] - gamma[6,7,0] - gamma[6,7,1] - gamma[6,7,9]
   - gamma[6,8,0] - gamma[6,8,1] - gamma[6,8,9] - gamma[6,9,0]
   - gamma[6,9,1] - gamma[6,9,9] - delta[6,0,0] - delta[6,0,1]
   - delta[6,0,2] - delta[6,0,4] - delta[6,0,5] - delta[6,0,6]
   - delta[6,0,7] - delta[6,0,8] - delta[6,1,0] - delta[6,1,1]
   - delta[6,1,2] - delta[6,1,4] - delta[6,1,5] - delta[6,1,6]
   - delta[6,1,7] - delta[6,1,8] - delta[6,2,0] - delta[6,2,1]
   - delta[6,2,2] - delta[6,2,4] - delta[6,2,5] - delta[6,2,6]
   - delta[6,2,7] - delta[6,2,8] - delta[6,3,0] - delta[6,3,1]
   - delta[6,3,2] - delta[6,3,4] - delta[6,3,5] - delta[6,3,6]
   - delta[6,3,7] - delta[6,3,8] - delta[6,4,0] - delta[6,4,1]
   - delta[6,4,2] - delta[6,4,4] - delta[6,4,5] - delta[6,4,6]
   - delta[6,4,7] - delta[6,4,8] - delta[6,5,0] - delta[6,5,1]
   - delta[6,5,2] - delta[6,5,4] - delta[6,5,5] - delta[6,5,6]
   - delta[6,5,7] - delta[6,5,8] - delta[6,6,0] - delta[6,6,1]
   - delta[6,6,2] - delta[6,6,4] - delta[6,6,5] - delta[6,6,6]
   - delta[6,6,7] - delta[6,6,8] - delta[6,7,0] - delta[6,7,1]
   - delta[6,7,2] - delta[6,7,4] - delta[6,7,5] - delta[6,7,6]
   - delta[6,7,7] - delta[6,7,8] - delta[6,8,0] - delta[6,8,1]
   - delta[6,8,2] - delta[6,8,4] - delta[6,8,5] - delta[6,8,6]
   - delta[6,8,7] - delta[6,8,8] - delta[6,9,0] - delta[6,9,1]
   - delta[6,9,2] - delta[6,9,4] - delta[6,9,5] - delta[6,9,6]
   - delta[6,9,7] - delta[6,9,8] - epsilon[6,0] - epsilon[6,1]
   - epsilon[6,2] - epsilon[6,3] - epsilon[6,4] - epsilon[6,5]
   - epsilon[6,6] - epsilon[6,7] - epsilon[6,8] - epsilon[6,9]
   - lambda[6,2] - lambda[6,3] - lambda[6,4] - lambda[6,5] - lambda[6,6]
   - lambda[6,7] - lambda[6,8] - mu[6,0] - mu[6,1] - mu[6,2] - mu[6,4]
   - mu[6,5] - mu[6,6] - mu[6,7] - mu[6,8] - nu[6] + Qk[6] = 0
 Q(k)[7]: - gamma[7,0,0] - gamma[7,0,1] - gamma[7,0,9] - gamma[7,1,0]
   - gamma[7,1,1] - gamma[7,1,9] - gamma[7,2,0] - gamma[7,2,1]
   - gamma[7,2,9] - gamma[7,3,0] - gamma[7,3,1] - gamma[7,3,9]
   - gamma[7,4,0] - gamma[7,4,1] - gamma[7,4,9] - gamma[7,5,0]
   - gamma[7,5,1] - gamma[7,5,9] - gamma[7,6,0] - gamma[7,6,1]
   - gamma[7,6,9] - gamma[7,7,0] - gamma[7,7,1] - gamma[7,7,9]
   - gamma[7,8,0] - gamma[7,8,1] - gamma[7,8,9] - gamma[7,9,0]
   - gamma[7,9,1] - gamma[7,9,9] - delta[7,0,0] - delta[7,0,1]
   - delta[7,0,2] - delta[7,0,3] - delta[7,0,4] - delta[7,0,6]
   - delta[7,0,8] - delta[7,0,9] - delta[7,1,0] - delta[7,1,1]
   - delta[7,1,2] - delta[7,1,3] - delta[7,1,4] - delta[7,1,6]
   - delta[7,1,8] - delta[7,1,9] - delta[7,2,0] - delta[7,2,1]
   - delta[7,2,2] - delta[7,2,3] - delta[7,2,4] - delta[7,2,6]
   - delta[7,2,8] - delta[7,2,9] - delta[7,3,0] - delta[7,3,1]
   - delta[7,3,2] - delta[7,3,3] - delta[7,3,4] - delta[7,3,6]
   - delta[7,3,8] - delta[7,3,9] - delta[7,4,0] - delta[7,4,1]
   - delta[7,4,2] - delta[7,4,3] - delta[7,4,4] - delta[7,4,6]
   - delta[7,4,8] - delta[7,4,9] - delta[7,5,0] - delta[7,5,1]
   - delta[7,5,2] - delta[7,5,3] - delta[7,5,4] - delta[7,5,6]
   - delta[7,5,8] - delta[7,5,9] - delta[7,6,0] - delta[7,6,1]
   - delta[7,6,2] - delta[7,6,3] - delta[7,6,4] - delta[7,6,6]
   - delta[7,6,8] - delta[7,6,9] - delta[7,7,0] - delta[7,7,1]
   - delta[7,7,2] - delta[7,7,3] - delta[7,7,4] - delta[7,7,6]
   - delta[7,7,8] - delta[7,7,9] - delta[7,8,0] - delta[7,8,1]
   - delta[7,8,2] - delta[7,8,3] - delta[7,8,4] - delta[7,8,6]
   - delta[7,8,8] - delta[7,8,9] - delta[7,9,0] - delta[7,9,1]
   - delta[7,9,2] - delta[7,9,3] - delta[7,9,4] - delta[7,9,6]
   - delta[7,9,8] - delta[7,9,9] - epsilon[7,0] - epsilon[7,1]
   - epsilon[7,2] - epsilon[7,3] - epsilon[7,4] - epsilon[7,5]
   - epsilon[7,6] - epsilon[7,7] - epsilon[7,8] - epsilon[7,9]
   - lambda[7,2] - lambda[7,3] - lambda[7,4] - lambda[7,5] - lambda[7,6]
   - lambda[7,7] - lambda[7,8] - mu[7,0] - mu[7,1] - mu[7,2] - mu[7,3]
   - mu[7,4] - mu[7,6] - mu[7,8] - mu[7,9] + Qk[7] = 0
 Q(k)[8]: - gamma[8,0,0] - gamma[8,0,1] - gamma[8,0,9] - gamma[8,1,0]
   - gamma[8,1,1] - gamma[8,1,9] - gamma[8,2,0] - gamma[8,2,1]
   - gamma[8,2,9] - gamma[8,3,0] - gamma[8,3,1] - gamma[8,3,9]
   - gamma[8,4,0] - gamma[8,4,1] - gamma[8,4,9] - gamma[8,5,0]
   - gamma[8,5,1] - gamma[8,5,9] - gamma[8,6,0] - gamma[8,6,1]
   - gamma[8,6,9] - gamma[8,7,0] - gamma[8,7,1] - gamma[8,7,9]
   - gamma[8,8,0] - gamma[8,8,1] - gamma[8,8,9] - gamma[8,9,0]
   - gamma[8,9,1] - gamma[8,9,9] - delta[8,0,0] - delta[8,0,1]
   - delta[8,0,2] - delta[8,0,4] - delta[8,0,5] - delta[8,0,6]
   - delta[8,0,7] - delta[8,0,8] - delta[8,1,0] - delta[8,1,1]
   - delta[8,1,2] - delta[8,1,4] - delta[8,1,5] - delta[8,1,6]
   - delta[8,1,7] - delta[8,1,8] - delta[8,2,0] - delta[8,2,1]
   - delta[8,2,2] - delta[8,2,4] - delta[8,2,5] - delta[8,2,6]
   - delta[8,2,7] - delta[8,2,8] - delta[8,3,0] - delta[8,3,1]
   - delta[8,3,2] - delta[8,3,4] - delta[8,3,5] - delta[8,3,6]
   - delta[8,3,7] - delta[8,3,8] - delta[8,4,0] - delta[8,4,1]
   - delta[8,4,2] - delta[8,4,4] - delta[8,4,5] - delta[8,4,6]
   - delta[8,4,7] - delta[8,4,8] - delta[8,5,0] - delta[8,5,1]
   - delta[8,5,2] - delta[8,5,4] - delta[8,5,5] - delta[8,5,6]
   - delta[8,5,7] - delta[8,5,8] - delta[8,6,0] - delta[8,6,1]
   - delta[8,6,2] - delta[8,6,4] - delta[8,6,5] - delta[8,6,6]
   - delta[8,6,7] - delta[8,6,8] - delta[8,7,0] - delta[8,7,1]
   - delta[8,7,2] - delta[8,7,4] - delta[8,7,5] - delta[8,7,6]
   - delta[8,7,7] - delta[8,7,8] - delta[8,8,0] - delta[8,8,1]
   - delta[8,8,2] - delta[8,8,4] - delta[8,8,5] - delta[8,8,6]
   - delta[8,8,7] - delta[8,8,8] - delta[8,9,0] - delta[8,9,1]
   - delta[8,9,2] - delta[8,9,4] - delta[8,9,5] - delta[8,9,6]
   - delta[8,9,7] - delta[8,9,8] - epsilon[8,0] - epsilon[8,1]
   - epsilon[8,2] - epsilon[8,3] - epsilon[8,4] - epsilon[8,5]
   - epsilon[8,6] - epsilon[8,7] - epsilon[8,8] - epsilon[8,9]
   - lambda[8,2] - lambda[8,3] - lambda[8,4] - lambda[8,5] - lambda[8,6]
   - lambda[8,7] - lambda[8,8] - mu[8,0] - mu[8,1] - mu[8,2] - mu[8,4]
   - mu[8,5] - mu[8,6] - mu[8,7] - mu[8,8] - nu[8] + Qk[8] = 0
 Q(k)[9]: - gamma[9,0,0] - gamma[9,0,1] - gamma[9,0,9] - gamma[9,1,0]
   - gamma[9,1,1] - gamma[9,1,9] - gamma[9,2,0] - gamma[9,2,1]
   - gamma[9,2,9] - gamma[9,3,0] - gamma[9,3,1] - gamma[9,3,9]
   - gamma[9,4,0] - gamma[9,4,1] - gamma[9,4,9] - gamma[9,5,0]
   - gamma[9,5,1] - gamma[9,5,9] - gamma[9,6,0] - gamma[9,6,1]
   - gamma[9,6,9] - gamma[9,7,0] - gamma[9,7,1] - gamma[9,7,9]
   - gamma[9,8,0] - gamma[9,8,1] - gamma[9,8,9] - gamma[9,9,0]
   - gamma[9,9,1] - gamma[9,9,9] - delta[9,0,1] - delta[9,0,2]
   - delta[9,0,3] - delta[9,0,5] - delta[9,0,6] - delta[9,0,7]
   - delta[9,0,8] - delta[9,0,9] - delta[9,1,1] - delta[9,1,2]
   - delta[9,1,3] - delta[9,1,5] - delta[9,1,6] - delta[9,1,7]
   - delta[9,1,8] - delta[9,1,9] - delta[9,2,1] - delta[9,2,2]
   - delta[9,2,3] - delta[9,2,5] - delta[9,2,6] - delta[9,2,7]
   - delta[9,2,8] - delta[9,2,9] - delta[9,3,1] - delta[9,3,2]
   - delta[9,3,3] - delta[9,3,5] - delta[9,3,6] - delta[9,3,7]
   - delta[9,3,8] - delta[9,3,9] - delta[9,4,1] - delta[9,4,2]
   - delta[9,4,3] - delta[9,4,5] - delta[9,4,6] - delta[9,4,7]
   - delta[9,4,8] - delta[9,4,9] - delta[9,5,1] - delta[9,5,2]
   - delta[9,5,3] - delta[9,5,5] - delta[9,5,6] - delta[9,5,7]
   - delta[9,5,8] - delta[9,5,9] - delta[9,6,1] - delta[9,6,2]
   - delta[9,6,3] - delta[9,6,5] - delta[9,6,6] - delta[9,6,7]
   - delta[9,6,8] - delta[9,6,9] - delta[9,7,1] - delta[9,7,2]
   - delta[9,7,3] - delta[9,7,5] - delta[9,7,6] - delta[9,7,7]
   - delta[9,7,8] - delta[9,7,9] - delta[9,8,1] - delta[9,8,2]
   - delta[9,8,3] - delta[9,8,5] - delta[9,8,6] - delta[9,8,7]
   - delta[9,8,8] - delta[9,8,9] - delta[9,9,1] - delta[9,9,2]
   - delta[9,9,3] - delta[9,9,5] - delta[9,9,6] - delta[9,9,7]
   - delta[9,9,8] - delta[9,9,9] - epsilon[9,0] - epsilon[9,1]
   - epsilon[9,2] - epsilon[9,3] - epsilon[9,4] - epsilon[9,5]
   - epsilon[9,6] - epsilon[9,7] - epsilon[9,8] - epsilon[9,9]
   - lambda[9,2] - lambda[9,3] - lambda[9,4] - lambda[9,5] - lambda[9,6]
   - lambda[9,7] - lambda[9,8] - mu[9,1] - mu[9,2] - mu[9,3] - mu[9,5]
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
 v[0,0,1]: 232 beta[0,0] + gamma[0,0,1] + delta[0,0,1] + epsilon[0,0] <= 0
 v[0,0,2]: 282 beta[0,0] + gamma[0,0,2] + delta[0,0,2] + epsilon[0,0] <= 0
 v[0,0,3]: 256 beta[0,0] + gamma[0,0,3] + delta[0,0,3] + epsilon[0,0] <= 0
 v[0,0,4]: 330 beta[0,0] + gamma[0,0,4] + delta[0,0,4] + epsilon[0,0] <= 0
 v[0,0,5]: 945 beta[0,0] + gamma[0,0,5] + delta[0,0,5] + epsilon[0,0] <= 0
 v[0,0,6]: 798 beta[0,0] + gamma[0,0,6] + delta[0,0,6] + epsilon[0,0] <= 0
 v[0,0,7]: 843 beta[0,0] + gamma[0,0,7] + delta[0,0,7] + epsilon[0,0] <= 0
 v[0,0,8]: 1045 beta[0,0] + gamma[0,0,8] + delta[0,0,8] + epsilon[0,0]
   <= 0
 v[0,0,9]: 856 beta[0,0] + gamma[0,0,9] + delta[0,0,9] + epsilon[0,0] <= 0
 v[0,1,0]: 318 beta[0,1] + gamma[0,1,0] + delta[0,1,0] + epsilon[0,1] <= 0
 v[0,1,1]: gamma[0,1,1] + delta[0,1,1] + epsilon[0,1] <= 0
 v[0,1,2]: 568 beta[0,1] + gamma[0,1,2] + delta[0,1,2] + epsilon[0,1] <= 0
 v[0,1,3]: 760 beta[0,1] + gamma[0,1,3] + delta[0,1,3] + epsilon[0,1] <= 0
 v[0,1,4]: 395 beta[0,1] + gamma[0,1,4] + delta[0,1,4] + epsilon[0,1] <= 0
 v[0,1,5]: 710 beta[0,1] + gamma[0,1,5] + delta[0,1,5] + epsilon[0,1] <= 0
 v[0,1,6]: 203 beta[0,1] + gamma[0,1,6] + delta[0,1,6] + epsilon[0,1] <= 0
 v[0,1,7]: 603 beta[0,1] + gamma[0,1,7] + delta[0,1,7] + epsilon[0,1] <= 0
 v[0,1,8]: 727 beta[0,1] + gamma[0,1,8] + delta[0,1,8] + epsilon[0,1] <= 0
 v[0,1,9]: 513 beta[0,1] + gamma[0,1,9] + delta[0,1,9] + epsilon[0,1] <= 0
 v[0,2,0]: 779 beta[0,2] + gamma[0,2,0] + delta[0,2,0] + epsilon[0,2] <= 0
 v[0,2,1]: 704 beta[0,2] + gamma[0,2,1] + delta[0,2,1] + epsilon[0,2] <= 0
 v[0,2,2]: gamma[0,2,2] + delta[0,2,2] + epsilon[0,2] <= 0
 v[0,2,3]: 555 beta[0,2] + gamma[0,2,3] + delta[0,2,3] + epsilon[0,2] <= 0
 v[0,2,4]: 393 beta[0,2] + gamma[0,2,4] + delta[0,2,4] + epsilon[0,2] <= 0
 v[0,2,5]: 596 beta[0,2] + gamma[0,2,5] + delta[0,2,5] + epsilon[0,2] <= 0
 v[0,2,6]: 824 beta[0,2] + gamma[0,2,6] + delta[0,2,6] + epsilon[0,2] <= 0
 v[0,2,7]: 569 beta[0,2] + gamma[0,2,7] + delta[0,2,7] + epsilon[0,2] <= 0
 v[0,2,8]: 521 beta[0,2] + gamma[0,2,8] + delta[0,2,8] + epsilon[0,2] <= 0
 v[0,2,9]: 1084 beta[0,2] + gamma[0,2,9] + delta[0,2,9] + epsilon[0,2]
   <= 0
 v[0,3,0]: 507 beta[0,3] + gamma[0,3,0] + delta[0,3,0] + epsilon[0,3] <= 0
 v[0,3,1]: 231 beta[0,3] + gamma[0,3,1] + delta[0,3,1] + epsilon[0,3] <= 0
 v[0,3,2]: 450 beta[0,3] + gamma[0,3,2] + delta[0,3,2] + epsilon[0,3] <= 0
 v[0,3,3]: gamma[0,3,3] + delta[0,3,3] + epsilon[0,3] <= 0
 v[0,3,4]: 869 beta[0,3] + gamma[0,3,4] + delta[0,3,4] + epsilon[0,3] <= 0
 v[0,3,5]: 364 beta[0,3] + gamma[0,3,5] + delta[0,3,5] + epsilon[0,3] <= 0
 v[0,3,6]: 292 beta[0,3] + gamma[0,3,6] + delta[0,3,6] + epsilon[0,3] <= 0
 v[0,3,7]: 1080 beta[0,3] + gamma[0,3,7] + delta[0,3,7] + epsilon[0,3]
   <= 0
 v[0,3,8]: 918 beta[0,3] + gamma[0,3,8] + delta[0,3,8] + epsilon[0,3] <= 0
 v[0,3,9]: 436 beta[0,3] + gamma[0,3,9] + delta[0,3,9] + epsilon[0,3] <= 0
 v[0,4,0]: 859 beta[0,4] + gamma[0,4,0] + delta[0,4,0] + epsilon[0,4] <= 0
 v[0,4,1]: 597 beta[0,4] + gamma[0,4,1] + delta[0,4,1] + epsilon[0,4] <= 0
 v[0,4,2]: 886 beta[0,4] + gamma[0,4,2] + delta[0,4,2] + epsilon[0,4] <= 0
 v[0,4,3]: 137 beta[0,4] + gamma[0,4,3] + delta[0,4,3] + epsilon[0,4] <= 0
 v[0,4,4]: gamma[0,4,4] + delta[0,4,4] + epsilon[0,4] <= 0
 v[0,4,5]: 903 beta[0,4] + gamma[0,4,5] + delta[0,4,5] + epsilon[0,4] <= 0
 v[0,4,6]: 233 beta[0,4] + gamma[0,4,6] + delta[0,4,6] + epsilon[0,4] <= 0
 v[0,4,7]: 695 beta[0,4] + gamma[0,4,7] + delta[0,4,7] + epsilon[0,4] <= 0
 v[0,4,8]: 531 beta[0,4] + gamma[0,4,8] + delta[0,4,8] + epsilon[0,4] <= 0
 v[0,4,9]: 967 beta[0,4] + gamma[0,4,9] + delta[0,4,9] + epsilon[0,4] <= 0
 v[0,5,0]: 817 beta[0,5] + gamma[0,5,0] + delta[0,5,0] + epsilon[0,5] <= 0
 v[0,5,1]: 849 beta[0,5] + gamma[0,5,1] + delta[0,5,1] + epsilon[0,5] <= 0
 v[0,5,2]: 848 beta[0,5] + gamma[0,5,2] + delta[0,5,2] + epsilon[0,5] <= 0
 v[0,5,3]: 601 beta[0,5] + gamma[0,5,3] + delta[0,5,3] + epsilon[0,5] <= 0
 v[0,5,4]: 893 beta[0,5] + gamma[0,5,4] + delta[0,5,4] + epsilon[0,5] <= 0
 v[0,5,5]: gamma[0,5,5] + delta[0,5,5] + epsilon[0,5] <= 0
 v[0,5,6]: 599 beta[0,5] + gamma[0,5,6] + delta[0,5,6] + epsilon[0,5] <= 0
 v[0,5,7]: 574 beta[0,5] + gamma[0,5,7] + delta[0,5,7] + epsilon[0,5] <= 0
 v[0,5,8]: 193 beta[0,5] + gamma[0,5,8] + delta[0,5,8] + epsilon[0,5] <= 0
 v[0,5,9]: 448 beta[0,5] + gamma[0,5,9] + delta[0,5,9] + epsilon[0,5] <= 0
 v[0,6,0]: 792 beta[0,6] + gamma[0,6,0] + delta[0,6,0] + epsilon[0,6] <= 0
 v[0,6,1]: 192 beta[0,6] + gamma[0,6,1] + delta[0,6,1] + epsilon[0,6] <= 0
 v[0,6,2]: 1051 beta[0,6] + gamma[0,6,2] + delta[0,6,2] + epsilon[0,6]
   <= 0
 v[0,6,3]: 543 beta[0,6] + gamma[0,6,3] + delta[0,6,3] + epsilon[0,6] <= 0
 v[0,6,4]: 975 beta[0,6] + gamma[0,6,4] + delta[0,6,4] + epsilon[0,6] <= 0
 v[0,6,5]: 820 beta[0,6] + gamma[0,6,5] + delta[0,6,5] + epsilon[0,6] <= 0
 v[0,6,6]: gamma[0,6,6] + delta[0,6,6] + epsilon[0,6] <= 0
 v[0,6,7]: 947 beta[0,6] + gamma[0,6,7] + delta[0,6,7] + epsilon[0,6] <= 0
 v[0,6,8]: 200 beta[0,6] + gamma[0,6,8] + delta[0,6,8] + epsilon[0,6] <= 0
 v[0,6,9]: 971 beta[0,6] + gamma[0,6,9] + delta[0,6,9] + epsilon[0,6] <= 0
 v[0,7,0]: 428 beta[0,7] + gamma[0,7,0] + delta[0,7,0] + epsilon[0,7] <= 0
 v[0,7,1]: 183 beta[0,7] + gamma[0,7,1] + delta[0,7,1] + epsilon[0,7] <= 0
 v[0,7,2]: 692 beta[0,7] + gamma[0,7,2] + delta[0,7,2] + epsilon[0,7] <= 0
 v[0,7,3]: 408 beta[0,7] + gamma[0,7,3] + delta[0,7,3] + epsilon[0,7] <= 0
 v[0,7,4]: 447 beta[0,7] + gamma[0,7,4] + delta[0,7,4] + epsilon[0,7] <= 0
 v[0,7,5]: 516 beta[0,7] + gamma[0,7,5] + delta[0,7,5] + epsilon[0,7] <= 0
 v[0,7,6]: 671 beta[0,7] + gamma[0,7,6] + delta[0,7,6] + epsilon[0,7] <= 0
 v[0,7,7]: gamma[0,7,7] + delta[0,7,7] + epsilon[0,7] <= 0
 v[0,7,8]: 150 beta[0,7] + gamma[0,7,8] + delta[0,7,8] + epsilon[0,7] <= 0
 v[0,7,9]: 639 beta[0,7] + gamma[0,7,9] + delta[0,7,9] + epsilon[0,7] <= 0
 v[0,8,0]: 114 beta[0,8] + gamma[0,8,0] + delta[0,8,0] + epsilon[0,8] <= 0
 v[0,8,1]: 306 beta[0,8] + gamma[0,8,1] + delta[0,8,1] + epsilon[0,8] <= 0
 v[0,8,2]: 692 beta[0,8] + gamma[0,8,2] + delta[0,8,2] + epsilon[0,8] <= 0
 v[0,8,3]: 801 beta[0,8] + gamma[0,8,3] + delta[0,8,3] + epsilon[0,8] <= 0
 v[0,8,4]: 521 beta[0,8] + gamma[0,8,4] + delta[0,8,4] + epsilon[0,8] <= 0
 v[0,8,5]: 157 beta[0,8] + gamma[0,8,5] + delta[0,8,5] + epsilon[0,8] <= 0
 v[0,8,6]: 231 beta[0,8] + gamma[0,8,6] + delta[0,8,6] + epsilon[0,8] <= 0
 v[0,8,7]: 825 beta[0,8] + gamma[0,8,7] + delta[0,8,7] + epsilon[0,8] <= 0
 v[0,8,8]: gamma[0,8,8] + delta[0,8,8] + epsilon[0,8] <= 0
 v[0,8,9]: 163 beta[0,8] + gamma[0,8,9] + delta[0,8,9] + epsilon[0,8] <= 0
 v[0,9,0]: 279 beta[0,9] + gamma[0,9,0] + delta[0,9,0] + epsilon[0,9] <= 0
 v[0,9,1]: 1004 beta[0,9] + gamma[0,9,1] + delta[0,9,1] + epsilon[0,9]
   <= 0
 v[0,9,2]: 653 beta[0,9] + gamma[0,9,2] + delta[0,9,2] + epsilon[0,9] <= 0
 v[0,9,3]: 449 beta[0,9] + gamma[0,9,3] + delta[0,9,3] + epsilon[0,9] <= 0
 v[0,9,4]: 991 beta[0,9] + gamma[0,9,4] + delta[0,9,4] + epsilon[0,9] <= 0
 v[0,9,5]: 654 beta[0,9] + gamma[0,9,5] + delta[0,9,5] + epsilon[0,9] <= 0
 v[0,9,6]: 560 beta[0,9] + gamma[0,9,6] + delta[0,9,6] + epsilon[0,9] <= 0
 v[0,9,7]: 338 beta[0,9] + gamma[0,9,7] + delta[0,9,7] + epsilon[0,9] <= 0
 v[0,9,8]: 464 beta[0,9] + gamma[0,9,8] + delta[0,9,8] + epsilon[0,9] <= 0
 v[0,9,9]: gamma[0,9,9] + delta[0,9,9] + epsilon[0,9] <= 0
 v[1,0,0]: gamma[1,0,0] + delta[1,0,0] + epsilon[1,0] <= 0
 v[1,0,1]: 762 beta[1,0] + gamma[1,0,1] + delta[1,0,1] + epsilon[1,0] <= 0
 v[1,0,2]: 859 beta[1,0] + gamma[1,0,2] + delta[1,0,2] + epsilon[1,0] <= 0
 v[1,0,3]: 545 beta[1,0] + gamma[1,0,3] + delta[1,0,3] + epsilon[1,0] <= 0
 v[1,0,4]: 746 beta[1,0] + gamma[1,0,4] + delta[1,0,4] + epsilon[1,0] <= 0
 v[1,0,5]: 911 beta[1,0] + gamma[1,0,5] + delta[1,0,5] + epsilon[1,0] <= 0
 v[1,0,6]: 565 beta[1,0] + gamma[1,0,6] + delta[1,0,6] + epsilon[1,0] <= 0
 v[1,0,7]: 244 beta[1,0] + gamma[1,0,7] + delta[1,0,7] + epsilon[1,0] <= 0
 v[1,0,8]: 623 beta[1,0] + gamma[1,0,8] + delta[1,0,8] + epsilon[1,0] <= 0
 v[1,0,9]: 273 beta[1,0] + gamma[1,0,9] + delta[1,0,9] + epsilon[1,0] <= 0
 v[1,1,0]: 1093 beta[1,1] + gamma[1,1,0] + delta[1,1,0] + epsilon[1,1]
   <= 0
 v[1,1,1]: gamma[1,1,1] + delta[1,1,1] + epsilon[1,1] <= 0
 v[1,1,2]: 816 beta[1,1] + gamma[1,1,2] + delta[1,1,2] + epsilon[1,1] <= 0
 v[1,1,3]: 897 beta[1,1] + gamma[1,1,3] + delta[1,1,3] + epsilon[1,1] <= 0
 v[1,1,4]: 1023 beta[1,1] + gamma[1,1,4] + delta[1,1,4] + epsilon[1,1]
   <= 0
 v[1,1,5]: 943 beta[1,1] + gamma[1,1,5] + delta[1,1,5] + epsilon[1,1] <= 0
 v[1,1,6]: 374 beta[1,1] + gamma[1,1,6] + delta[1,1,6] + epsilon[1,1] <= 0
 v[1,1,7]: 844 beta[1,1] + gamma[1,1,7] + delta[1,1,7] + epsilon[1,1] <= 0
 v[1,1,8]: 933 beta[1,1] + gamma[1,1,8] + delta[1,1,8] + epsilon[1,1] <= 0
 v[1,1,9]: 163 beta[1,1] + gamma[1,1,9] + delta[1,1,9] + epsilon[1,1] <= 0
 v[1,2,0]: 1096 beta[1,2] + gamma[1,2,0] + delta[1,2,0] + epsilon[1,2]
   <= 0
 v[1,2,1]: 253 beta[1,2] + gamma[1,2,1] + delta[1,2,1] + epsilon[1,2] <= 0
 v[1,2,2]: gamma[1,2,2] + delta[1,2,2] + epsilon[1,2] <= 0
 v[1,2,3]: 183 beta[1,2] + gamma[1,2,3] + delta[1,2,3] + epsilon[1,2] <= 0
 v[1,2,4]: 122 beta[1,2] + gamma[1,2,4] + delta[1,2,4] + epsilon[1,2] <= 0
 v[1,2,5]: 218 beta[1,2] + gamma[1,2,5] + delta[1,2,5] + epsilon[1,2] <= 0
 v[1,2,6]: 940 beta[1,2] + gamma[1,2,6] + delta[1,2,6] + epsilon[1,2] <= 0
 v[1,2,7]: 1009 beta[1,2] + gamma[1,2,7] + delta[1,2,7] + epsilon[1,2]
   <= 0
 v[1,2,8]: 129 beta[1,2] + gamma[1,2,8] + delta[1,2,8] + epsilon[1,2] <= 0
 v[1,2,9]: 251 beta[1,2] + gamma[1,2,9] + delta[1,2,9] + epsilon[1,2] <= 0
 v[1,3,0]: 429 beta[1,3] + gamma[1,3,0] + delta[1,3,0] + epsilon[1,3] <= 0
 v[1,3,1]: 941 beta[1,3] + gamma[1,3,1] + delta[1,3,1] + epsilon[1,3] <= 0
 v[1,3,2]: 767 beta[1,3] + gamma[1,3,2] + delta[1,3,2] + epsilon[1,3] <= 0
 v[1,3,3]: gamma[1,3,3] + delta[1,3,3] + epsilon[1,3] <= 0
 v[1,3,4]: 775 beta[1,3] + gamma[1,3,4] + delta[1,3,4] + epsilon[1,3] <= 0
 v[1,3,5]: 768 beta[1,3] + gamma[1,3,5] + delta[1,3,5] + epsilon[1,3] <= 0
 v[1,3,6]: 1029 beta[1,3] + gamma[1,3,6] + delta[1,3,6] + epsilon[1,3]
   <= 0
 v[1,3,7]: 424 beta[1,3] + gamma[1,3,7] + delta[1,3,7] + epsilon[1,3] <= 0
 v[1,3,8]: 390 beta[1,3] + gamma[1,3,8] + delta[1,3,8] + epsilon[1,3] <= 0
 v[1,3,9]: 580 beta[1,3] + gamma[1,3,9] + delta[1,3,9] + epsilon[1,3] <= 0
 v[1,4,0]: 875 beta[1,4] + gamma[1,4,0] + delta[1,4,0] + epsilon[1,4] <= 0
 v[1,4,1]: 435 beta[1,4] + gamma[1,4,1] + delta[1,4,1] + epsilon[1,4] <= 0
 v[1,4,2]: 211 beta[1,4] + gamma[1,4,2] + delta[1,4,2] + epsilon[1,4] <= 0
 v[1,4,3]: 1010 beta[1,4] + gamma[1,4,3] + delta[1,4,3] + epsilon[1,4]
   <= 0
 v[1,4,4]: gamma[1,4,4] + delta[1,4,4] + epsilon[1,4] <= 0
 v[1,4,5]: 1027 beta[1,4] + gamma[1,4,5] + delta[1,4,5] + epsilon[1,4]
   <= 0
 v[1,4,6]: 341 beta[1,4] + gamma[1,4,6] + delta[1,4,6] + epsilon[1,4] <= 0
 v[1,4,7]: 308 beta[1,4] + gamma[1,4,7] + delta[1,4,7] + epsilon[1,4] <= 0
 v[1,4,8]: 294 beta[1,4] + gamma[1,4,8] + delta[1,4,8] + epsilon[1,4] <= 0
 v[1,4,9]: 415 beta[1,4] + gamma[1,4,9] + delta[1,4,9] + epsilon[1,4] <= 0
 v[1,5,0]: 526 beta[1,5] + gamma[1,5,0] + delta[1,5,0] + epsilon[1,5] <= 0
 v[1,5,1]: 1048 beta[1,5] + gamma[1,5,1] + delta[1,5,1] + epsilon[1,5]
   <= 0
 v[1,5,2]: 127 beta[1,5] + gamma[1,5,2] + delta[1,5,2] + epsilon[1,5] <= 0
 v[1,5,3]: 335 beta[1,5] + gamma[1,5,3] + delta[1,5,3] + epsilon[1,5] <= 0
 v[1,5,4]: 371 beta[1,5] + gamma[1,5,4] + delta[1,5,4] + epsilon[1,5] <= 0
 v[1,5,5]: gamma[1,5,5] + delta[1,5,5] + epsilon[1,5] <= 0
 v[1,5,6]: 171 beta[1,5] + gamma[1,5,6] + delta[1,5,6] + epsilon[1,5] <= 0
 v[1,5,7]: 1001 beta[1,5] + gamma[1,5,7] + delta[1,5,7] + epsilon[1,5]
   <= 0
 v[1,5,8]: 989 beta[1,5] + gamma[1,5,8] + delta[1,5,8] + epsilon[1,5] <= 0
 v[1,5,9]: 653 beta[1,5] + gamma[1,5,9] + delta[1,5,9] + epsilon[1,5] <= 0
 v[1,6,0]: 209 beta[1,6] + gamma[1,6,0] + delta[1,6,0] + epsilon[1,6] <= 0
 v[1,6,1]: 1054 beta[1,6] + gamma[1,6,1] + delta[1,6,1] + epsilon[1,6]
   <= 0
 v[1,6,2]: 538 beta[1,6] + gamma[1,6,2] + delta[1,6,2] + epsilon[1,6] <= 0
 v[1,6,3]: 116 beta[1,6] + gamma[1,6,3] + delta[1,6,3] + epsilon[1,6] <= 0
 v[1,6,4]: 943 beta[1,6] + gamma[1,6,4] + delta[1,6,4] + epsilon[1,6] <= 0
 v[1,6,5]: 185 beta[1,6] + gamma[1,6,5] + delta[1,6,5] + epsilon[1,6] <= 0
 v[1,6,6]: gamma[1,6,6] + delta[1,6,6] + epsilon[1,6] <= 0
 v[1,6,7]: 172 beta[1,6] + gamma[1,6,7] + delta[1,6,7] + epsilon[1,6] <= 0
 v[1,6,8]: 779 beta[1,6] + gamma[1,6,8] + delta[1,6,8] + epsilon[1,6] <= 0
 v[1,6,9]: 221 beta[1,6] + gamma[1,6,9] + delta[1,6,9] + epsilon[1,6] <= 0
 v[1,7,0]: 871 beta[1,7] + gamma[1,7,0] + delta[1,7,0] + epsilon[1,7] <= 0
 v[1,7,1]: 339 beta[1,7] + gamma[1,7,1] + delta[1,7,1] + epsilon[1,7] <= 0
 v[1,7,2]: 208 beta[1,7] + gamma[1,7,2] + delta[1,7,2] + epsilon[1,7] <= 0
 v[1,7,3]: 182 beta[1,7] + gamma[1,7,3] + delta[1,7,3] + epsilon[1,7] <= 0
 v[1,7,4]: 633 beta[1,7] + gamma[1,7,4] + delta[1,7,4] + epsilon[1,7] <= 0
 v[1,7,5]: 727 beta[1,7] + gamma[1,7,5] + delta[1,7,5] + epsilon[1,7] <= 0
 v[1,7,6]: 979 beta[1,7] + gamma[1,7,6] + delta[1,7,6] + epsilon[1,7] <= 0
 v[1,7,7]: gamma[1,7,7] + delta[1,7,7] + epsilon[1,7] <= 0
 v[1,7,8]: 438 beta[1,7] + gamma[1,7,8] + delta[1,7,8] + epsilon[1,7] <= 0
 v[1,7,9]: 762 beta[1,7] + gamma[1,7,9] + delta[1,7,9] + epsilon[1,7] <= 0
 v[1,8,0]: 533 beta[1,8] + gamma[1,8,0] + delta[1,8,0] + epsilon[1,8] <= 0
 v[1,8,1]: 310 beta[1,8] + gamma[1,8,1] + delta[1,8,1] + epsilon[1,8] <= 0
 v[1,8,2]: 234 beta[1,8] + gamma[1,8,2] + delta[1,8,2] + epsilon[1,8] <= 0
 v[1,8,3]: 162 beta[1,8] + gamma[1,8,3] + delta[1,8,3] + epsilon[1,8] <= 0
 v[1,8,4]: 440 beta[1,8] + gamma[1,8,4] + delta[1,8,4] + epsilon[1,8] <= 0
 v[1,8,5]: 471 beta[1,8] + gamma[1,8,5] + delta[1,8,5] + epsilon[1,8] <= 0
 v[1,8,6]: 894 beta[1,8] + gamma[1,8,6] + delta[1,8,6] + epsilon[1,8] <= 0
 v[1,8,7]: 211 beta[1,8] + gamma[1,8,7] + delta[1,8,7] + epsilon[1,8] <= 0
 v[1,8,8]: gamma[1,8,8] + delta[1,8,8] + epsilon[1,8] <= 0
 v[1,8,9]: 907 beta[1,8] + gamma[1,8,9] + delta[1,8,9] + epsilon[1,8] <= 0
 v[1,9,0]: 420 beta[1,9] + gamma[1,9,0] + delta[1,9,0] + epsilon[1,9] <= 0
 v[1,9,1]: 405 beta[1,9] + gamma[1,9,1] + delta[1,9,1] + epsilon[1,9] <= 0
 v[1,9,2]: 652 beta[1,9] + gamma[1,9,2] + delta[1,9,2] + epsilon[1,9] <= 0
 v[1,9,3]: 709 beta[1,9] + gamma[1,9,3] + delta[1,9,3] + epsilon[1,9] <= 0
 v[1,9,4]: 132 beta[1,9] + gamma[1,9,4] + delta[1,9,4] + epsilon[1,9] <= 0
 v[1,9,5]: 247 beta[1,9] + gamma[1,9,5] + delta[1,9,5] + epsilon[1,9] <= 0
 v[1,9,6]: 752 beta[1,9] + gamma[1,9,6] + delta[1,9,6] + epsilon[1,9] <= 0
 v[1,9,7]: 916 beta[1,9] + gamma[1,9,7] + delta[1,9,7] + epsilon[1,9] <= 0
 v[1,9,8]: 244 beta[1,9] + gamma[1,9,8] + delta[1,9,8] + epsilon[1,9] <= 0
 v[1,9,9]: gamma[1,9,9] + delta[1,9,9] + epsilon[1,9] <= 0
 v[2,0,0]: gamma[2,0,0] + delta[2,0,0] + epsilon[2,0] <= 0
 v[2,0,1]: 818 beta[2,0] + gamma[2,0,1] + delta[2,0,1] + epsilon[2,0] <= 0
 v[2,0,2]: 102 beta[2,0] + gamma[2,0,2] + delta[2,0,2] + epsilon[2,0] <= 0
 v[2,0,3]: 810 beta[2,0] + gamma[2,0,3] + delta[2,0,3] + epsilon[2,0] <= 0
 v[2,0,4]: 385 beta[2,0] + gamma[2,0,4] + delta[2,0,4] + epsilon[2,0] <= 0
 v[2,0,5]: 915 beta[2,0] + gamma[2,0,5] + delta[2,0,5] + epsilon[2,0] <= 0
 v[2,0,6]: 862 beta[2,0] + gamma[2,0,6] + delta[2,0,6] + epsilon[2,0] <= 0
 v[2,0,7]: 485 beta[2,0] + gamma[2,0,7] + delta[2,0,7] + epsilon[2,0] <= 0
 v[2,0,8]: 842 beta[2,0] + gamma[2,0,8] + delta[2,0,8] + epsilon[2,0] <= 0
 v[2,0,9]: 590 beta[2,0] + gamma[2,0,9] + delta[2,0,9] + epsilon[2,0] <= 0
 v[2,1,0]: 115 beta[2,1] + gamma[2,1,0] + delta[2,1,0] + epsilon[2,1] <= 0
 v[2,1,1]: gamma[2,1,1] + delta[2,1,1] + epsilon[2,1] <= 0
 v[2,1,2]: 751 beta[2,1] + gamma[2,1,2] + delta[2,1,2] + epsilon[2,1] <= 0
 v[2,1,3]: 379 beta[2,1] + gamma[2,1,3] + delta[2,1,3] + epsilon[2,1] <= 0
 v[2,1,4]: 682 beta[2,1] + gamma[2,1,4] + delta[2,1,4] + epsilon[2,1] <= 0
 v[2,1,5]: 1099 beta[2,1] + gamma[2,1,5] + delta[2,1,5] + epsilon[2,1]
   <= 0
 v[2,1,6]: 833 beta[2,1] + gamma[2,1,6] + delta[2,1,6] + epsilon[2,1] <= 0
 v[2,1,7]: 583 beta[2,1] + gamma[2,1,7] + delta[2,1,7] + epsilon[2,1] <= 0
 v[2,1,8]: 730 beta[2,1] + gamma[2,1,8] + delta[2,1,8] + epsilon[2,1] <= 0
 v[2,1,9]: 889 beta[2,1] + gamma[2,1,9] + delta[2,1,9] + epsilon[2,1] <= 0
 v[2,2,0]: 805 beta[2,2] + gamma[2,2,0] + delta[2,2,0] + epsilon[2,2] <= 0
 v[2,2,1]: 498 beta[2,2] + gamma[2,2,1] + delta[2,2,1] + epsilon[2,2] <= 0
 v[2,2,2]: gamma[2,2,2] + delta[2,2,2] + epsilon[2,2] <= 0
 v[2,2,3]: 385 beta[2,2] + gamma[2,2,3] + delta[2,2,3] + epsilon[2,2] <= 0
 v[2,2,4]: 191 beta[2,2] + gamma[2,2,4] + delta[2,2,4] + epsilon[2,2] <= 0
 v[2,2,5]: 411 beta[2,2] + gamma[2,2,5] + delta[2,2,5] + epsilon[2,2] <= 0
 v[2,2,6]: 498 beta[2,2] + gamma[2,2,6] + delta[2,2,6] + epsilon[2,2] <= 0
 v[2,2,7]: 481 beta[2,2] + gamma[2,2,7] + delta[2,2,7] + epsilon[2,2] <= 0
 v[2,2,8]: 760 beta[2,2] + gamma[2,2,8] + delta[2,2,8] + epsilon[2,2] <= 0
 v[2,2,9]: 436 beta[2,2] + gamma[2,2,9] + delta[2,2,9] + epsilon[2,2] <= 0
 v[2,3,0]: 511 beta[2,3] + gamma[2,3,0] + delta[2,3,0] + epsilon[2,3] <= 0
 v[2,3,1]: 513 beta[2,3] + gamma[2,3,1] + delta[2,3,1] + epsilon[2,3] <= 0
 v[2,3,2]: 864 beta[2,3] + gamma[2,3,2] + delta[2,3,2] + epsilon[2,3] <= 0
 v[2,3,3]: gamma[2,3,3] + delta[2,3,3] + epsilon[2,3] <= 0
 v[2,3,4]: 831 beta[2,3] + gamma[2,3,4] + delta[2,3,4] + epsilon[2,3] <= 0
 v[2,3,5]: 577 beta[2,3] + gamma[2,3,5] + delta[2,3,5] + epsilon[2,3] <= 0
 v[2,3,6]: 490 beta[2,3] + gamma[2,3,6] + delta[2,3,6] + epsilon[2,3] <= 0
 v[2,3,7]: 213 beta[2,3] + gamma[2,3,7] + delta[2,3,7] + epsilon[2,3] <= 0
 v[2,3,8]: 934 beta[2,3] + gamma[2,3,8] + delta[2,3,8] + epsilon[2,3] <= 0
 v[2,3,9]: 739 beta[2,3] + gamma[2,3,9] + delta[2,3,9] + epsilon[2,3] <= 0
 v[2,4,0]: 625 beta[2,4] + gamma[2,4,0] + delta[2,4,0] + epsilon[2,4] <= 0
 v[2,4,1]: 507 beta[2,4] + gamma[2,4,1] + delta[2,4,1] + epsilon[2,4] <= 0
 v[2,4,2]: 813 beta[2,4] + gamma[2,4,2] + delta[2,4,2] + epsilon[2,4] <= 0
 v[2,4,3]: 202 beta[2,4] + gamma[2,4,3] + delta[2,4,3] + epsilon[2,4] <= 0
 v[2,4,4]: gamma[2,4,4] + delta[2,4,4] + epsilon[2,4] <= 0
 v[2,4,5]: 499 beta[2,4] + gamma[2,4,5] + delta[2,4,5] + epsilon[2,4] <= 0
 v[2,4,6]: 447 beta[2,4] + gamma[2,4,6] + delta[2,4,6] + epsilon[2,4] <= 0
 v[2,4,7]: 243 beta[2,4] + gamma[2,4,7] + delta[2,4,7] + epsilon[2,4] <= 0
 v[2,4,8]: 294 beta[2,4] + gamma[2,4,8] + delta[2,4,8] + epsilon[2,4] <= 0
 v[2,4,9]: 785 beta[2,4] + gamma[2,4,9] + delta[2,4,9] + epsilon[2,4] <= 0
 v[2,5,0]: 229 beta[2,5] + gamma[2,5,0] + delta[2,5,0] + epsilon[2,5] <= 0
 v[2,5,1]: 846 beta[2,5] + gamma[2,5,1] + delta[2,5,1] + epsilon[2,5] <= 0
 v[2,5,2]: 602 beta[2,5] + gamma[2,5,2] + delta[2,5,2] + epsilon[2,5] <= 0
 v[2,5,3]: 541 beta[2,5] + gamma[2,5,3] + delta[2,5,3] + epsilon[2,5] <= 0
 v[2,5,4]: 482 beta[2,5] + gamma[2,5,4] + delta[2,5,4] + epsilon[2,5] <= 0
 v[2,5,5]: gamma[2,5,5] + delta[2,5,5] + epsilon[2,5] <= 0
 v[2,5,6]: 863 beta[2,5] + gamma[2,5,6] + delta[2,5,6] + epsilon[2,5] <= 0
 v[2,5,7]: 523 beta[2,5] + gamma[2,5,7] + delta[2,5,7] + epsilon[2,5] <= 0
 v[2,5,8]: 465 beta[2,5] + gamma[2,5,8] + delta[2,5,8] + epsilon[2,5] <= 0
 v[2,5,9]: 922 beta[2,5] + gamma[2,5,9] + delta[2,5,9] + epsilon[2,5] <= 0
 v[2,6,0]: 1092 beta[2,6] + gamma[2,6,0] + delta[2,6,0] + epsilon[2,6]
   <= 0
 v[2,6,1]: 554 beta[2,6] + gamma[2,6,1] + delta[2,6,1] + epsilon[2,6] <= 0
 v[2,6,2]: 1051 beta[2,6] + gamma[2,6,2] + delta[2,6,2] + epsilon[2,6]
   <= 0
 v[2,6,3]: 112 beta[2,6] + gamma[2,6,3] + delta[2,6,3] + epsilon[2,6] <= 0
 v[2,6,4]: 880 beta[2,6] + gamma[2,6,4] + delta[2,6,4] + epsilon[2,6] <= 0
 v[2,6,5]: 1020 beta[2,6] + gamma[2,6,5] + delta[2,6,5] + epsilon[2,6]
   <= 0
 v[2,6,6]: gamma[2,6,6] + delta[2,6,6] + epsilon[2,6] <= 0
 v[2,6,7]: 1043 beta[2,6] + gamma[2,6,7] + delta[2,6,7] + epsilon[2,6]
   <= 0
 v[2,6,8]: 366 beta[2,6] + gamma[2,6,8] + delta[2,6,8] + epsilon[2,6] <= 0
 v[2,6,9]: 1015 beta[2,6] + gamma[2,6,9] + delta[2,6,9] + epsilon[2,6]
   <= 0
 v[2,7,0]: 743 beta[2,7] + gamma[2,7,0] + delta[2,7,0] + epsilon[2,7] <= 0
 v[2,7,1]: 886 beta[2,7] + gamma[2,7,1] + delta[2,7,1] + epsilon[2,7] <= 0
 v[2,7,2]: 253 beta[2,7] + gamma[2,7,2] + delta[2,7,2] + epsilon[2,7] <= 0
 v[2,7,3]: 723 beta[2,7] + gamma[2,7,3] + delta[2,7,3] + epsilon[2,7] <= 0
 v[2,7,4]: 441 beta[2,7] + gamma[2,7,4] + delta[2,7,4] + epsilon[2,7] <= 0
 v[2,7,5]: 746 beta[2,7] + gamma[2,7,5] + delta[2,7,5] + epsilon[2,7] <= 0
 v[2,7,6]: 257 beta[2,7] + gamma[2,7,6] + delta[2,7,6] + epsilon[2,7] <= 0
 v[2,7,7]: gamma[2,7,7] + delta[2,7,7] + epsilon[2,7] <= 0
 v[2,7,8]: 1005 beta[2,7] + gamma[2,7,8] + delta[2,7,8] + epsilon[2,7]
   <= 0
 v[2,7,9]: 102 beta[2,7] + gamma[2,7,9] + delta[2,7,9] + epsilon[2,7] <= 0
 v[2,8,0]: 310 beta[2,8] + gamma[2,8,0] + delta[2,8,0] + epsilon[2,8] <= 0
 v[2,8,1]: 1023 beta[2,8] + gamma[2,8,1] + delta[2,8,1] + epsilon[2,8]
   <= 0
 v[2,8,2]: 991 beta[2,8] + gamma[2,8,2] + delta[2,8,2] + epsilon[2,8] <= 0
 v[2,8,3]: 511 beta[2,8] + gamma[2,8,3] + delta[2,8,3] + epsilon[2,8] <= 0
 v[2,8,4]: 970 beta[2,8] + gamma[2,8,4] + delta[2,8,4] + epsilon[2,8] <= 0
 v[2,8,5]: 261 beta[2,8] + gamma[2,8,5] + delta[2,8,5] + epsilon[2,8] <= 0
 v[2,8,6]: 140 beta[2,8] + gamma[2,8,6] + delta[2,8,6] + epsilon[2,8] <= 0
 v[2,8,7]: 103 beta[2,8] + gamma[2,8,7] + delta[2,8,7] + epsilon[2,8] <= 0
 v[2,8,8]: gamma[2,8,8] + delta[2,8,8] + epsilon[2,8] <= 0
 v[2,8,9]: 214 beta[2,8] + gamma[2,8,9] + delta[2,8,9] + epsilon[2,8] <= 0
 v[2,9,0]: 901 beta[2,9] + gamma[2,9,0] + delta[2,9,0] + epsilon[2,9] <= 0
 v[2,9,1]: 687 beta[2,9] + gamma[2,9,1] + delta[2,9,1] + epsilon[2,9] <= 0
 v[2,9,2]: 471 beta[2,9] + gamma[2,9,2] + delta[2,9,2] + epsilon[2,9] <= 0
 v[2,9,3]: 817 beta[2,9] + gamma[2,9,3] + delta[2,9,3] + epsilon[2,9] <= 0
 v[2,9,4]: 399 beta[2,9] + gamma[2,9,4] + delta[2,9,4] + epsilon[2,9] <= 0
 v[2,9,5]: 917 beta[2,9] + gamma[2,9,5] + delta[2,9,5] + epsilon[2,9] <= 0
 v[2,9,6]: 357 beta[2,9] + gamma[2,9,6] + delta[2,9,6] + epsilon[2,9] <= 0
 v[2,9,7]: 475 beta[2,9] + gamma[2,9,7] + delta[2,9,7] + epsilon[2,9] <= 0
 v[2,9,8]: 867 beta[2,9] + gamma[2,9,8] + delta[2,9,8] + epsilon[2,9] <= 0
 v[2,9,9]: gamma[2,9,9] + delta[2,9,9] + epsilon[2,9] <= 0
 v[3,0,0]: gamma[3,0,0] + delta[3,0,0] + epsilon[3,0] <= 0
 v[3,0,1]: 715 beta[3,0] + gamma[3,0,1] + delta[3,0,1] + epsilon[3,0] <= 0
 v[3,0,2]: 232 beta[3,0] + gamma[3,0,2] + delta[3,0,2] + epsilon[3,0] <= 0
 v[3,0,3]: 337 beta[3,0] + gamma[3,0,3] + delta[3,0,3] + epsilon[3,0] <= 0
 v[3,0,4]: 133 beta[3,0] + gamma[3,0,4] + delta[3,0,4] + epsilon[3,0] <= 0
 v[3,0,5]: 1034 beta[3,0] + gamma[3,0,5] + delta[3,0,5] + epsilon[3,0]
   <= 0
 v[3,0,6]: 183 beta[3,0] + gamma[3,0,6] + delta[3,0,6] + epsilon[3,0] <= 0
 v[3,0,7]: 683 beta[3,0] + gamma[3,0,7] + delta[3,0,7] + epsilon[3,0] <= 0
 v[3,0,8]: 1007 beta[3,0] + gamma[3,0,8] + delta[3,0,8] + epsilon[3,0]
   <= 0
 v[3,0,9]: 645 beta[3,0] + gamma[3,0,9] + delta[3,0,9] + epsilon[3,0] <= 0
 v[3,1,0]: 457 beta[3,1] + gamma[3,1,0] + delta[3,1,0] + epsilon[3,1] <= 0
 v[3,1,1]: gamma[3,1,1] + delta[3,1,1] + epsilon[3,1] <= 0
 v[3,1,2]: 620 beta[3,1] + gamma[3,1,2] + delta[3,1,2] + epsilon[3,1] <= 0
 v[3,1,3]: 773 beta[3,1] + gamma[3,1,3] + delta[3,1,3] + epsilon[3,1] <= 0
 v[3,1,4]: 1037 beta[3,1] + gamma[3,1,4] + delta[3,1,4] + epsilon[3,1]
   <= 0
 v[3,1,5]: 485 beta[3,1] + gamma[3,1,5] + delta[3,1,5] + epsilon[3,1] <= 0
 v[3,1,6]: 659 beta[3,1] + gamma[3,1,6] + delta[3,1,6] + epsilon[3,1] <= 0
 v[3,1,7]: 371 beta[3,1] + gamma[3,1,7] + delta[3,1,7] + epsilon[3,1] <= 0
 v[3,1,8]: 263 beta[3,1] + gamma[3,1,8] + delta[3,1,8] + epsilon[3,1] <= 0
 v[3,1,9]: 819 beta[3,1] + gamma[3,1,9] + delta[3,1,9] + epsilon[3,1] <= 0
 v[3,2,0]: 562 beta[3,2] + gamma[3,2,0] + delta[3,2,0] + epsilon[3,2] <= 0
 v[3,2,1]: 251 beta[3,2] + gamma[3,2,1] + delta[3,2,1] + epsilon[3,2] <= 0
 v[3,2,2]: gamma[3,2,2] + delta[3,2,2] + epsilon[3,2] <= 0
 v[3,2,3]: 621 beta[3,2] + gamma[3,2,3] + delta[3,2,3] + epsilon[3,2] <= 0
 v[3,2,4]: 539 beta[3,2] + gamma[3,2,4] + delta[3,2,4] + epsilon[3,2] <= 0
 v[3,2,5]: 624 beta[3,2] + gamma[3,2,5] + delta[3,2,5] + epsilon[3,2] <= 0
 v[3,2,6]: 624 beta[3,2] + gamma[3,2,6] + delta[3,2,6] + epsilon[3,2] <= 0
 v[3,2,7]: 1019 beta[3,2] + gamma[3,2,7] + delta[3,2,7] + epsilon[3,2]
   <= 0
 v[3,2,8]: 723 beta[3,2] + gamma[3,2,8] + delta[3,2,8] + epsilon[3,2] <= 0
 v[3,2,9]: 287 beta[3,2] + gamma[3,2,9] + delta[3,2,9] + epsilon[3,2] <= 0
 v[3,3,0]: 1019 beta[3,3] + gamma[3,3,0] + delta[3,3,0] + epsilon[3,3]
   <= 0
 v[3,3,1]: 164 beta[3,3] + gamma[3,3,1] + delta[3,3,1] + epsilon[3,3] <= 0
 v[3,3,2]: 619 beta[3,3] + gamma[3,3,2] + delta[3,3,2] + epsilon[3,3] <= 0
 v[3,3,3]: gamma[3,3,3] + delta[3,3,3] + epsilon[3,3] <= 0
 v[3,3,4]: 260 beta[3,3] + gamma[3,3,4] + delta[3,3,4] + epsilon[3,3] <= 0
 v[3,3,5]: 327 beta[3,3] + gamma[3,3,5] + delta[3,3,5] + epsilon[3,3] <= 0
 v[3,3,6]: 1100 beta[3,3] + gamma[3,3,6] + delta[3,3,6] + epsilon[3,3]
   <= 0
 v[3,3,7]: 1014 beta[3,3] + gamma[3,3,7] + delta[3,3,7] + epsilon[3,3]
   <= 0
 v[3,3,8]: 330 beta[3,3] + gamma[3,3,8] + delta[3,3,8] + epsilon[3,3] <= 0
 v[3,3,9]: 664 beta[3,3] + gamma[3,3,9] + delta[3,3,9] + epsilon[3,3] <= 0
 v[3,4,0]: 503 beta[3,4] + gamma[3,4,0] + delta[3,4,0] + epsilon[3,4] <= 0
 v[3,4,1]: 278 beta[3,4] + gamma[3,4,1] + delta[3,4,1] + epsilon[3,4] <= 0
 v[3,4,2]: 474 beta[3,4] + gamma[3,4,2] + delta[3,4,2] + epsilon[3,4] <= 0
 v[3,4,3]: 466 beta[3,4] + gamma[3,4,3] + delta[3,4,3] + epsilon[3,4] <= 0
 v[3,4,4]: gamma[3,4,4] + delta[3,4,4] + epsilon[3,4] <= 0
 v[3,4,5]: 787 beta[3,4] + gamma[3,4,5] + delta[3,4,5] + epsilon[3,4] <= 0
 v[3,4,6]: 142 beta[3,4] + gamma[3,4,6] + delta[3,4,6] + epsilon[3,4] <= 0
 v[3,4,7]: 485 beta[3,4] + gamma[3,4,7] + delta[3,4,7] + epsilon[3,4] <= 0
 v[3,4,8]: 849 beta[3,4] + gamma[3,4,8] + delta[3,4,8] + epsilon[3,4] <= 0
 v[3,4,9]: 342 beta[3,4] + gamma[3,4,9] + delta[3,4,9] + epsilon[3,4] <= 0
 v[3,5,0]: 323 beta[3,5] + gamma[3,5,0] + delta[3,5,0] + epsilon[3,5] <= 0
 v[3,5,1]: 908 beta[3,5] + gamma[3,5,1] + delta[3,5,1] + epsilon[3,5] <= 0
 v[3,5,2]: 957 beta[3,5] + gamma[3,5,2] + delta[3,5,2] + epsilon[3,5] <= 0
 v[3,5,3]: 706 beta[3,5] + gamma[3,5,3] + delta[3,5,3] + epsilon[3,5] <= 0
 v[3,5,4]: 843 beta[3,5] + gamma[3,5,4] + delta[3,5,4] + epsilon[3,5] <= 0
 v[3,5,5]: gamma[3,5,5] + delta[3,5,5] + epsilon[3,5] <= 0
 v[3,5,6]: 661 beta[3,5] + gamma[3,5,6] + delta[3,5,6] + epsilon[3,5] <= 0
 v[3,5,7]: 753 beta[3,5] + gamma[3,5,7] + delta[3,5,7] + epsilon[3,5] <= 0
 v[3,5,8]: 685 beta[3,5] + gamma[3,5,8] + delta[3,5,8] + epsilon[3,5] <= 0
 v[3,5,9]: 280 beta[3,5] + gamma[3,5,9] + delta[3,5,9] + epsilon[3,5] <= 0
 v[3,6,0]: 371 beta[3,6] + gamma[3,6,0] + delta[3,6,0] + epsilon[3,6] <= 0
 v[3,6,1]: 401 beta[3,6] + gamma[3,6,1] + delta[3,6,1] + epsilon[3,6] <= 0
 v[3,6,2]: 516 beta[3,6] + gamma[3,6,2] + delta[3,6,2] + epsilon[3,6] <= 0
 v[3,6,3]: 671 beta[3,6] + gamma[3,6,3] + delta[3,6,3] + epsilon[3,6] <= 0
 v[3,6,4]: 199 beta[3,6] + gamma[3,6,4] + delta[3,6,4] + epsilon[3,6] <= 0
 v[3,6,5]: 765 beta[3,6] + gamma[3,6,5] + delta[3,6,5] + epsilon[3,6] <= 0
 v[3,6,6]: gamma[3,6,6] + delta[3,6,6] + epsilon[3,6] <= 0
 v[3,6,7]: 520 beta[3,6] + gamma[3,6,7] + delta[3,6,7] + epsilon[3,6] <= 0
 v[3,6,8]: 646 beta[3,6] + gamma[3,6,8] + delta[3,6,8] + epsilon[3,6] <= 0
 v[3,6,9]: 720 beta[3,6] + gamma[3,6,9] + delta[3,6,9] + epsilon[3,6] <= 0
 v[3,7,0]: 814 beta[3,7] + gamma[3,7,0] + delta[3,7,0] + epsilon[3,7] <= 0
 v[3,7,1]: 164 beta[3,7] + gamma[3,7,1] + delta[3,7,1] + epsilon[3,7] <= 0
 v[3,7,2]: 357 beta[3,7] + gamma[3,7,2] + delta[3,7,2] + epsilon[3,7] <= 0
 v[3,7,3]: 372 beta[3,7] + gamma[3,7,3] + delta[3,7,3] + epsilon[3,7] <= 0
 v[3,7,4]: 1042 beta[3,7] + gamma[3,7,4] + delta[3,7,4] + epsilon[3,7]
   <= 0
 v[3,7,5]: 791 beta[3,7] + gamma[3,7,5] + delta[3,7,5] + epsilon[3,7] <= 0
 v[3,7,6]: 630 beta[3,7] + gamma[3,7,6] + delta[3,7,6] + epsilon[3,7] <= 0
 v[3,7,7]: gamma[3,7,7] + delta[3,7,7] + epsilon[3,7] <= 0
 v[3,7,8]: 940 beta[3,7] + gamma[3,7,8] + delta[3,7,8] + epsilon[3,7] <= 0
 v[3,7,9]: 714 beta[3,7] + gamma[3,7,9] + delta[3,7,9] + epsilon[3,7] <= 0
 v[3,8,0]: 1083 beta[3,8] + gamma[3,8,0] + delta[3,8,0] + epsilon[3,8]
   <= 0
 v[3,8,1]: 919 beta[3,8] + gamma[3,8,1] + delta[3,8,1] + epsilon[3,8] <= 0
 v[3,8,2]: 1069 beta[3,8] + gamma[3,8,2] + delta[3,8,2] + epsilon[3,8]
   <= 0
 v[3,8,3]: 212 beta[3,8] + gamma[3,8,3] + delta[3,8,3] + epsilon[3,8] <= 0
 v[3,8,4]: 999 beta[3,8] + gamma[3,8,4] + delta[3,8,4] + epsilon[3,8] <= 0
 v[3,8,5]: 162 beta[3,8] + gamma[3,8,5] + delta[3,8,5] + epsilon[3,8] <= 0
 v[3,8,6]: 416 beta[3,8] + gamma[3,8,6] + delta[3,8,6] + epsilon[3,8] <= 0
 v[3,8,7]: 768 beta[3,8] + gamma[3,8,7] + delta[3,8,7] + epsilon[3,8] <= 0
 v[3,8,8]: gamma[3,8,8] + delta[3,8,8] + epsilon[3,8] <= 0
 v[3,8,9]: 965 beta[3,8] + gamma[3,8,9] + delta[3,8,9] + epsilon[3,8] <= 0
 v[3,9,0]: 795 beta[3,9] + gamma[3,9,0] + delta[3,9,0] + epsilon[3,9] <= 0
 v[3,9,1]: 893 beta[3,9] + gamma[3,9,1] + delta[3,9,1] + epsilon[3,9] <= 0
 v[3,9,2]: 658 beta[3,9] + gamma[3,9,2] + delta[3,9,2] + epsilon[3,9] <= 0
 v[3,9,3]: 270 beta[3,9] + gamma[3,9,3] + delta[3,9,3] + epsilon[3,9] <= 0
 v[3,9,4]: 668 beta[3,9] + gamma[3,9,4] + delta[3,9,4] + epsilon[3,9] <= 0
 v[3,9,5]: 1053 beta[3,9] + gamma[3,9,5] + delta[3,9,5] + epsilon[3,9]
   <= 0
 v[3,9,6]: 650 beta[3,9] + gamma[3,9,6] + delta[3,9,6] + epsilon[3,9] <= 0
 v[3,9,7]: 205 beta[3,9] + gamma[3,9,7] + delta[3,9,7] + epsilon[3,9] <= 0
 v[3,9,8]: 620 beta[3,9] + gamma[3,9,8] + delta[3,9,8] + epsilon[3,9] <= 0
 v[3,9,9]: gamma[3,9,9] + delta[3,9,9] + epsilon[3,9] <= 0
 v[4,0,0]: gamma[4,0,0] + delta[4,0,0] + epsilon[4,0] <= 0
 v[4,0,1]: 1080 beta[4,0] + gamma[4,0,1] + delta[4,0,1] + epsilon[4,0]
   <= 0
 v[4,0,2]: 682 beta[4,0] + gamma[4,0,2] + delta[4,0,2] + epsilon[4,0] <= 0
 v[4,0,3]: 1075 beta[4,0] + gamma[4,0,3] + delta[4,0,3] + epsilon[4,0]
   <= 0
 v[4,0,4]: 1019 beta[4,0] + gamma[4,0,4] + delta[4,0,4] + epsilon[4,0]
   <= 0
 v[4,0,5]: 424 beta[4,0] + gamma[4,0,5] + delta[4,0,5] + epsilon[4,0] <= 0
 v[4,0,6]: 790 beta[4,0] + gamma[4,0,6] + delta[4,0,6] + epsilon[4,0] <= 0
 v[4,0,7]: 138 beta[4,0] + gamma[4,0,7] + delta[4,0,7] + epsilon[4,0] <= 0
 v[4,0,8]: 761 beta[4,0] + gamma[4,0,8] + delta[4,0,8] + epsilon[4,0] <= 0
 v[4,0,9]: 574 beta[4,0] + gamma[4,0,9] + delta[4,0,9] + epsilon[4,0] <= 0
 v[4,1,0]: 552 beta[4,1] + gamma[4,1,0] + delta[4,1,0] + epsilon[4,1] <= 0
 v[4,1,1]: gamma[4,1,1] + delta[4,1,1] + epsilon[4,1] <= 0
 v[4,1,2]: 413 beta[4,1] + gamma[4,1,2] + delta[4,1,2] + epsilon[4,1] <= 0
 v[4,1,3]: 378 beta[4,1] + gamma[4,1,3] + delta[4,1,3] + epsilon[4,1] <= 0
 v[4,1,4]: 185 beta[4,1] + gamma[4,1,4] + delta[4,1,4] + epsilon[4,1] <= 0
 v[4,1,5]: 418 beta[4,1] + gamma[4,1,5] + delta[4,1,5] + epsilon[4,1] <= 0
 v[4,1,6]: 425 beta[4,1] + gamma[4,1,6] + delta[4,1,6] + epsilon[4,1] <= 0
 v[4,1,7]: 150 beta[4,1] + gamma[4,1,7] + delta[4,1,7] + epsilon[4,1] <= 0
 v[4,1,8]: 735 beta[4,1] + gamma[4,1,8] + delta[4,1,8] + epsilon[4,1] <= 0
 v[4,1,9]: 846 beta[4,1] + gamma[4,1,9] + delta[4,1,9] + epsilon[4,1] <= 0
 v[4,2,0]: 410 beta[4,2] + gamma[4,2,0] + delta[4,2,0] + epsilon[4,2] <= 0
 v[4,2,1]: 732 beta[4,2] + gamma[4,2,1] + delta[4,2,1] + epsilon[4,2] <= 0
 v[4,2,2]: gamma[4,2,2] + delta[4,2,2] + epsilon[4,2] <= 0
 v[4,2,3]: 991 beta[4,2] + gamma[4,2,3] + delta[4,2,3] + epsilon[4,2] <= 0
 v[4,2,4]: 810 beta[4,2] + gamma[4,2,4] + delta[4,2,4] + epsilon[4,2] <= 0
 v[4,2,5]: 862 beta[4,2] + gamma[4,2,5] + delta[4,2,5] + epsilon[4,2] <= 0
 v[4,2,6]: 985 beta[4,2] + gamma[4,2,6] + delta[4,2,6] + epsilon[4,2] <= 0
 v[4,2,7]: 150 beta[4,2] + gamma[4,2,7] + delta[4,2,7] + epsilon[4,2] <= 0
 v[4,2,8]: 392 beta[4,2] + gamma[4,2,8] + delta[4,2,8] + epsilon[4,2] <= 0
 v[4,2,9]: 930 beta[4,2] + gamma[4,2,9] + delta[4,2,9] + epsilon[4,2] <= 0
 v[4,3,0]: 857 beta[4,3] + gamma[4,3,0] + delta[4,3,0] + epsilon[4,3] <= 0
 v[4,3,1]: 577 beta[4,3] + gamma[4,3,1] + delta[4,3,1] + epsilon[4,3] <= 0
 v[4,3,2]: 904 beta[4,3] + gamma[4,3,2] + delta[4,3,2] + epsilon[4,3] <= 0
 v[4,3,3]: gamma[4,3,3] + delta[4,3,3] + epsilon[4,3] <= 0
 v[4,3,4]: 881 beta[4,3] + gamma[4,3,4] + delta[4,3,4] + epsilon[4,3] <= 0
 v[4,3,5]: 984 beta[4,3] + gamma[4,3,5] + delta[4,3,5] + epsilon[4,3] <= 0
 v[4,3,6]: 856 beta[4,3] + gamma[4,3,6] + delta[4,3,6] + epsilon[4,3] <= 0
 v[4,3,7]: 520 beta[4,3] + gamma[4,3,7] + delta[4,3,7] + epsilon[4,3] <= 0
 v[4,3,8]: 476 beta[4,3] + gamma[4,3,8] + delta[4,3,8] + epsilon[4,3] <= 0
 v[4,3,9]: 1052 beta[4,3] + gamma[4,3,9] + delta[4,3,9] + epsilon[4,3]
   <= 0
 v[4,4,0]: 595 beta[4,4] + gamma[4,4,0] + delta[4,4,0] + epsilon[4,4] <= 0
 v[4,4,1]: 846 beta[4,4] + gamma[4,4,1] + delta[4,4,1] + epsilon[4,4] <= 0
 v[4,4,2]: 219 beta[4,4] + gamma[4,4,2] + delta[4,4,2] + epsilon[4,4] <= 0
 v[4,4,3]: 643 beta[4,4] + gamma[4,4,3] + delta[4,4,3] + epsilon[4,4] <= 0
 v[4,4,4]: gamma[4,4,4] + delta[4,4,4] + epsilon[4,4] <= 0
 v[4,4,5]: 111 beta[4,4] + gamma[4,4,5] + delta[4,4,5] + epsilon[4,4] <= 0
 v[4,4,6]: 511 beta[4,4] + gamma[4,4,6] + delta[4,4,6] + epsilon[4,4] <= 0
 v[4,4,7]: 990 beta[4,4] + gamma[4,4,7] + delta[4,4,7] + epsilon[4,4] <= 0
 v[4,4,8]: 785 beta[4,4] + gamma[4,4,8] + delta[4,4,8] + epsilon[4,4] <= 0
 v[4,4,9]: 1037 beta[4,4] + gamma[4,4,9] + delta[4,4,9] + epsilon[4,4]
   <= 0
 v[4,5,0]: 978 beta[4,5] + gamma[4,5,0] + delta[4,5,0] + epsilon[4,5] <= 0
 v[4,5,1]: 1097 beta[4,5] + gamma[4,5,1] + delta[4,5,1] + epsilon[4,5]
   <= 0
 v[4,5,2]: 846 beta[4,5] + gamma[4,5,2] + delta[4,5,2] + epsilon[4,5] <= 0
 v[4,5,3]: 740 beta[4,5] + gamma[4,5,3] + delta[4,5,3] + epsilon[4,5] <= 0
 v[4,5,4]: 999 beta[4,5] + gamma[4,5,4] + delta[4,5,4] + epsilon[4,5] <= 0
 v[4,5,5]: gamma[4,5,5] + delta[4,5,5] + epsilon[4,5] <= 0
 v[4,5,6]: 644 beta[4,5] + gamma[4,5,6] + delta[4,5,6] + epsilon[4,5] <= 0
 v[4,5,7]: 618 beta[4,5] + gamma[4,5,7] + delta[4,5,7] + epsilon[4,5] <= 0
 v[4,5,8]: 1099 beta[4,5] + gamma[4,5,8] + delta[4,5,8] + epsilon[4,5]
   <= 0
 v[4,5,9]: 226 beta[4,5] + gamma[4,5,9] + delta[4,5,9] + epsilon[4,5] <= 0
 v[4,6,0]: 379 beta[4,6] + gamma[4,6,0] + delta[4,6,0] + epsilon[4,6] <= 0
 v[4,6,1]: 898 beta[4,6] + gamma[4,6,1] + delta[4,6,1] + epsilon[4,6] <= 0
 v[4,6,2]: 242 beta[4,6] + gamma[4,6,2] + delta[4,6,2] + epsilon[4,6] <= 0
 v[4,6,3]: 674 beta[4,6] + gamma[4,6,3] + delta[4,6,3] + epsilon[4,6] <= 0
 v[4,6,4]: 827 beta[4,6] + gamma[4,6,4] + delta[4,6,4] + epsilon[4,6] <= 0
 v[4,6,5]: 1076 beta[4,6] + gamma[4,6,5] + delta[4,6,5] + epsilon[4,6]
   <= 0
 v[4,6,6]: gamma[4,6,6] + delta[4,6,6] + epsilon[4,6] <= 0
 v[4,6,7]: 423 beta[4,6] + gamma[4,6,7] + delta[4,6,7] + epsilon[4,6] <= 0
 v[4,6,8]: 200 beta[4,6] + gamma[4,6,8] + delta[4,6,8] + epsilon[4,6] <= 0
 v[4,6,9]: 522 beta[4,6] + gamma[4,6,9] + delta[4,6,9] + epsilon[4,6] <= 0
 v[4,7,0]: 912 beta[4,7] + gamma[4,7,0] + delta[4,7,0] + epsilon[4,7] <= 0
 v[4,7,1]: 546 beta[4,7] + gamma[4,7,1] + delta[4,7,1] + epsilon[4,7] <= 0
 v[4,7,2]: 212 beta[4,7] + gamma[4,7,2] + delta[4,7,2] + epsilon[4,7] <= 0
 v[4,7,3]: 728 beta[4,7] + gamma[4,7,3] + delta[4,7,3] + epsilon[4,7] <= 0
 v[4,7,4]: 658 beta[4,7] + gamma[4,7,4] + delta[4,7,4] + epsilon[4,7] <= 0
 v[4,7,5]: 371 beta[4,7] + gamma[4,7,5] + delta[4,7,5] + epsilon[4,7] <= 0
 v[4,7,6]: 531 beta[4,7] + gamma[4,7,6] + delta[4,7,6] + epsilon[4,7] <= 0
 v[4,7,7]: gamma[4,7,7] + delta[4,7,7] + epsilon[4,7] <= 0
 v[4,7,8]: 157 beta[4,7] + gamma[4,7,8] + delta[4,7,8] + epsilon[4,7] <= 0
 v[4,7,9]: 367 beta[4,7] + gamma[4,7,9] + delta[4,7,9] + epsilon[4,7] <= 0
 v[4,8,0]: 451 beta[4,8] + gamma[4,8,0] + delta[4,8,0] + epsilon[4,8] <= 0
 v[4,8,1]: 770 beta[4,8] + gamma[4,8,1] + delta[4,8,1] + epsilon[4,8] <= 0
 v[4,8,2]: 540 beta[4,8] + gamma[4,8,2] + delta[4,8,2] + epsilon[4,8] <= 0
 v[4,8,3]: 1072 beta[4,8] + gamma[4,8,3] + delta[4,8,3] + epsilon[4,8]
   <= 0
 v[4,8,4]: 846 beta[4,8] + gamma[4,8,4] + delta[4,8,4] + epsilon[4,8] <= 0
 v[4,8,5]: 759 beta[4,8] + gamma[4,8,5] + delta[4,8,5] + epsilon[4,8] <= 0
 v[4,8,6]: 825 beta[4,8] + gamma[4,8,6] + delta[4,8,6] + epsilon[4,8] <= 0
 v[4,8,7]: 1063 beta[4,8] + gamma[4,8,7] + delta[4,8,7] + epsilon[4,8]
   <= 0
 v[4,8,8]: gamma[4,8,8] + delta[4,8,8] + epsilon[4,8] <= 0
 v[4,8,9]: 816 beta[4,8] + gamma[4,8,9] + delta[4,8,9] + epsilon[4,8] <= 0
 v[4,9,0]: 1052 beta[4,9] + gamma[4,9,0] + delta[4,9,0] + epsilon[4,9]
   <= 0
 v[4,9,1]: 189 beta[4,9] + gamma[4,9,1] + delta[4,9,1] + epsilon[4,9] <= 0
 v[4,9,2]: 888 beta[4,9] + gamma[4,9,2] + delta[4,9,2] + epsilon[4,9] <= 0
 v[4,9,3]: 906 beta[4,9] + gamma[4,9,3] + delta[4,9,3] + epsilon[4,9] <= 0
 v[4,9,4]: 321 beta[4,9] + gamma[4,9,4] + delta[4,9,4] + epsilon[4,9] <= 0
 v[4,9,5]: 812 beta[4,9] + gamma[4,9,5] + delta[4,9,5] + epsilon[4,9] <= 0
 v[4,9,6]: 846 beta[4,9] + gamma[4,9,6] + delta[4,9,6] + epsilon[4,9] <= 0
 v[4,9,7]: 588 beta[4,9] + gamma[4,9,7] + delta[4,9,7] + epsilon[4,9] <= 0
 v[4,9,8]: 417 beta[4,9] + gamma[4,9,8] + delta[4,9,8] + epsilon[4,9] <= 0
 v[4,9,9]: gamma[4,9,9] + delta[4,9,9] + epsilon[4,9] <= 0
 v[5,0,0]: gamma[5,0,0] + delta[5,0,0] + epsilon[5,0] <= 0
 v[5,0,1]: 1058 beta[5,0] + gamma[5,0,1] + delta[5,0,1] + epsilon[5,0]
   <= 0
 v[5,0,2]: 544 beta[5,0] + gamma[5,0,2] + delta[5,0,2] + epsilon[5,0] <= 0
 v[5,0,3]: 714 beta[5,0] + gamma[5,0,3] + delta[5,0,3] + epsilon[5,0] <= 0
 v[5,0,4]: 1037 beta[5,0] + gamma[5,0,4] + delta[5,0,4] + epsilon[5,0]
   <= 0
 v[5,0,5]: 734 beta[5,0] + gamma[5,0,5] + delta[5,0,5] + epsilon[5,0] <= 0
 v[5,0,6]: 536 beta[5,0] + gamma[5,0,6] + delta[5,0,6] + epsilon[5,0] <= 0
 v[5,0,7]: 934 beta[5,0] + gamma[5,0,7] + delta[5,0,7] + epsilon[5,0] <= 0
 v[5,0,8]: 678 beta[5,0] + gamma[5,0,8] + delta[5,0,8] + epsilon[5,0] <= 0
 v[5,0,9]: 1009 beta[5,0] + gamma[5,0,9] + delta[5,0,9] + epsilon[5,0]
   <= 0
 v[5,1,0]: 944 beta[5,1] + gamma[5,1,0] + delta[5,1,0] + epsilon[5,1] <= 0
 v[5,1,1]: gamma[5,1,1] + delta[5,1,1] + epsilon[5,1] <= 0
 v[5,1,2]: 1053 beta[5,1] + gamma[5,1,2] + delta[5,1,2] + epsilon[5,1]
   <= 0
 v[5,1,3]: 944 beta[5,1] + gamma[5,1,3] + delta[5,1,3] + epsilon[5,1] <= 0
 v[5,1,4]: 646 beta[5,1] + gamma[5,1,4] + delta[5,1,4] + epsilon[5,1] <= 0
 v[5,1,5]: 354 beta[5,1] + gamma[5,1,5] + delta[5,1,5] + epsilon[5,1] <= 0
 v[5,1,6]: 281 beta[5,1] + gamma[5,1,6] + delta[5,1,6] + epsilon[5,1] <= 0
 v[5,1,7]: 786 beta[5,1] + gamma[5,1,7] + delta[5,1,7] + epsilon[5,1] <= 0
 v[5,1,8]: 1060 beta[5,1] + gamma[5,1,8] + delta[5,1,8] + epsilon[5,1]
   <= 0
 v[5,1,9]: 615 beta[5,1] + gamma[5,1,9] + delta[5,1,9] + epsilon[5,1] <= 0
 v[5,2,0]: 803 beta[5,2] + gamma[5,2,0] + delta[5,2,0] + epsilon[5,2] <= 0
 v[5,2,1]: 522 beta[5,2] + gamma[5,2,1] + delta[5,2,1] + epsilon[5,2] <= 0
 v[5,2,2]: gamma[5,2,2] + delta[5,2,2] + epsilon[5,2] <= 0
 v[5,2,3]: 967 beta[5,2] + gamma[5,2,3] + delta[5,2,3] + epsilon[5,2] <= 0
 v[5,2,4]: 962 beta[5,2] + gamma[5,2,4] + delta[5,2,4] + epsilon[5,2] <= 0
 v[5,2,5]: 720 beta[5,2] + gamma[5,2,5] + delta[5,2,5] + epsilon[5,2] <= 0
 v[5,2,6]: 341 beta[5,2] + gamma[5,2,6] + delta[5,2,6] + epsilon[5,2] <= 0
 v[5,2,7]: 434 beta[5,2] + gamma[5,2,7] + delta[5,2,7] + epsilon[5,2] <= 0
 v[5,2,8]: 271 beta[5,2] + gamma[5,2,8] + delta[5,2,8] + epsilon[5,2] <= 0
 v[5,2,9]: 523 beta[5,2] + gamma[5,2,9] + delta[5,2,9] + epsilon[5,2] <= 0
 v[5,3,0]: 717 beta[5,3] + gamma[5,3,0] + delta[5,3,0] + epsilon[5,3] <= 0
 v[5,3,1]: 694 beta[5,3] + gamma[5,3,1] + delta[5,3,1] + epsilon[5,3] <= 0
 v[5,3,2]: 1059 beta[5,3] + gamma[5,3,2] + delta[5,3,2] + epsilon[5,3]
   <= 0
 v[5,3,3]: gamma[5,3,3] + delta[5,3,3] + epsilon[5,3] <= 0
 v[5,3,4]: 409 beta[5,3] + gamma[5,3,4] + delta[5,3,4] + epsilon[5,3] <= 0
 v[5,3,5]: 931 beta[5,3] + gamma[5,3,5] + delta[5,3,5] + epsilon[5,3] <= 0
 v[5,3,6]: 143 beta[5,3] + gamma[5,3,6] + delta[5,3,6] + epsilon[5,3] <= 0
 v[5,3,7]: 202 beta[5,3] + gamma[5,3,7] + delta[5,3,7] + epsilon[5,3] <= 0
 v[5,3,8]: 954 beta[5,3] + gamma[5,3,8] + delta[5,3,8] + epsilon[5,3] <= 0
 v[5,3,9]: 882 beta[5,3] + gamma[5,3,9] + delta[5,3,9] + epsilon[5,3] <= 0
 v[5,4,0]: 482 beta[5,4] + gamma[5,4,0] + delta[5,4,0] + epsilon[5,4] <= 0
 v[5,4,1]: 839 beta[5,4] + gamma[5,4,1] + delta[5,4,1] + epsilon[5,4] <= 0
 v[5,4,2]: 256 beta[5,4] + gamma[5,4,2] + delta[5,4,2] + epsilon[5,4] <= 0
 v[5,4,3]: 926 beta[5,4] + gamma[5,4,3] + delta[5,4,3] + epsilon[5,4] <= 0
 v[5,4,4]: gamma[5,4,4] + delta[5,4,4] + epsilon[5,4] <= 0
 v[5,4,5]: 1026 beta[5,4] + gamma[5,4,5] + delta[5,4,5] + epsilon[5,4]
   <= 0
 v[5,4,6]: 239 beta[5,4] + gamma[5,4,6] + delta[5,4,6] + epsilon[5,4] <= 0
 v[5,4,7]: 914 beta[5,4] + gamma[5,4,7] + delta[5,4,7] + epsilon[5,4] <= 0
 v[5,4,8]: 681 beta[5,4] + gamma[5,4,8] + delta[5,4,8] + epsilon[5,4] <= 0
 v[5,4,9]: 172 beta[5,4] + gamma[5,4,9] + delta[5,4,9] + epsilon[5,4] <= 0
 v[5,5,0]: 702 beta[5,5] + gamma[5,5,0] + delta[5,5,0] + epsilon[5,5] <= 0
 v[5,5,1]: 929 beta[5,5] + gamma[5,5,1] + delta[5,5,1] + epsilon[5,5] <= 0
 v[5,5,2]: 498 beta[5,5] + gamma[5,5,2] + delta[5,5,2] + epsilon[5,5] <= 0
 v[5,5,3]: 473 beta[5,5] + gamma[5,5,3] + delta[5,5,3] + epsilon[5,5] <= 0
 v[5,5,4]: 1088 beta[5,5] + gamma[5,5,4] + delta[5,5,4] + epsilon[5,5]
   <= 0
 v[5,5,5]: gamma[5,5,5] + delta[5,5,5] + epsilon[5,5] <= 0
 v[5,5,6]: 215 beta[5,5] + gamma[5,5,6] + delta[5,5,6] + epsilon[5,5] <= 0
 v[5,5,7]: 371 beta[5,5] + gamma[5,5,7] + delta[5,5,7] + epsilon[5,5] <= 0
 v[5,5,8]: 660 beta[5,5] + gamma[5,5,8] + delta[5,5,8] + epsilon[5,5] <= 0
 v[5,5,9]: 665 beta[5,5] + gamma[5,5,9] + delta[5,5,9] + epsilon[5,5] <= 0
 v[5,6,0]: 398 beta[5,6] + gamma[5,6,0] + delta[5,6,0] + epsilon[5,6] <= 0
 v[5,6,1]: 576 beta[5,6] + gamma[5,6,1] + delta[5,6,1] + epsilon[5,6] <= 0
 v[5,6,2]: 866 beta[5,6] + gamma[5,6,2] + delta[5,6,2] + epsilon[5,6] <= 0
 v[5,6,3]: 118 beta[5,6] + gamma[5,6,3] + delta[5,6,3] + epsilon[5,6] <= 0
 v[5,6,4]: 305 beta[5,6] + gamma[5,6,4] + delta[5,6,4] + epsilon[5,6] <= 0
 v[5,6,5]: 645 beta[5,6] + gamma[5,6,5] + delta[5,6,5] + epsilon[5,6] <= 0
 v[5,6,6]: gamma[5,6,6] + delta[5,6,6] + epsilon[5,6] <= 0
 v[5,6,7]: 1043 beta[5,6] + gamma[5,6,7] + delta[5,6,7] + epsilon[5,6]
   <= 0
 v[5,6,8]: 210 beta[5,6] + gamma[5,6,8] + delta[5,6,8] + epsilon[5,6] <= 0
 v[5,6,9]: 747 beta[5,6] + gamma[5,6,9] + delta[5,6,9] + epsilon[5,6] <= 0
 v[5,7,0]: 582 beta[5,7] + gamma[5,7,0] + delta[5,7,0] + epsilon[5,7] <= 0
 v[5,7,1]: 887 beta[5,7] + gamma[5,7,1] + delta[5,7,1] + epsilon[5,7] <= 0
 v[5,7,2]: 1016 beta[5,7] + gamma[5,7,2] + delta[5,7,2] + epsilon[5,7]
   <= 0
 v[5,7,3]: 1002 beta[5,7] + gamma[5,7,3] + delta[5,7,3] + epsilon[5,7]
   <= 0
 v[5,7,4]: 451 beta[5,7] + gamma[5,7,4] + delta[5,7,4] + epsilon[5,7] <= 0
 v[5,7,5]: 734 beta[5,7] + gamma[5,7,5] + delta[5,7,5] + epsilon[5,7] <= 0
 v[5,7,6]: 433 beta[5,7] + gamma[5,7,6] + delta[5,7,6] + epsilon[5,7] <= 0
 v[5,7,7]: gamma[5,7,7] + delta[5,7,7] + epsilon[5,7] <= 0
 v[5,7,8]: 373 beta[5,7] + gamma[5,7,8] + delta[5,7,8] + epsilon[5,7] <= 0
 v[5,7,9]: 157 beta[5,7] + gamma[5,7,9] + delta[5,7,9] + epsilon[5,7] <= 0
 v[5,8,0]: 743 beta[5,8] + gamma[5,8,0] + delta[5,8,0] + epsilon[5,8] <= 0
 v[5,8,1]: 1017 beta[5,8] + gamma[5,8,1] + delta[5,8,1] + epsilon[5,8]
   <= 0
 v[5,8,2]: 520 beta[5,8] + gamma[5,8,2] + delta[5,8,2] + epsilon[5,8] <= 0
 v[5,8,3]: 501 beta[5,8] + gamma[5,8,3] + delta[5,8,3] + epsilon[5,8] <= 0
 v[5,8,4]: 1025 beta[5,8] + gamma[5,8,4] + delta[5,8,4] + epsilon[5,8]
   <= 0
 v[5,8,5]: 321 beta[5,8] + gamma[5,8,5] + delta[5,8,5] + epsilon[5,8] <= 0
 v[5,8,6]: 1036 beta[5,8] + gamma[5,8,6] + delta[5,8,6] + epsilon[5,8]
   <= 0
 v[5,8,7]: 363 beta[5,8] + gamma[5,8,7] + delta[5,8,7] + epsilon[5,8] <= 0
 v[5,8,8]: gamma[5,8,8] + delta[5,8,8] + epsilon[5,8] <= 0
 v[5,8,9]: 205 beta[5,8] + gamma[5,8,9] + delta[5,8,9] + epsilon[5,8] <= 0
 v[5,9,0]: 342 beta[5,9] + gamma[5,9,0] + delta[5,9,0] + epsilon[5,9] <= 0
 v[5,9,1]: 516 beta[5,9] + gamma[5,9,1] + delta[5,9,1] + epsilon[5,9] <= 0
 v[5,9,2]: 733 beta[5,9] + gamma[5,9,2] + delta[5,9,2] + epsilon[5,9] <= 0
 v[5,9,3]: 155 beta[5,9] + gamma[5,9,3] + delta[5,9,3] + epsilon[5,9] <= 0
 v[5,9,4]: 630 beta[5,9] + gamma[5,9,4] + delta[5,9,4] + epsilon[5,9] <= 0
 v[5,9,5]: 637 beta[5,9] + gamma[5,9,5] + delta[5,9,5] + epsilon[5,9] <= 0
 v[5,9,6]: 440 beta[5,9] + gamma[5,9,6] + delta[5,9,6] + epsilon[5,9] <= 0
 v[5,9,7]: 411 beta[5,9] + gamma[5,9,7] + delta[5,9,7] + epsilon[5,9] <= 0
 v[5,9,8]: 642 beta[5,9] + gamma[5,9,8] + delta[5,9,8] + epsilon[5,9] <= 0
 v[5,9,9]: gamma[5,9,9] + delta[5,9,9] + epsilon[5,9] <= 0
 v[6,0,0]: gamma[6,0,0] + delta[6,0,0] + epsilon[6,0] <= 0
 v[6,0,1]: 224 beta[6,0] + gamma[6,0,1] + delta[6,0,1] + epsilon[6,0] <= 0
 v[6,0,2]: 337 beta[6,0] + gamma[6,0,2] + delta[6,0,2] + epsilon[6,0] <= 0
 v[6,0,3]: 827 beta[6,0] + gamma[6,0,3] + delta[6,0,3] + epsilon[6,0] <= 0
 v[6,0,4]: 1024 beta[6,0] + gamma[6,0,4] + delta[6,0,4] + epsilon[6,0]
   <= 0
 v[6,0,5]: 289 beta[6,0] + gamma[6,0,5] + delta[6,0,5] + epsilon[6,0] <= 0
 v[6,0,6]: 996 beta[6,0] + gamma[6,0,6] + delta[6,0,6] + epsilon[6,0] <= 0
 v[6,0,7]: 939 beta[6,0] + gamma[6,0,7] + delta[6,0,7] + epsilon[6,0] <= 0
 v[6,0,8]: 994 beta[6,0] + gamma[6,0,8] + delta[6,0,8] + epsilon[6,0] <= 0
 v[6,0,9]: 649 beta[6,0] + gamma[6,0,9] + delta[6,0,9] + epsilon[6,0] <= 0
 v[6,1,0]: 1006 beta[6,1] + gamma[6,1,0] + delta[6,1,0] + epsilon[6,1]
   <= 0
 v[6,1,1]: gamma[6,1,1] + delta[6,1,1] + epsilon[6,1] <= 0
 v[6,1,2]: 269 beta[6,1] + gamma[6,1,2] + delta[6,1,2] + epsilon[6,1] <= 0
 v[6,1,3]: 911 beta[6,1] + gamma[6,1,3] + delta[6,1,3] + epsilon[6,1] <= 0
 v[6,1,4]: 882 beta[6,1] + gamma[6,1,4] + delta[6,1,4] + epsilon[6,1] <= 0
 v[6,1,5]: 738 beta[6,1] + gamma[6,1,5] + delta[6,1,5] + epsilon[6,1] <= 0
 v[6,1,6]: 995 beta[6,1] + gamma[6,1,6] + delta[6,1,6] + epsilon[6,1] <= 0
 v[6,1,7]: 158 beta[6,1] + gamma[6,1,7] + delta[6,1,7] + epsilon[6,1] <= 0
 v[6,1,8]: 676 beta[6,1] + gamma[6,1,8] + delta[6,1,8] + epsilon[6,1] <= 0
 v[6,1,9]: 468 beta[6,1] + gamma[6,1,9] + delta[6,1,9] + epsilon[6,1] <= 0
 v[6,2,0]: 274 beta[6,2] + gamma[6,2,0] + delta[6,2,0] + epsilon[6,2] <= 0
 v[6,2,1]: 604 beta[6,2] + gamma[6,2,1] + delta[6,2,1] + epsilon[6,2] <= 0
 v[6,2,2]: gamma[6,2,2] + delta[6,2,2] + epsilon[6,2] <= 0
 v[6,2,3]: 819 beta[6,2] + gamma[6,2,3] + delta[6,2,3] + epsilon[6,2] <= 0
 v[6,2,4]: 1072 beta[6,2] + gamma[6,2,4] + delta[6,2,4] + epsilon[6,2]
   <= 0
 v[6,2,5]: 129 beta[6,2] + gamma[6,2,5] + delta[6,2,5] + epsilon[6,2] <= 0
 v[6,2,6]: 1042 beta[6,2] + gamma[6,2,6] + delta[6,2,6] + epsilon[6,2]
   <= 0
 v[6,2,7]: 709 beta[6,2] + gamma[6,2,7] + delta[6,2,7] + epsilon[6,2] <= 0
 v[6,2,8]: 762 beta[6,2] + gamma[6,2,8] + delta[6,2,8] + epsilon[6,2] <= 0
 v[6,2,9]: 930 beta[6,2] + gamma[6,2,9] + delta[6,2,9] + epsilon[6,2] <= 0
 v[6,3,0]: 635 beta[6,3] + gamma[6,3,0] + delta[6,3,0] + epsilon[6,3] <= 0
 v[6,3,1]: 965 beta[6,3] + gamma[6,3,1] + delta[6,3,1] + epsilon[6,3] <= 0
 v[6,3,2]: 669 beta[6,3] + gamma[6,3,2] + delta[6,3,2] + epsilon[6,3] <= 0
 v[6,3,3]: gamma[6,3,3] + delta[6,3,3] + epsilon[6,3] <= 0
 v[6,3,4]: 441 beta[6,3] + gamma[6,3,4] + delta[6,3,4] + epsilon[6,3] <= 0
 v[6,3,5]: 267 beta[6,3] + gamma[6,3,5] + delta[6,3,5] + epsilon[6,3] <= 0
 v[6,3,6]: 867 beta[6,3] + gamma[6,3,6] + delta[6,3,6] + epsilon[6,3] <= 0
 v[6,3,7]: 336 beta[6,3] + gamma[6,3,7] + delta[6,3,7] + epsilon[6,3] <= 0
 v[6,3,8]: 338 beta[6,3] + gamma[6,3,8] + delta[6,3,8] + epsilon[6,3] <= 0
 v[6,3,9]: 894 beta[6,3] + gamma[6,3,9] + delta[6,3,9] + epsilon[6,3] <= 0
 v[6,4,0]: 307 beta[6,4] + gamma[6,4,0] + delta[6,4,0] + epsilon[6,4] <= 0
 v[6,4,1]: 636 beta[6,4] + gamma[6,4,1] + delta[6,4,1] + epsilon[6,4] <= 0
 v[6,4,2]: 247 beta[6,4] + gamma[6,4,2] + delta[6,4,2] + epsilon[6,4] <= 0
 v[6,4,3]: 788 beta[6,4] + gamma[6,4,3] + delta[6,4,3] + epsilon[6,4] <= 0
 v[6,4,4]: gamma[6,4,4] + delta[6,4,4] + epsilon[6,4] <= 0
 v[6,4,5]: 1050 beta[6,4] + gamma[6,4,5] + delta[6,4,5] + epsilon[6,4]
   <= 0
 v[6,4,6]: 705 beta[6,4] + gamma[6,4,6] + delta[6,4,6] + epsilon[6,4] <= 0
 v[6,4,7]: 673 beta[6,4] + gamma[6,4,7] + delta[6,4,7] + epsilon[6,4] <= 0
 v[6,4,8]: 485 beta[6,4] + gamma[6,4,8] + delta[6,4,8] + epsilon[6,4] <= 0
 v[6,4,9]: 133 beta[6,4] + gamma[6,4,9] + delta[6,4,9] + epsilon[6,4] <= 0
 v[6,5,0]: 725 beta[6,5] + gamma[6,5,0] + delta[6,5,0] + epsilon[6,5] <= 0
 v[6,5,1]: 876 beta[6,5] + gamma[6,5,1] + delta[6,5,1] + epsilon[6,5] <= 0
 v[6,5,2]: 546 beta[6,5] + gamma[6,5,2] + delta[6,5,2] + epsilon[6,5] <= 0
 v[6,5,3]: 1026 beta[6,5] + gamma[6,5,3] + delta[6,5,3] + epsilon[6,5]
   <= 0
 v[6,5,4]: 133 beta[6,5] + gamma[6,5,4] + delta[6,5,4] + epsilon[6,5] <= 0
 v[6,5,5]: gamma[6,5,5] + delta[6,5,5] + epsilon[6,5] <= 0
 v[6,5,6]: 423 beta[6,5] + gamma[6,5,6] + delta[6,5,6] + epsilon[6,5] <= 0
 v[6,5,7]: 227 beta[6,5] + gamma[6,5,7] + delta[6,5,7] + epsilon[6,5] <= 0
 v[6,5,8]: 1100 beta[6,5] + gamma[6,5,8] + delta[6,5,8] + epsilon[6,5]
   <= 0
 v[6,5,9]: 1052 beta[6,5] + gamma[6,5,9] + delta[6,5,9] + epsilon[6,5]
   <= 0
 v[6,6,0]: 871 beta[6,6] + gamma[6,6,0] + delta[6,6,0] + epsilon[6,6] <= 0
 v[6,6,1]: 147 beta[6,6] + gamma[6,6,1] + delta[6,6,1] + epsilon[6,6] <= 0
 v[6,6,2]: 294 beta[6,6] + gamma[6,6,2] + delta[6,6,2] + epsilon[6,6] <= 0
 v[6,6,3]: 157 beta[6,6] + gamma[6,6,3] + delta[6,6,3] + epsilon[6,6] <= 0
 v[6,6,4]: 684 beta[6,6] + gamma[6,6,4] + delta[6,6,4] + epsilon[6,6] <= 0
 v[6,6,5]: 163 beta[6,6] + gamma[6,6,5] + delta[6,6,5] + epsilon[6,6] <= 0
 v[6,6,6]: gamma[6,6,6] + delta[6,6,6] + epsilon[6,6] <= 0
 v[6,6,7]: 269 beta[6,6] + gamma[6,6,7] + delta[6,6,7] + epsilon[6,6] <= 0
 v[6,6,8]: 532 beta[6,6] + gamma[6,6,8] + delta[6,6,8] + epsilon[6,6] <= 0
 v[6,6,9]: 102 beta[6,6] + gamma[6,6,9] + delta[6,6,9] + epsilon[6,6] <= 0
 v[6,7,0]: 265 beta[6,7] + gamma[6,7,0] + delta[6,7,0] + epsilon[6,7] <= 0
 v[6,7,1]: 124 beta[6,7] + gamma[6,7,1] + delta[6,7,1] + epsilon[6,7] <= 0
 v[6,7,2]: 966 beta[6,7] + gamma[6,7,2] + delta[6,7,2] + epsilon[6,7] <= 0
 v[6,7,3]: 640 beta[6,7] + gamma[6,7,3] + delta[6,7,3] + epsilon[6,7] <= 0
 v[6,7,4]: 450 beta[6,7] + gamma[6,7,4] + delta[6,7,4] + epsilon[6,7] <= 0
 v[6,7,5]: 604 beta[6,7] + gamma[6,7,5] + delta[6,7,5] + epsilon[6,7] <= 0
 v[6,7,6]: 404 beta[6,7] + gamma[6,7,6] + delta[6,7,6] + epsilon[6,7] <= 0
 v[6,7,7]: gamma[6,7,7] + delta[6,7,7] + epsilon[6,7] <= 0
 v[6,7,8]: 755 beta[6,7] + gamma[6,7,8] + delta[6,7,8] + epsilon[6,7] <= 0
 v[6,7,9]: 218 beta[6,7] + gamma[6,7,9] + delta[6,7,9] + epsilon[6,7] <= 0
 v[6,8,0]: 531 beta[6,8] + gamma[6,8,0] + delta[6,8,0] + epsilon[6,8] <= 0
 v[6,8,1]: 681 beta[6,8] + gamma[6,8,1] + delta[6,8,1] + epsilon[6,8] <= 0
 v[6,8,2]: 210 beta[6,8] + gamma[6,8,2] + delta[6,8,2] + epsilon[6,8] <= 0
 v[6,8,3]: 106 beta[6,8] + gamma[6,8,3] + delta[6,8,3] + epsilon[6,8] <= 0
 v[6,8,4]: 145 beta[6,8] + gamma[6,8,4] + delta[6,8,4] + epsilon[6,8] <= 0
 v[6,8,5]: 407 beta[6,8] + gamma[6,8,5] + delta[6,8,5] + epsilon[6,8] <= 0
 v[6,8,6]: 272 beta[6,8] + gamma[6,8,6] + delta[6,8,6] + epsilon[6,8] <= 0
 v[6,8,7]: 345 beta[6,8] + gamma[6,8,7] + delta[6,8,7] + epsilon[6,8] <= 0
 v[6,8,8]: gamma[6,8,8] + delta[6,8,8] + epsilon[6,8] <= 0
 v[6,8,9]: 301 beta[6,8] + gamma[6,8,9] + delta[6,8,9] + epsilon[6,8] <= 0
 v[6,9,0]: 275 beta[6,9] + gamma[6,9,0] + delta[6,9,0] + epsilon[6,9] <= 0
 v[6,9,1]: 1082 beta[6,9] + gamma[6,9,1] + delta[6,9,1] + epsilon[6,9]
   <= 0
 v[6,9,2]: 1033 beta[6,9] + gamma[6,9,2] + delta[6,9,2] + epsilon[6,9]
   <= 0
 v[6,9,3]: 192 beta[6,9] + gamma[6,9,3] + delta[6,9,3] + epsilon[6,9] <= 0
 v[6,9,4]: 762 beta[6,9] + gamma[6,9,4] + delta[6,9,4] + epsilon[6,9] <= 0
 v[6,9,5]: 347 beta[6,9] + gamma[6,9,5] + delta[6,9,5] + epsilon[6,9] <= 0
 v[6,9,6]: 798 beta[6,9] + gamma[6,9,6] + delta[6,9,6] + epsilon[6,9] <= 0
 v[6,9,7]: 687 beta[6,9] + gamma[6,9,7] + delta[6,9,7] + epsilon[6,9] <= 0
 v[6,9,8]: 1011 beta[6,9] + gamma[6,9,8] + delta[6,9,8] + epsilon[6,9]
   <= 0
 v[6,9,9]: gamma[6,9,9] + delta[6,9,9] + epsilon[6,9] <= 0
 v[7,0,0]: gamma[7,0,0] + delta[7,0,0] + epsilon[7,0] <= 0
 v[7,0,1]: 303 beta[7,0] + gamma[7,0,1] + delta[7,0,1] + epsilon[7,0] <= 0
 v[7,0,2]: 561 beta[7,0] + gamma[7,0,2] + delta[7,0,2] + epsilon[7,0] <= 0
 v[7,0,3]: 824 beta[7,0] + gamma[7,0,3] + delta[7,0,3] + epsilon[7,0] <= 0
 v[7,0,4]: 947 beta[7,0] + gamma[7,0,4] + delta[7,0,4] + epsilon[7,0] <= 0
 v[7,0,5]: 342 beta[7,0] + gamma[7,0,5] + delta[7,0,5] + epsilon[7,0] <= 0
 v[7,0,6]: 386 beta[7,0] + gamma[7,0,6] + delta[7,0,6] + epsilon[7,0] <= 0
 v[7,0,7]: 656 beta[7,0] + gamma[7,0,7] + delta[7,0,7] + epsilon[7,0] <= 0
 v[7,0,8]: 619 beta[7,0] + gamma[7,0,8] + delta[7,0,8] + epsilon[7,0] <= 0
 v[7,0,9]: 161 beta[7,0] + gamma[7,0,9] + delta[7,0,9] + epsilon[7,0] <= 0
 v[7,1,0]: 584 beta[7,1] + gamma[7,1,0] + delta[7,1,0] + epsilon[7,1] <= 0
 v[7,1,1]: gamma[7,1,1] + delta[7,1,1] + epsilon[7,1] <= 0
 v[7,1,2]: 170 beta[7,1] + gamma[7,1,2] + delta[7,1,2] + epsilon[7,1] <= 0
 v[7,1,3]: 103 beta[7,1] + gamma[7,1,3] + delta[7,1,3] + epsilon[7,1] <= 0
 v[7,1,4]: 1006 beta[7,1] + gamma[7,1,4] + delta[7,1,4] + epsilon[7,1]
   <= 0
 v[7,1,5]: 308 beta[7,1] + gamma[7,1,5] + delta[7,1,5] + epsilon[7,1] <= 0
 v[7,1,6]: 631 beta[7,1] + gamma[7,1,6] + delta[7,1,6] + epsilon[7,1] <= 0
 v[7,1,7]: 308 beta[7,1] + gamma[7,1,7] + delta[7,1,7] + epsilon[7,1] <= 0
 v[7,1,8]: 185 beta[7,1] + gamma[7,1,8] + delta[7,1,8] + epsilon[7,1] <= 0
 v[7,1,9]: 1021 beta[7,1] + gamma[7,1,9] + delta[7,1,9] + epsilon[7,1]
   <= 0
 v[7,2,0]: 1002 beta[7,2] + gamma[7,2,0] + delta[7,2,0] + epsilon[7,2]
   <= 0
 v[7,2,1]: 809 beta[7,2] + gamma[7,2,1] + delta[7,2,1] + epsilon[7,2] <= 0
 v[7,2,2]: gamma[7,2,2] + delta[7,2,2] + epsilon[7,2] <= 0
 v[7,2,3]: 389 beta[7,2] + gamma[7,2,3] + delta[7,2,3] + epsilon[7,2] <= 0
 v[7,2,4]: 597 beta[7,2] + gamma[7,2,4] + delta[7,2,4] + epsilon[7,2] <= 0
 v[7,2,5]: 591 beta[7,2] + gamma[7,2,5] + delta[7,2,5] + epsilon[7,2] <= 0
 v[7,2,6]: 165 beta[7,2] + gamma[7,2,6] + delta[7,2,6] + epsilon[7,2] <= 0
 v[7,2,7]: 285 beta[7,2] + gamma[7,2,7] + delta[7,2,7] + epsilon[7,2] <= 0
 v[7,2,8]: 882 beta[7,2] + gamma[7,2,8] + delta[7,2,8] + epsilon[7,2] <= 0
 v[7,2,9]: 180 beta[7,2] + gamma[7,2,9] + delta[7,2,9] + epsilon[7,2] <= 0
 v[7,3,0]: 353 beta[7,3] + gamma[7,3,0] + delta[7,3,0] + epsilon[7,3] <= 0
 v[7,3,1]: 197 beta[7,3] + gamma[7,3,1] + delta[7,3,1] + epsilon[7,3] <= 0
 v[7,3,2]: 811 beta[7,3] + gamma[7,3,2] + delta[7,3,2] + epsilon[7,3] <= 0
 v[7,3,3]: gamma[7,3,3] + delta[7,3,3] + epsilon[7,3] <= 0
 v[7,3,4]: 714 beta[7,3] + gamma[7,3,4] + delta[7,3,4] + epsilon[7,3] <= 0
 v[7,3,5]: 121 beta[7,3] + gamma[7,3,5] + delta[7,3,5] + epsilon[7,3] <= 0
 v[7,3,6]: 176 beta[7,3] + gamma[7,3,6] + delta[7,3,6] + epsilon[7,3] <= 0
 v[7,3,7]: 820 beta[7,3] + gamma[7,3,7] + delta[7,3,7] + epsilon[7,3] <= 0
 v[7,3,8]: 465 beta[7,3] + gamma[7,3,8] + delta[7,3,8] + epsilon[7,3] <= 0
 v[7,3,9]: 748 beta[7,3] + gamma[7,3,9] + delta[7,3,9] + epsilon[7,3] <= 0
 v[7,4,0]: 838 beta[7,4] + gamma[7,4,0] + delta[7,4,0] + epsilon[7,4] <= 0
 v[7,4,1]: 179 beta[7,4] + gamma[7,4,1] + delta[7,4,1] + epsilon[7,4] <= 0
 v[7,4,2]: 418 beta[7,4] + gamma[7,4,2] + delta[7,4,2] + epsilon[7,4] <= 0
 v[7,4,3]: 1073 beta[7,4] + gamma[7,4,3] + delta[7,4,3] + epsilon[7,4]
   <= 0
 v[7,4,4]: gamma[7,4,4] + delta[7,4,4] + epsilon[7,4] <= 0
 v[7,4,5]: 1049 beta[7,4] + gamma[7,4,5] + delta[7,4,5] + epsilon[7,4]
   <= 0
 v[7,4,6]: 177 beta[7,4] + gamma[7,4,6] + delta[7,4,6] + epsilon[7,4] <= 0
 v[7,4,7]: 195 beta[7,4] + gamma[7,4,7] + delta[7,4,7] + epsilon[7,4] <= 0
 v[7,4,8]: 483 beta[7,4] + gamma[7,4,8] + delta[7,4,8] + epsilon[7,4] <= 0
 v[7,4,9]: 414 beta[7,4] + gamma[7,4,9] + delta[7,4,9] + epsilon[7,4] <= 0
 v[7,5,0]: 121 beta[7,5] + gamma[7,5,0] + delta[7,5,0] + epsilon[7,5] <= 0
 v[7,5,1]: 202 beta[7,5] + gamma[7,5,1] + delta[7,5,1] + epsilon[7,5] <= 0
 v[7,5,2]: 556 beta[7,5] + gamma[7,5,2] + delta[7,5,2] + epsilon[7,5] <= 0
 v[7,5,3]: 482 beta[7,5] + gamma[7,5,3] + delta[7,5,3] + epsilon[7,5] <= 0
 v[7,5,4]: 531 beta[7,5] + gamma[7,5,4] + delta[7,5,4] + epsilon[7,5] <= 0
 v[7,5,5]: gamma[7,5,5] + delta[7,5,5] + epsilon[7,5] <= 0
 v[7,5,6]: 139 beta[7,5] + gamma[7,5,6] + delta[7,5,6] + epsilon[7,5] <= 0
 v[7,5,7]: 705 beta[7,5] + gamma[7,5,7] + delta[7,5,7] + epsilon[7,5] <= 0
 v[7,5,8]: 1047 beta[7,5] + gamma[7,5,8] + delta[7,5,8] + epsilon[7,5]
   <= 0
 v[7,5,9]: 207 beta[7,5] + gamma[7,5,9] + delta[7,5,9] + epsilon[7,5] <= 0
 v[7,6,0]: 1045 beta[7,6] + gamma[7,6,0] + delta[7,6,0] + epsilon[7,6]
   <= 0
 v[7,6,1]: 1008 beta[7,6] + gamma[7,6,1] + delta[7,6,1] + epsilon[7,6]
   <= 0
 v[7,6,2]: 762 beta[7,6] + gamma[7,6,2] + delta[7,6,2] + epsilon[7,6] <= 0
 v[7,6,3]: 217 beta[7,6] + gamma[7,6,3] + delta[7,6,3] + epsilon[7,6] <= 0
 v[7,6,4]: 159 beta[7,6] + gamma[7,6,4] + delta[7,6,4] + epsilon[7,6] <= 0
 v[7,6,5]: 835 beta[7,6] + gamma[7,6,5] + delta[7,6,5] + epsilon[7,6] <= 0
 v[7,6,6]: gamma[7,6,6] + delta[7,6,6] + epsilon[7,6] <= 0
 v[7,6,7]: 652 beta[7,6] + gamma[7,6,7] + delta[7,6,7] + epsilon[7,6] <= 0
 v[7,6,8]: 593 beta[7,6] + gamma[7,6,8] + delta[7,6,8] + epsilon[7,6] <= 0
 v[7,6,9]: 720 beta[7,6] + gamma[7,6,9] + delta[7,6,9] + epsilon[7,6] <= 0
 v[7,7,0]: 595 beta[7,7] + gamma[7,7,0] + delta[7,7,0] + epsilon[7,7] <= 0
 v[7,7,1]: 688 beta[7,7] + gamma[7,7,1] + delta[7,7,1] + epsilon[7,7] <= 0
 v[7,7,2]: 1007 beta[7,7] + gamma[7,7,2] + delta[7,7,2] + epsilon[7,7]
   <= 0
 v[7,7,3]: 463 beta[7,7] + gamma[7,7,3] + delta[7,7,3] + epsilon[7,7] <= 0
 v[7,7,4]: 953 beta[7,7] + gamma[7,7,4] + delta[7,7,4] + epsilon[7,7] <= 0
 v[7,7,5]: 839 beta[7,7] + gamma[7,7,5] + delta[7,7,5] + epsilon[7,7] <= 0
 v[7,7,6]: 669 beta[7,7] + gamma[7,7,6] + delta[7,7,6] + epsilon[7,7] <= 0
 v[7,7,7]: gamma[7,7,7] + delta[7,7,7] + epsilon[7,7] <= 0
 v[7,7,8]: 887 beta[7,7] + gamma[7,7,8] + delta[7,7,8] + epsilon[7,7] <= 0
 v[7,7,9]: 896 beta[7,7] + gamma[7,7,9] + delta[7,7,9] + epsilon[7,7] <= 0
 v[7,8,0]: 477 beta[7,8] + gamma[7,8,0] + delta[7,8,0] + epsilon[7,8] <= 0
 v[7,8,1]: 139 beta[7,8] + gamma[7,8,1] + delta[7,8,1] + epsilon[7,8] <= 0
 v[7,8,2]: 368 beta[7,8] + gamma[7,8,2] + delta[7,8,2] + epsilon[7,8] <= 0
 v[7,8,3]: 152 beta[7,8] + gamma[7,8,3] + delta[7,8,3] + epsilon[7,8] <= 0
 v[7,8,4]: 769 beta[7,8] + gamma[7,8,4] + delta[7,8,4] + epsilon[7,8] <= 0
 v[7,8,5]: 792 beta[7,8] + gamma[7,8,5] + delta[7,8,5] + epsilon[7,8] <= 0
 v[7,8,6]: 661 beta[7,8] + gamma[7,8,6] + delta[7,8,6] + epsilon[7,8] <= 0
 v[7,8,7]: 249 beta[7,8] + gamma[7,8,7] + delta[7,8,7] + epsilon[7,8] <= 0
 v[7,8,8]: gamma[7,8,8] + delta[7,8,8] + epsilon[7,8] <= 0
 v[7,8,9]: 457 beta[7,8] + gamma[7,8,9] + delta[7,8,9] + epsilon[7,8] <= 0
 v[7,9,0]: 511 beta[7,9] + gamma[7,9,0] + delta[7,9,0] + epsilon[7,9] <= 0
 v[7,9,1]: 1013 beta[7,9] + gamma[7,9,1] + delta[7,9,1] + epsilon[7,9]
   <= 0
 v[7,9,2]: 371 beta[7,9] + gamma[7,9,2] + delta[7,9,2] + epsilon[7,9] <= 0
 v[7,9,3]: 240 beta[7,9] + gamma[7,9,3] + delta[7,9,3] + epsilon[7,9] <= 0
 v[7,9,4]: 987 beta[7,9] + gamma[7,9,4] + delta[7,9,4] + epsilon[7,9] <= 0
 v[7,9,5]: 221 beta[7,9] + gamma[7,9,5] + delta[7,9,5] + epsilon[7,9] <= 0
 v[7,9,6]: 215 beta[7,9] + gamma[7,9,6] + delta[7,9,6] + epsilon[7,9] <= 0
 v[7,9,7]: 452 beta[7,9] + gamma[7,9,7] + delta[7,9,7] + epsilon[7,9] <= 0
 v[7,9,8]: 759 beta[7,9] + gamma[7,9,8] + delta[7,9,8] + epsilon[7,9] <= 0
 v[7,9,9]: gamma[7,9,9] + delta[7,9,9] + epsilon[7,9] <= 0
 v[8,0,0]: gamma[8,0,0] + delta[8,0,0] + epsilon[8,0] <= 0
 v[8,0,1]: 160 beta[8,0] + gamma[8,0,1] + delta[8,0,1] + epsilon[8,0] <= 0
 v[8,0,2]: 636 beta[8,0] + gamma[8,0,2] + delta[8,0,2] + epsilon[8,0] <= 0
 v[8,0,3]: 943 beta[8,0] + gamma[8,0,3] + delta[8,0,3] + epsilon[8,0] <= 0
 v[8,0,4]: 162 beta[8,0] + gamma[8,0,4] + delta[8,0,4] + epsilon[8,0] <= 0
 v[8,0,5]: 920 beta[8,0] + gamma[8,0,5] + delta[8,0,5] + epsilon[8,0] <= 0
 v[8,0,6]: 534 beta[8,0] + gamma[8,0,6] + delta[8,0,6] + epsilon[8,0] <= 0
 v[8,0,7]: 425 beta[8,0] + gamma[8,0,7] + delta[8,0,7] + epsilon[8,0] <= 0
 v[8,0,8]: 382 beta[8,0] + gamma[8,0,8] + delta[8,0,8] + epsilon[8,0] <= 0
 v[8,0,9]: 563 beta[8,0] + gamma[8,0,9] + delta[8,0,9] + epsilon[8,0] <= 0
 v[8,1,0]: 587 beta[8,1] + gamma[8,1,0] + delta[8,1,0] + epsilon[8,1] <= 0
 v[8,1,1]: gamma[8,1,1] + delta[8,1,1] + epsilon[8,1] <= 0
 v[8,1,2]: 199 beta[8,1] + gamma[8,1,2] + delta[8,1,2] + epsilon[8,1] <= 0
 v[8,1,3]: 532 beta[8,1] + gamma[8,1,3] + delta[8,1,3] + epsilon[8,1] <= 0
 v[8,1,4]: 897 beta[8,1] + gamma[8,1,4] + delta[8,1,4] + epsilon[8,1] <= 0
 v[8,1,5]: 846 beta[8,1] + gamma[8,1,5] + delta[8,1,5] + epsilon[8,1] <= 0
 v[8,1,6]: 309 beta[8,1] + gamma[8,1,6] + delta[8,1,6] + epsilon[8,1] <= 0
 v[8,1,7]: 785 beta[8,1] + gamma[8,1,7] + delta[8,1,7] + epsilon[8,1] <= 0
 v[8,1,8]: 343 beta[8,1] + gamma[8,1,8] + delta[8,1,8] + epsilon[8,1] <= 0
 v[8,1,9]: 501 beta[8,1] + gamma[8,1,9] + delta[8,1,9] + epsilon[8,1] <= 0
 v[8,2,0]: 931 beta[8,2] + gamma[8,2,0] + delta[8,2,0] + epsilon[8,2] <= 0
 v[8,2,1]: 362 beta[8,2] + gamma[8,2,1] + delta[8,2,1] + epsilon[8,2] <= 0
 v[8,2,2]: gamma[8,2,2] + delta[8,2,2] + epsilon[8,2] <= 0
 v[8,2,3]: 953 beta[8,2] + gamma[8,2,3] + delta[8,2,3] + epsilon[8,2] <= 0
 v[8,2,4]: 528 beta[8,2] + gamma[8,2,4] + delta[8,2,4] + epsilon[8,2] <= 0
 v[8,2,5]: 554 beta[8,2] + gamma[8,2,5] + delta[8,2,5] + epsilon[8,2] <= 0
 v[8,2,6]: 194 beta[8,2] + gamma[8,2,6] + delta[8,2,6] + epsilon[8,2] <= 0
 v[8,2,7]: 728 beta[8,2] + gamma[8,2,7] + delta[8,2,7] + epsilon[8,2] <= 0
 v[8,2,8]: 787 beta[8,2] + gamma[8,2,8] + delta[8,2,8] + epsilon[8,2] <= 0
 v[8,2,9]: 1041 beta[8,2] + gamma[8,2,9] + delta[8,2,9] + epsilon[8,2]
   <= 0
 v[8,3,0]: 317 beta[8,3] + gamma[8,3,0] + delta[8,3,0] + epsilon[8,3] <= 0
 v[8,3,1]: 124 beta[8,3] + gamma[8,3,1] + delta[8,3,1] + epsilon[8,3] <= 0
 v[8,3,2]: 965 beta[8,3] + gamma[8,3,2] + delta[8,3,2] + epsilon[8,3] <= 0
 v[8,3,3]: gamma[8,3,3] + delta[8,3,3] + epsilon[8,3] <= 0
 v[8,3,4]: 912 beta[8,3] + gamma[8,3,4] + delta[8,3,4] + epsilon[8,3] <= 0
 v[8,3,5]: 757 beta[8,3] + gamma[8,3,5] + delta[8,3,5] + epsilon[8,3] <= 0
 v[8,3,6]: 539 beta[8,3] + gamma[8,3,6] + delta[8,3,6] + epsilon[8,3] <= 0
 v[8,3,7]: 814 beta[8,3] + gamma[8,3,7] + delta[8,3,7] + epsilon[8,3] <= 0
 v[8,3,8]: 949 beta[8,3] + gamma[8,3,8] + delta[8,3,8] + epsilon[8,3] <= 0
 v[8,3,9]: 365 beta[8,3] + gamma[8,3,9] + delta[8,3,9] + epsilon[8,3] <= 0
 v[8,4,0]: 591 beta[8,4] + gamma[8,4,0] + delta[8,4,0] + epsilon[8,4] <= 0
 v[8,4,1]: 437 beta[8,4] + gamma[8,4,1] + delta[8,4,1] + epsilon[8,4] <= 0
 v[8,4,2]: 849 beta[8,4] + gamma[8,4,2] + delta[8,4,2] + epsilon[8,4] <= 0
 v[8,4,3]: 452 beta[8,4] + gamma[8,4,3] + delta[8,4,3] + epsilon[8,4] <= 0
 v[8,4,4]: gamma[8,4,4] + delta[8,4,4] + epsilon[8,4] <= 0
 v[8,4,5]: 867 beta[8,4] + gamma[8,4,5] + delta[8,4,5] + epsilon[8,4] <= 0
 v[8,4,6]: 735 beta[8,4] + gamma[8,4,6] + delta[8,4,6] + epsilon[8,4] <= 0
 v[8,4,7]: 189 beta[8,4] + gamma[8,4,7] + delta[8,4,7] + epsilon[8,4] <= 0
 v[8,4,8]: 946 beta[8,4] + gamma[8,4,8] + delta[8,4,8] + epsilon[8,4] <= 0
 v[8,4,9]: 694 beta[8,4] + gamma[8,4,9] + delta[8,4,9] + epsilon[8,4] <= 0
 v[8,5,0]: 1079 beta[8,5] + gamma[8,5,0] + delta[8,5,0] + epsilon[8,5]
   <= 0
 v[8,5,1]: 923 beta[8,5] + gamma[8,5,1] + delta[8,5,1] + epsilon[8,5] <= 0
 v[8,5,2]: 125 beta[8,5] + gamma[8,5,2] + delta[8,5,2] + epsilon[8,5] <= 0
 v[8,5,3]: 228 beta[8,5] + gamma[8,5,3] + delta[8,5,3] + epsilon[8,5] <= 0
 v[8,5,4]: 173 beta[8,5] + gamma[8,5,4] + delta[8,5,4] + epsilon[8,5] <= 0
 v[8,5,5]: gamma[8,5,5] + delta[8,5,5] + epsilon[8,5] <= 0
 v[8,5,6]: 394 beta[8,5] + gamma[8,5,6] + delta[8,5,6] + epsilon[8,5] <= 0
 v[8,5,7]: 581 beta[8,5] + gamma[8,5,7] + delta[8,5,7] + epsilon[8,5] <= 0
 v[8,5,8]: 388 beta[8,5] + gamma[8,5,8] + delta[8,5,8] + epsilon[8,5] <= 0
 v[8,5,9]: 195 beta[8,5] + gamma[8,5,9] + delta[8,5,9] + epsilon[8,5] <= 0
 v[8,6,0]: 150 beta[8,6] + gamma[8,6,0] + delta[8,6,0] + epsilon[8,6] <= 0
 v[8,6,1]: 565 beta[8,6] + gamma[8,6,1] + delta[8,6,1] + epsilon[8,6] <= 0
 v[8,6,2]: 215 beta[8,6] + gamma[8,6,2] + delta[8,6,2] + epsilon[8,6] <= 0
 v[8,6,3]: 565 beta[8,6] + gamma[8,6,3] + delta[8,6,3] + epsilon[8,6] <= 0
 v[8,6,4]: 485 beta[8,6] + gamma[8,6,4] + delta[8,6,4] + epsilon[8,6] <= 0
 v[8,6,5]: 104 beta[8,6] + gamma[8,6,5] + delta[8,6,5] + epsilon[8,6] <= 0
 v[8,6,6]: gamma[8,6,6] + delta[8,6,6] + epsilon[8,6] <= 0
 v[8,6,7]: 420 beta[8,6] + gamma[8,6,7] + delta[8,6,7] + epsilon[8,6] <= 0
 v[8,6,8]: 836 beta[8,6] + gamma[8,6,8] + delta[8,6,8] + epsilon[8,6] <= 0
 v[8,6,9]: 521 beta[8,6] + gamma[8,6,9] + delta[8,6,9] + epsilon[8,6] <= 0
 v[8,7,0]: 402 beta[8,7] + gamma[8,7,0] + delta[8,7,0] + epsilon[8,7] <= 0
 v[8,7,1]: 104 beta[8,7] + gamma[8,7,1] + delta[8,7,1] + epsilon[8,7] <= 0
 v[8,7,2]: 250 beta[8,7] + gamma[8,7,2] + delta[8,7,2] + epsilon[8,7] <= 0
 v[8,7,3]: 696 beta[8,7] + gamma[8,7,3] + delta[8,7,3] + epsilon[8,7] <= 0
 v[8,7,4]: 937 beta[8,7] + gamma[8,7,4] + delta[8,7,4] + epsilon[8,7] <= 0
 v[8,7,5]: 589 beta[8,7] + gamma[8,7,5] + delta[8,7,5] + epsilon[8,7] <= 0
 v[8,7,6]: 302 beta[8,7] + gamma[8,7,6] + delta[8,7,6] + epsilon[8,7] <= 0
 v[8,7,7]: gamma[8,7,7] + delta[8,7,7] + epsilon[8,7] <= 0
 v[8,7,8]: 678 beta[8,7] + gamma[8,7,8] + delta[8,7,8] + epsilon[8,7] <= 0
 v[8,7,9]: 995 beta[8,7] + gamma[8,7,9] + delta[8,7,9] + epsilon[8,7] <= 0
 v[8,8,0]: 818 beta[8,8] + gamma[8,8,0] + delta[8,8,0] + epsilon[8,8] <= 0
 v[8,8,1]: 309 beta[8,8] + gamma[8,8,1] + delta[8,8,1] + epsilon[8,8] <= 0
 v[8,8,2]: 232 beta[8,8] + gamma[8,8,2] + delta[8,8,2] + epsilon[8,8] <= 0
 v[8,8,3]: 783 beta[8,8] + gamma[8,8,3] + delta[8,8,3] + epsilon[8,8] <= 0
 v[8,8,4]: 675 beta[8,8] + gamma[8,8,4] + delta[8,8,4] + epsilon[8,8] <= 0
 v[8,8,5]: 973 beta[8,8] + gamma[8,8,5] + delta[8,8,5] + epsilon[8,8] <= 0
 v[8,8,6]: 962 beta[8,8] + gamma[8,8,6] + delta[8,8,6] + epsilon[8,8] <= 0
 v[8,8,7]: 802 beta[8,8] + gamma[8,8,7] + delta[8,8,7] + epsilon[8,8] <= 0
 v[8,8,8]: gamma[8,8,8] + delta[8,8,8] + epsilon[8,8] <= 0
 v[8,8,9]: 968 beta[8,8] + gamma[8,8,9] + delta[8,8,9] + epsilon[8,8] <= 0
 v[8,9,0]: 115 beta[8,9] + gamma[8,9,0] + delta[8,9,0] + epsilon[8,9] <= 0
 v[8,9,1]: 791 beta[8,9] + gamma[8,9,1] + delta[8,9,1] + epsilon[8,9] <= 0
 v[8,9,2]: 257 beta[8,9] + gamma[8,9,2] + delta[8,9,2] + epsilon[8,9] <= 0
 v[8,9,3]: 160 beta[8,9] + gamma[8,9,3] + delta[8,9,3] + epsilon[8,9] <= 0
 v[8,9,4]: 1046 beta[8,9] + gamma[8,9,4] + delta[8,9,4] + epsilon[8,9]
   <= 0
 v[8,9,5]: 296 beta[8,9] + gamma[8,9,5] + delta[8,9,5] + epsilon[8,9] <= 0
 v[8,9,6]: 1097 beta[8,9] + gamma[8,9,6] + delta[8,9,6] + epsilon[8,9]
   <= 0
 v[8,9,7]: 188 beta[8,9] + gamma[8,9,7] + delta[8,9,7] + epsilon[8,9] <= 0
 v[8,9,8]: 632 beta[8,9] + gamma[8,9,8] + delta[8,9,8] + epsilon[8,9] <= 0
 v[8,9,9]: gamma[8,9,9] + delta[8,9,9] + epsilon[8,9] <= 0
 v[9,0,0]: gamma[9,0,0] + delta[9,0,0] + epsilon[9,0] <= 0
 v[9,0,1]: 487 beta[9,0] + gamma[9,0,1] + delta[9,0,1] + epsilon[9,0] <= 0
 v[9,0,2]: 542 beta[9,0] + gamma[9,0,2] + delta[9,0,2] + epsilon[9,0] <= 0
 v[9,0,3]: 1028 beta[9,0] + gamma[9,0,3] + delta[9,0,3] + epsilon[9,0]
   <= 0
 v[9,0,4]: 160 beta[9,0] + gamma[9,0,4] + delta[9,0,4] + epsilon[9,0] <= 0
 v[9,0,5]: 1003 beta[9,0] + gamma[9,0,5] + delta[9,0,5] + epsilon[9,0]
   <= 0
 v[9,0,6]: 1058 beta[9,0] + gamma[9,0,6] + delta[9,0,6] + epsilon[9,0]
   <= 0
 v[9,0,7]: 882 beta[9,0] + gamma[9,0,7] + delta[9,0,7] + epsilon[9,0] <= 0
 v[9,0,8]: 604 beta[9,0] + gamma[9,0,8] + delta[9,0,8] + epsilon[9,0] <= 0
 v[9,0,9]: 484 beta[9,0] + gamma[9,0,9] + delta[9,0,9] + epsilon[9,0] <= 0
 v[9,1,0]: 169 beta[9,1] + gamma[9,1,0] + delta[9,1,0] + epsilon[9,1] <= 0
 v[9,1,1]: gamma[9,1,1] + delta[9,1,1] + epsilon[9,1] <= 0
 v[9,1,2]: 704 beta[9,1] + gamma[9,1,2] + delta[9,1,2] + epsilon[9,1] <= 0
 v[9,1,3]: 433 beta[9,1] + gamma[9,1,3] + delta[9,1,3] + epsilon[9,1] <= 0
 v[9,1,4]: 857 beta[9,1] + gamma[9,1,4] + delta[9,1,4] + epsilon[9,1] <= 0
 v[9,1,5]: 575 beta[9,1] + gamma[9,1,5] + delta[9,1,5] + epsilon[9,1] <= 0
 v[9,1,6]: 1004 beta[9,1] + gamma[9,1,6] + delta[9,1,6] + epsilon[9,1]
   <= 0
 v[9,1,7]: 827 beta[9,1] + gamma[9,1,7] + delta[9,1,7] + epsilon[9,1] <= 0
 v[9,1,8]: 128 beta[9,1] + gamma[9,1,8] + delta[9,1,8] + epsilon[9,1] <= 0
 v[9,1,9]: 993 beta[9,1] + gamma[9,1,9] + delta[9,1,9] + epsilon[9,1] <= 0
 v[9,2,0]: 126 beta[9,2] + gamma[9,2,0] + delta[9,2,0] + epsilon[9,2] <= 0
 v[9,2,1]: 317 beta[9,2] + gamma[9,2,1] + delta[9,2,1] + epsilon[9,2] <= 0
 v[9,2,2]: gamma[9,2,2] + delta[9,2,2] + epsilon[9,2] <= 0
 v[9,2,3]: 1061 beta[9,2] + gamma[9,2,3] + delta[9,2,3] + epsilon[9,2]
   <= 0
 v[9,2,4]: 324 beta[9,2] + gamma[9,2,4] + delta[9,2,4] + epsilon[9,2] <= 0
 v[9,2,5]: 437 beta[9,2] + gamma[9,2,5] + delta[9,2,5] + epsilon[9,2] <= 0
 v[9,2,6]: 812 beta[9,2] + gamma[9,2,6] + delta[9,2,6] + epsilon[9,2] <= 0
 v[9,2,7]: 514 beta[9,2] + gamma[9,2,7] + delta[9,2,7] + epsilon[9,2] <= 0
 v[9,2,8]: 392 beta[9,2] + gamma[9,2,8] + delta[9,2,8] + epsilon[9,2] <= 0
 v[9,2,9]: 509 beta[9,2] + gamma[9,2,9] + delta[9,2,9] + epsilon[9,2] <= 0
 v[9,3,0]: 1085 beta[9,3] + gamma[9,3,0] + delta[9,3,0] + epsilon[9,3]
   <= 0
 v[9,3,1]: 660 beta[9,3] + gamma[9,3,1] + delta[9,3,1] + epsilon[9,3] <= 0
 v[9,3,2]: 144 beta[9,3] + gamma[9,3,2] + delta[9,3,2] + epsilon[9,3] <= 0
 v[9,3,3]: gamma[9,3,3] + delta[9,3,3] + epsilon[9,3] <= 0
 v[9,3,4]: 994 beta[9,3] + gamma[9,3,4] + delta[9,3,4] + epsilon[9,3] <= 0
 v[9,3,5]: 138 beta[9,3] + gamma[9,3,5] + delta[9,3,5] + epsilon[9,3] <= 0
 v[9,3,6]: 1068 beta[9,3] + gamma[9,3,6] + delta[9,3,6] + epsilon[9,3]
   <= 0
 v[9,3,7]: 890 beta[9,3] + gamma[9,3,7] + delta[9,3,7] + epsilon[9,3] <= 0
 v[9,3,8]: 864 beta[9,3] + gamma[9,3,8] + delta[9,3,8] + epsilon[9,3] <= 0
 v[9,3,9]: 348 beta[9,3] + gamma[9,3,9] + delta[9,3,9] + epsilon[9,3] <= 0
 v[9,4,0]: 570 beta[9,4] + gamma[9,4,0] + delta[9,4,0] + epsilon[9,4] <= 0
 v[9,4,1]: 1001 beta[9,4] + gamma[9,4,1] + delta[9,4,1] + epsilon[9,4]
   <= 0
 v[9,4,2]: 211 beta[9,4] + gamma[9,4,2] + delta[9,4,2] + epsilon[9,4] <= 0
 v[9,4,3]: 776 beta[9,4] + gamma[9,4,3] + delta[9,4,3] + epsilon[9,4] <= 0
 v[9,4,4]: gamma[9,4,4] + delta[9,4,4] + epsilon[9,4] <= 0
 v[9,4,5]: 449 beta[9,4] + gamma[9,4,5] + delta[9,4,5] + epsilon[9,4] <= 0
 v[9,4,6]: 991 beta[9,4] + gamma[9,4,6] + delta[9,4,6] + epsilon[9,4] <= 0
 v[9,4,7]: 215 beta[9,4] + gamma[9,4,7] + delta[9,4,7] + epsilon[9,4] <= 0
 v[9,4,8]: 135 beta[9,4] + gamma[9,4,8] + delta[9,4,8] + epsilon[9,4] <= 0
 v[9,4,9]: 367 beta[9,4] + gamma[9,4,9] + delta[9,4,9] + epsilon[9,4] <= 0
 v[9,5,0]: 523 beta[9,5] + gamma[9,5,0] + delta[9,5,0] + epsilon[9,5] <= 0
 v[9,5,1]: 105 beta[9,5] + gamma[9,5,1] + delta[9,5,1] + epsilon[9,5] <= 0
 v[9,5,2]: 174 beta[9,5] + gamma[9,5,2] + delta[9,5,2] + epsilon[9,5] <= 0
 v[9,5,3]: 500 beta[9,5] + gamma[9,5,3] + delta[9,5,3] + epsilon[9,5] <= 0
 v[9,5,4]: 714 beta[9,5] + gamma[9,5,4] + delta[9,5,4] + epsilon[9,5] <= 0
 v[9,5,5]: gamma[9,5,5] + delta[9,5,5] + epsilon[9,5] <= 0
 v[9,5,6]: 621 beta[9,5] + gamma[9,5,6] + delta[9,5,6] + epsilon[9,5] <= 0
 v[9,5,7]: 647 beta[9,5] + gamma[9,5,7] + delta[9,5,7] + epsilon[9,5] <= 0
 v[9,5,8]: 464 beta[9,5] + gamma[9,5,8] + delta[9,5,8] + epsilon[9,5] <= 0
 v[9,5,9]: 312 beta[9,5] + gamma[9,5,9] + delta[9,5,9] + epsilon[9,5] <= 0
 v[9,6,0]: 607 beta[9,6] + gamma[9,6,0] + delta[9,6,0] + epsilon[9,6] <= 0
 v[9,6,1]: 1041 beta[9,6] + gamma[9,6,1] + delta[9,6,1] + epsilon[9,6]
   <= 0
 v[9,6,2]: 201 beta[9,6] + gamma[9,6,2] + delta[9,6,2] + epsilon[9,6] <= 0
 v[9,6,3]: 1091 beta[9,6] + gamma[9,6,3] + delta[9,6,3] + epsilon[9,6]
   <= 0
 v[9,6,4]: 851 beta[9,6] + gamma[9,6,4] + delta[9,6,4] + epsilon[9,6] <= 0
 v[9,6,5]: 1076 beta[9,6] + gamma[9,6,5] + delta[9,6,5] + epsilon[9,6]
   <= 0
 v[9,6,6]: gamma[9,6,6] + delta[9,6,6] + epsilon[9,6] <= 0
 v[9,6,7]: 316 beta[9,6] + gamma[9,6,7] + delta[9,6,7] + epsilon[9,6] <= 0
 v[9,6,8]: 1011 beta[9,6] + gamma[9,6,8] + delta[9,6,8] + epsilon[9,6]
   <= 0
 v[9,6,9]: 157 beta[9,6] + gamma[9,6,9] + delta[9,6,9] + epsilon[9,6] <= 0
 v[9,7,0]: 518 beta[9,7] + gamma[9,7,0] + delta[9,7,0] + epsilon[9,7] <= 0
 v[9,7,1]: 121 beta[9,7] + gamma[9,7,1] + delta[9,7,1] + epsilon[9,7] <= 0
 v[9,7,2]: 989 beta[9,7] + gamma[9,7,2] + delta[9,7,2] + epsilon[9,7] <= 0
 v[9,7,3]: 391 beta[9,7] + gamma[9,7,3] + delta[9,7,3] + epsilon[9,7] <= 0
 v[9,7,4]: 127 beta[9,7] + gamma[9,7,4] + delta[9,7,4] + epsilon[9,7] <= 0
 v[9,7,5]: 421 beta[9,7] + gamma[9,7,5] + delta[9,7,5] + epsilon[9,7] <= 0
 v[9,7,6]: 378 beta[9,7] + gamma[9,7,6] + delta[9,7,6] + epsilon[9,7] <= 0
 v[9,7,7]: gamma[9,7,7] + delta[9,7,7] + epsilon[9,7] <= 0
 v[9,7,8]: 508 beta[9,7] + gamma[9,7,8] + delta[9,7,8] + epsilon[9,7] <= 0
 v[9,7,9]: 948 beta[9,7] + gamma[9,7,9] + delta[9,7,9] + epsilon[9,7] <= 0
 v[9,8,0]: 547 beta[9,8] + gamma[9,8,0] + delta[9,8,0] + epsilon[9,8] <= 0
 v[9,8,1]: 527 beta[9,8] + gamma[9,8,1] + delta[9,8,1] + epsilon[9,8] <= 0
 v[9,8,2]: 297 beta[9,8] + gamma[9,8,2] + delta[9,8,2] + epsilon[9,8] <= 0
 v[9,8,3]: 1091 beta[9,8] + gamma[9,8,3] + delta[9,8,3] + epsilon[9,8]
   <= 0
 v[9,8,4]: 908 beta[9,8] + gamma[9,8,4] + delta[9,8,4] + epsilon[9,8] <= 0
 v[9,8,5]: 914 beta[9,8] + gamma[9,8,5] + delta[9,8,5] + epsilon[9,8] <= 0
 v[9,8,6]: 863 beta[9,8] + gamma[9,8,6] + delta[9,8,6] + epsilon[9,8] <= 0
 v[9,8,7]: 785 beta[9,8] + gamma[9,8,7] + delta[9,8,7] + epsilon[9,8] <= 0
 v[9,8,8]: gamma[9,8,8] + delta[9,8,8] + epsilon[9,8] <= 0
 v[9,8,9]: 1004 beta[9,8] + gamma[9,8,9] + delta[9,8,9] + epsilon[9,8]
   <= 0
 v[9,9,0]: 455 beta[9,9] + gamma[9,9,0] + delta[9,9,0] + epsilon[9,9] <= 0
 v[9,9,1]: 674 beta[9,9] + gamma[9,9,1] + delta[9,9,1] + epsilon[9,9] <= 0
 v[9,9,2]: 782 beta[9,9] + gamma[9,9,2] + delta[9,9,2] + epsilon[9,9] <= 0
 v[9,9,3]: 300 beta[9,9] + gamma[9,9,3] + delta[9,9,3] + epsilon[9,9] <= 0
 v[9,9,4]: 123 beta[9,9] + gamma[9,9,4] + delta[9,9,4] + epsilon[9,9] <= 0
 v[9,9,5]: 632 beta[9,9] + gamma[9,9,5] + delta[9,9,5] + epsilon[9,9] <= 0
 v[9,9,6]: 266 beta[9,9] + gamma[9,9,6] + delta[9,9,6] + epsilon[9,9] <= 0
 v[9,9,7]: 970 beta[9,9] + gamma[9,9,7] + delta[9,9,7] + epsilon[9,9] <= 0
 v[9,9,8]: 685 beta[9,9] + gamma[9,9,8] + delta[9,9,8] + epsilon[9,9] <= 0
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
 epsilon[0,0] free
 epsilon[0,1] free
 epsilon[0,2] free
 epsilon[0,3] free
 epsilon[0,4] free
 epsilon[0,5] free
 epsilon[0,6] free
 epsilon[0,7] free
 epsilon[0,8] free
 epsilon[0,9] free
 epsilon[1,0] free
 epsilon[1,1] free
 epsilon[1,2] free
 epsilon[1,3] free
 epsilon[1,4] free
 epsilon[1,5] free
 epsilon[1,6] free
 epsilon[1,7] free
 epsilon[1,8] free
 epsilon[1,9] free
 epsilon[2,0] free
 epsilon[2,1] free
 epsilon[2,2] free
 epsilon[2,3] free
 epsilon[2,4] free
 epsilon[2,5] free
 epsilon[2,6] free
 epsilon[2,7] free
 epsilon[2,8] free
 epsilon[2,9] free
 epsilon[3,0] free
 epsilon[3,1] free
 epsilon[3,2] free
 epsilon[3,3] free
 epsilon[3,4] free
 epsilon[3,5] free
 epsilon[3,6] free
 epsilon[3,7] free
 epsilon[3,8] free
 epsilon[3,9] free
 epsilon[4,0] free
 epsilon[4,1] free
 epsilon[4,2] free
 epsilon[4,3] free
 epsilon[4,4] free
 epsilon[4,5] free
 epsilon[4,6] free
 epsilon[4,7] free
 epsilon[4,8] free
 epsilon[4,9] free
 epsilon[5,0] free
 epsilon[5,1] free
 epsilon[5,2] free
 epsilon[5,3] free
 epsilon[5,4] free
 epsilon[5,5] free
 epsilon[5,6] free
 epsilon[5,7] free
 epsilon[5,8] free
 epsilon[5,9] free
 epsilon[6,0] free
 epsilon[6,1] free
 epsilon[6,2] free
 epsilon[6,3] free
 epsilon[6,4] free
 epsilon[6,5] free
 epsilon[6,6] free
 epsilon[6,7] free
 epsilon[6,8] free
 epsilon[6,9] free
 epsilon[7,0] free
 epsilon[7,1] free
 epsilon[7,2] free
 epsilon[7,3] free
 epsilon[7,4] free
 epsilon[7,5] free
 epsilon[7,6] free
 epsilon[7,7] free
 epsilon[7,8] free
 epsilon[7,9] free
 epsilon[8,0] free
 epsilon[8,1] free
 epsilon[8,2] free
 epsilon[8,3] free
 epsilon[8,4] free
 epsilon[8,5] free
 epsilon[8,6] free
 epsilon[8,7] free
 epsilon[8,8] free
 epsilon[8,9] free
 epsilon[9,0] free
 epsilon[9,1] free
 epsilon[9,2] free
 epsilon[9,3] free
 epsilon[9,4] free
 epsilon[9,5] free
 epsilon[9,6] free
 epsilon[9,7] free
 epsilon[9,8] free
 epsilon[9,9] free
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
 nu[0] free
 nu[1] free
 nu[2] free
 nu[3] free
 nu[4] free
 nu[5] free
 nu[6] free
 nu[7] free
 nu[8] free
 nu[9] free
End
