\ Model p-center
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0 y[0] + 0 y[1] + 0 y[2] + L
Subject To
 Capacity[0]: x[0,0] + x[0,1] + x[0,2] - L <= 0
 Capacity[1]: x[1,0] + x[1,1] + x[1,2] - L <= 0
 Capacity[2]: x[2,0] + x[2,1] + x[2,2] - L <= 0
Bounds
Binaries
 x[0,0] x[0,1] x[0,2] x[1,0] x[1,1] x[1,2] x[2,0] x[2,1] x[2,2] y[0] y[1]
 y[2]
End
