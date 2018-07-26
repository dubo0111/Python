\ Model p-center
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  L
Subject To
 epigraph[0]: 2000 x[0,1] + 5000 x[0,2] - L <= 0
 epigraph[1]: 6000 x[1,0] + 8000 x[1,2] - L <= 0
 epigraph[2]: 4500 x[2,0] + 7500 x[2,1] - L <= 0
 p: y[0] + y[1] + y[2] = 1
 x<y[0,0]: x[0,0] - y[0] <= 0
 x<y[0,1]: x[0,1] - y[1] <= 0
 x<y[0,2]: x[0,2] - y[2] <= 0
 x<y[1,0]: x[1,0] - y[0] <= 0
 x<y[1,1]: x[1,1] - y[1] <= 0
 x<y[1,2]: x[1,2] - y[2] <= 0
 x<y[2,0]: x[2,0] - y[0] <= 0
 x<y[2,1]: x[2,1] - y[1] <= 0
 x<y[2,2]: x[2,2] - y[2] <= 0
 sumx[0]: x[0,0] + x[0,1] + x[0,2] = 1
 sumx[1]: x[1,0] + x[1,1] + x[1,2] = 1
 sumx[2]: x[2,0] + x[2,1] + x[2,2] = 1
Bounds
Binaries
 x[0,0] x[0,1] x[0,2] x[1,0] x[1,1] x[1,2] x[2,0] x[2,1] x[2,2] y[0] y[1]
 y[2]
End
