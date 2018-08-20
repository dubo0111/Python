\ Model master model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.5 L
Subject To
 epigraph[0]: 517 x[0,1] + 100 x[0,2] - L <= 0
 epigraph[1]: 247 x[1,0] + 286 x[1,2] - L <= 0
 epigraph[2]: 497 x[2,0] + 519 x[2,1] - L <= 0
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
