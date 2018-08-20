\ Model p-center
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  L
Subject To
 epigraph[0]: 878 x[0,1] + 141 x[0,2] + 699 x[0,3] - L <= 0
 epigraph[1]: 399 x[1,0] + 1097 x[1,2] + 726 x[1,3] - L <= 0
 epigraph[2]: 349 x[2,0] + 102 x[2,1] + 135 x[2,3] - L <= 0
 epigraph[3]: 446 x[3,0] + 663 x[3,1] + 425 x[3,2] - L <= 0
 p: y[0] + y[1] + y[2] + y[3] = 1
 x<y[0,0]: x[0,0] - y[0] <= 0
 x<y[0,1]: x[0,1] - y[1] <= 0
 x<y[0,2]: x[0,2] - y[2] <= 0
 x<y[0,3]: x[0,3] - y[3] <= 0
 x<y[1,0]: x[1,0] - y[0] <= 0
 x<y[1,1]: x[1,1] - y[1] <= 0
 x<y[1,2]: x[1,2] - y[2] <= 0
 x<y[1,3]: x[1,3] - y[3] <= 0
 x<y[2,0]: x[2,0] - y[0] <= 0
 x<y[2,1]: x[2,1] - y[1] <= 0
 x<y[2,2]: x[2,2] - y[2] <= 0
 x<y[2,3]: x[2,3] - y[3] <= 0
 x<y[3,0]: x[3,0] - y[0] <= 0
 x<y[3,1]: x[3,1] - y[1] <= 0
 x<y[3,2]: x[3,2] - y[2] <= 0
 x<y[3,3]: x[3,3] - y[3] <= 0
 sumx[0]: x[0,0] + x[0,1] + x[0,2] + x[0,3] = 1
 sumx[1]: x[1,0] + x[1,1] + x[1,2] + x[1,3] = 1
 sumx[2]: x[2,0] + x[2,1] + x[2,2] + x[2,3] = 1
 sumx[3]: x[3,0] + x[3,1] + x[3,2] + x[3,3] = 1
Bounds
Binaries
 y[0] y[1] y[2] y[3]
End
