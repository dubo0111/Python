import LIP0 as one
import LIP as two


ni=6
nk=ni
rnd = 38
a1=0.2


L0,L0S,cd,y,x,u,v,L03,coordinate = two.two_stage(ni,nk,rnd,a1)
L1,L1S,_,y1,x1,u1,v1,L13 = one.one_stage(ni,nk,rnd,a1)

#for n in range(10000):
#    L0,L0S = two.two_stage(ni,nk,rnd,a1)
#    L1,L1S = one.one_stage(ni,nk,rnd,a1)
#    if L0 > L1 and L0S < L1S and abs(L0-L0S)>0.01:
#        print('$$$$$$$$$$  ',rnd)
#        break
#    rnd += 1

## reserve
#    rnd:38
