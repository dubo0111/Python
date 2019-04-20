import numpy as np

def primal_integral(conv1,conv12,conv13,conv14,timeline = 0):
    #
    def integral_compute(t,time,gap):
        pi = 0
        # time_count = [x for x in time if x <= t]
        for n in range(t-1):
            time = np.array(time)
            a = np.argmax(time>n)
            b = np.argmax(time>n+1)
            if a == 0:
                a = 1
            if b == 0:
                b = 1
            pi += (gap[a-1]+gap[b-1])/2
        return pi
    # get optimality
    best = [conv1[0][-1],conv12[0][-1],conv13[0][-1],conv14[0][-1]]
    index_min = min(range(len(best)), key=best.__getitem__)
    opt_val = best[index_min]
    # get min time
    time = [conv1[2][-1],conv12[2][-1],conv13[2][-1],conv14[2][-1]]
    index_mint = min(range(len(time)), key=time.__getitem__)
    min_time = time[index_mint]
    # process conv
    gap1 = [(conv1[0][n]-opt_val)/conv1[0][n] for n in range(len(conv1[0]))]
    gap2 = [(conv12[0][n]-opt_val)/conv12[0][n] for n in range(len(conv12[0]))]
    gap3 = [(conv13[0][n]-opt_val)/conv13[0][n] for n in range(len(conv13[0]))]
    gap4 = [(conv14[0][n]-opt_val)/conv14[0][n] for n in range(len(conv14[0]))]
    if timeline == 0:
        timeline = [min_time]
    pi_final = [[] for n in range(4)]
    # integral
    for t in timeline:
        pi_final[0].append(integral_compute(t,conv1[2],gap1))
        pi_final[1].append(integral_compute(t,conv12[2],gap2))
        pi_final[2].append(integral_compute(t,conv13[2],gap3))
        pi_final[3].append(integral_compute(t,conv14[2],gap4))
    return pi_final,opt_val
