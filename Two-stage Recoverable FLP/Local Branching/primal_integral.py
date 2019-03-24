import numpy as np

def primal_integral(conv1,conv12,conv13,conv14,timeline = 0):
    #
    def integral_compute(t,time,gap):
        PI = 0
        time_count = [x for x in time if x <= t]
        for n in range(len(time_count)-1):
            PI += (gap[n] + gap[n+1])*(time_count[n+1]-time_count[n])/2
        return PI
    # get optimality
    best = [conv1[0][-1],conv12[0][-1],conv13[0][-1],conv14[0][-1]]
    index_min = min(range(len(best)), key=best.__getitem__)
    opt_val = best[index_min]
    # get min time
    time = [conv1[2][-1],conv12[2][-1],conv13[2][-1],conv14[2][-1]]
    index_min = min(range(len(time)), key=time.__getitem__)
    min_time = time[index_min]
    # process conv
    gap1 = [(conv1[0][n]-opt_val)/conv1[0][n] for n in range(len(conv1[0]))]
    gap2 = [(conv12[0][n]-opt_val)/conv12[0][n] for n in range(len(conv12[0]))]
    gap3 = [(conv13[0][n]-opt_val)/conv13[0][n] for n in range(len(conv13[0]))]
    gap4 = [(conv14[0][n]-opt_val)/conv14[0][n] for n in range(len(conv14[0]))]
    if timeline == 0:
        timeline = [min_time]
    PI_final = [[] for n in range(4)]
    # integral
    for t in timeline:
        PI_final[0].append(integral_compute(t,conv1[2],gap1))
        PI_final[1].append(integral_compute(t,conv12[2],gap2))
        PI_final[2].append(integral_compute(t,conv13[2],gap3))
        PI_final[3].append(integral_compute(t,conv14[2],gap4))
    return PI_final
