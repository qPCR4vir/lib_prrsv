import string

def Compute_GC (seq):
    GC = 0
    total =0
    for i in range(0, len(seq)):
        if string.upper(seq[i]) in ['A','T','G','C']:
            total = total +1
            if string.upper(seq[i]) in ['G', 'C']:
                GC = GC+1

    return [total, GC]

















