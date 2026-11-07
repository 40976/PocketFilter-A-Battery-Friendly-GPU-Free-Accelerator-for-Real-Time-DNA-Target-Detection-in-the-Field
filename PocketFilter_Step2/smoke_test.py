import numpy as np
from golden_model import score_preprocessed
# Independent dense implementation of the upstream recurrence (NOT textbook DTW).
def dense(ref, seq):
    n,m=len(seq),len(ref)
    c=np.zeros((n,m)); run=np.zeros((n,m))
    c[:,0]=np.cumsum(np.abs(seq-ref[0]))
    for j in range(1,m):
        for i in range(1,n):
            diag=c[i-1,j-1]-10*run[i-1,j-1]
            if diag < c[i-1,j]:
                c[i,j]=diag+abs(seq[i]-ref[j])
            else:
                c[i,j]=c[i-1,j]+abs(seq[i]-ref[j])
                run[i,j]=min(10,run[i,j-1]+1)
    return c.min(axis=1)
rng=np.random.default_rng(42)
for _ in range(30):
    ref=rng.integers(0,257,17).astype(float)
    seq=rng.integers(0,257,11).astype(float)
    lengths=np.array([2,5,11])
    np.testing.assert_array_equal(score_preprocessed(ref,seq,lengths),dense(ref,seq)[lengths-1])
print("PASS: 30 random dense-recurrence comparisons, including full-length boundary")
