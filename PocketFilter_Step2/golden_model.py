"""Adapted from TimD1/SquiggleFilter; see README for semantic changes."""
import numpy as np
try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        return lambda fn: fn

def discrete_normalize(seq, bits=8, minval=-4, maxval=4):
    ''' Approximate normalization which converts signal to integer of desired precision. '''
    mean = int(np.mean(seq))
    mean_avg_dev = int(np.mean(np.abs(seq - mean)))
    if mean_avg_dev == 0:
        raise ValueError("Zero integer mean absolute deviation")
    norm_seq = (seq - mean) / mean_avg_dev
    norm_seq[norm_seq < minval] = minval # threshold
    norm_seq[norm_seq > maxval] = maxval 
    norm_seq = ( (norm_seq - minval) * (2**(bits)/(maxval-minval)) ).astype(int)
    return norm_seq


@njit()
def _sdtw(seq, ref, prefix_lengths):
    ''' Returns minimum alignment score for subsequence DTW. '''
    
    # initialize column vectors
    prev_consec = np.zeros((len(seq)))
    curr_consec = np.zeros((len(seq)))
    prev_cost = np.zeros((len(seq)))
    curr_cost = np.zeros((len(seq)))
    min_cost = np.zeros((len(seq)))
    prev_cost[0] = abs(seq[0]-ref[0])
    min_cost[0] = abs(seq[0]-ref[0])
    for i in range(1, len(seq)):
        prev_cost[i] = prev_cost[i-1] + abs(seq[i]-ref[0])
        min_cost[i] = min_cost[i-1] + abs(seq[i]-ref[0])
    
    # compute entire cost matrix
    for j in range(1, len(ref)):
        bonus = 10
        for i in range(1, len(seq)):
            move = prev_cost[i-1] - prev_consec[i-1]*bonus < curr_cost[i-1]
            if move:
                curr_consec[i] = 0
                curr_cost[i] = prev_cost[i-1] - prev_consec[i-1]*bonus + abs(seq[i]-ref[j])
            else:
                curr_consec[i] = min(10, prev_consec[i] + 1)
                curr_cost[i] = curr_cost[i-1] + abs(seq[i]-ref[j])
            min_cost[i] = min(min_cost[i], curr_cost[i])
        for i in range(len(seq)):
            prev_consec[i] = curr_consec[i]
            curr_consec[i] = 0
            prev_cost[i] = curr_cost[i]
            curr_cost[i] = 0
    
    # return cost of optimal alignment
    cost_mins = np.zeros((len(prefix_lengths),))
    for i, l in enumerate(prefix_lengths):
        if l <= len(seq):
            cost_mins[i] = min_cost[l-1]
    return cost_mins

def score_preprocessed(reference, query, prefix_lengths=None):
    """Inputs must ALREADY have upstream-compatible normalization."""
    ref = np.asarray(reference, dtype=np.float64)
    seq = np.asarray(query, dtype=np.float64)
    if ref.ndim != 1 or seq.ndim != 1 or len(ref) < 1 or len(seq) < 2:
        raise ValueError("Expected nonempty 1D reference and query with >=2 samples")
    if not np.isfinite(ref).all() or not np.isfinite(seq).all():
        raise ValueError("Nonfinite signal")
    lengths = np.asarray([len(seq)] if prefix_lengths is None else prefix_lengths)
    if lengths.ndim != 1 or not np.issubdtype(lengths.dtype, np.integer):
        raise ValueError("prefix_lengths must be a 1D integer array")
    if np.any(lengths < 2) or np.any(lengths > len(seq)):
        raise ValueError("Prefix length outside query")
    return _sdtw(seq, ref, lengths.astype(np.int64))


def classify_preprocessed(reference, query, threshold):
    """Fixed query length; threshold calibrated on independent labeled reads."""
    if not np.isfinite(threshold):
        raise ValueError("Threshold must be finite")
    cost = float(score_preprocessed(reference, query)[0])
    return cost, "keep" if cost < threshold else "eject"
