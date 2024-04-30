from dataclasses import dataclass
from typing import List

@dataclass
class SequenceInfo:
    sequence: List[int]
    n_traj: int # number of trajectories sampled so far
    v: float # V_s^n in paper  (eqn 14) (performance measure)
    queue_probs: List[float] # stores "probability" that aircraft will immediately go into service at end of remaining travel time - upsilon in paper (page 18)
    # Rob thinks this upsilon is redundant as Rob set lambda to zero in paper which means that if one random trajectory where aircraft enters
    # service immediately, this triggers release from pool.
    w: float # W_s^n in paper (eqn 15)
    age: int # number of times sequence has been passed to Repopulate
