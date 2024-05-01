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
    def add_observation(self, cost: float, xi_list: List[bool]):
        """Uses latest (possibly sampled) observation update state of information on sequence
            
            Arguments
            ----------
            cost: overall sampled cost of sequence
            xi_list: list indicating for each flight in sequence whether or not it went straight into service
        """
        self.n_traj += 1
        gam=1/self.n_traj
        self.v = (1-gam)*self.v + gam*cost
        self.w += cost**2
        for (i, straight_into_service) in enumerate(xi_list):
            self.queue_probs[i] = (1 - gam) * self.queue_probs[i] + gam * straight_into_service

