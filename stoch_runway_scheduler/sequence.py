from __future__ import annotations
from dataclasses import dataclass, field
import math
from typing import List

@dataclass
class SequenceInfo:
    sequence: List[int]
    n_traj: int = field(init=False) # number of trajectories sampled so far
    v: float = field(init=False) # V_s^n in paper  (eqn 14) (performance measure)
    queue_probs: List[float] = field(init=False) # stores "probability" that aircraft will immediately go into service at end of remaining travel time - upsilon in paper (page 18)
    w: float = field(init=False) # W_s^n in paper (eqn 15)
    age: int = field(init=False) # number of times sequence has been passed to Repopulate

    def __post_init__(self):
        self.n_traj = 0
        self.v = 0
        self.w = 0
        self.age = 0
        self.queue_probs = len(self.sequence)*[0.0]

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

    def mean_cost(self):
        "Mean cost of sampled trajectories for sequence"
        return self.v

    def var_cost(self):
        "Variance of cost of sampled trajectories for sequence"
        n = self.n_traj
        mn = self.mean_cost()
        return (self.w-(n*mn**2)) / (n - 1)

    def reset(self):
        "Re-initialise all data on sequence except for age."
        self.n_traj = 0
        self.v = 0
        self.w = 0
        for (i, qp) in enumerate(self.queue_probs):
            self.queue_probs[i] = 0

    @staticmethod
    def rank_and_select(info_list: List[SequenceInfo]) -> List[int]:
        """Uses rank and select procedure to identify sequences which should be removed.

        Arguments
        ---------
        info_list: list of SequenceInfo objects to analyse

        Returns
        -------
        to_remove: list of indices of sequences in input list which should be removed
        """
        
        t_val = 1.96 # 97.5th percentile of normal dist
        keep = len(info_list) * [True]
        # Ranking and selection
        for (i, info) in enumerate(info_list):
            if keep[i]:
                mn1 = info.mean_cost()
                n1 = info.n_traj
                var1 = info.var_cost()

                # This is in Section 3.1 (2B of algorithm) of the paper (Equations 14-17)
                for (j, infoj) in enumerate(info_list):
                    if keep[j]:
                        mn2 = infoj.v
                        n2 = infoj.n_traj
                        var2 = infoj.var_cost()
                        w_val = math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1 > mn2 + w_val:
                            keep[i] = False # This marks flights for removal
                            break

                        elif mn2 > mn1 + w_val:
                            keep[j] = False
        
        return [i for (i, k) in enumerate(keep) if not k]