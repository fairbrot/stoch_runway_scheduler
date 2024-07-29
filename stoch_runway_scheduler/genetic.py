from typing import List
import logging
import numpy as np 

from .utils import FlightInfo, Cost
from .weather import StochasticWeatherProcess
from .sequence import SequenceInfo
from .separation import StochasticSeparation, landing_time
from .simulate import simulate_sequences

# Params internal to Genetic
# GA_Checksize is reset to GA_Check_Increment at beginning of each full run
# It is also reset to GA_Check_Increment when a flight is released to the queue
# and when set of sequences is repopulated. It is increased or reset after
# running rank_and_select in Genetic
# 
# GA_Check_Increment is initialised to GA_LoopSize / 10 at the beginning of every
# full run but then isn't changed.
#
# GA_LoopSize is initialized at beginning of each run depending on what
# policy is being used but then never changes
# 
# GA_Counter is incremented every time Genetic is run (inside Genetic).
# It is initialised to 0 at beginning of every full run.
# It is reset to zero after running Populate and Repopulate
# It is also incremented by 1 after Genetic_determ is run
#
# basecost is an input which isn't changed by Genetic
#
# tau is a constant input used for simulating sequences
# 
# Max_Lookahead is no longer used
# 
# cost_fn is a fixed input
# 
# S_min is fixed as is wiener_sig

# Outputs

# Ac_added is list of flights to be released to queue
# counter is number of trajectories sampled for a sequence
# from which flights are released. 0 if no flights are released.
# It doesn't seem to be used anywhere


# JF: this is the main sim heuristic
def Genetic(Ac_Info: list[FlightInfo], Ac_queue: list[int], tm, sep: StochasticSeparation, prev_class, GA_Info: list[SequenceInfo], GA_counter, basecost, weather: StochasticWeatherProcess, tau: int, cost_fn: Cost, GA_Check_Increment: int, S_min: int, wiener_sig: float):

    # Simulate costs of sequences in GA_Info
    costs, xi_lists = simulate_sequences(GA_Info, tm, Ac_Info, Ac_queue, tau, sep, weather, wiener_sig, prev_class, cost_fn)
    for info, cost, xi_list in zip(GA_Info, costs, xi_lists):
        info.add_observation(cost + basecost, xi_list)

    GA_Info.sort(key=lambda x: x.v)
    GA_counter += 1

    # Run rank and select if sufficient number of iterations have passed
    if GA_counter % GA_Check_Increment == 0:
        to_remove = SequenceInfo.rank_and_select(GA_Info)
        to_remove.sort(reverse=True)
        for i in to_remove:
            info = GA_Info.pop(i)

    # Mark some flights in best sequence for release
    Ac_added = []

    n_rel = GA_Check_Increment
    assert len(GA_Info) > 0
    if GA_counter > n_rel or (len(GA_Info) <= S_min):
        perm = GA_Info[0]
        for (j, AC) in enumerate(perm.sequence):
            if perm.queue_probs[j] <= 0: # JF Question: in paper this check is only done on first element - modifying this to be the case could allow for simplifications
                break
            Ac_added.append(AC)

    return Ac_added, GA_counter
