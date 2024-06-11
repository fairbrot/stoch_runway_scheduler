from typing import List, Tuple, Optional
import numpy as np
from .gamma import sample_cond_gamma

def Gamma_GetServ(k: int, Time_Sep: List[List[int]], prev_class: int, cur_class: int, weather_state: int, w_rho: float, serv_time: Optional[float] = None) -> float:
    """
    Sample minimum separation time between current flight and previous one.

    This is used when flight is not already in service.

    Arguments:
    ---------
    k: Erlang shape parameter for service time
    Time_Sep: separation times between different classes
    prev_class: weight class of previous flight
    cur_class: weight class of current flight
    weather_state: code for current state of weather (0, 1 or 2)
    w_rho: multiplier for service time in case of bad weather
    serv_time: if provided, this value is appropriately scaled and used as the service time for the flight.
                If not provided, service time is sampled directly from Gamma distribution.
    """

    rate = k / (Time_Sep[prev_class][cur_class]/60)
    if weather_state == 1:
        rate *= 1/w_rho

    # Rob To Check: Can you check this function is correct? It accounts for both cases where a normalised service is provided and isn't provided
    # Note that serv_time when provided comes from a Gamma(k, 1) distribution. Is it right that we multiply the rate above by k?
    return np.random.gamma(k, 1/rate) if serv_time is None else serv_time/rate # service time


def Gamma_Conditional_GetServ(k: int, Time_Sep: List[List[int]], t_elapsed, prev_class, cur_class, weather_state, w_rho: float) -> float:
    """
    Sample minimum separation time between current flight and previous one,
    conditional on a given amount of time having elapsed from previous completion.


    Arguments:
    ---------
    k: Erlang shape parameter for service time
    Time_Sep: separation times between different classes
    t_elapsed: amount of time elapsed since previous completion
    prev_class: weight class of previous flight
    cur_class: weight class of current flight
    weather_state: code for current state of weather (0, 1 or 2)
    w_rho: multiplier for service time in case of bad weather
    """
    rate=k/(Time_Sep[prev_class][cur_class]/60)
    if weather_state==1:
        rate*=1/w_rho #=0.5

    cond_serv = sample_cond_gamma(t_elapsed, k, rate)
    return cond_serv

def landing_time(prev_comp: float, min_sep: float,
                 rel_time: float, trav_time: float) -> Tuple[float, int]:
    """
    Calculates landing time and whether flight went straight into service on reaching runway.

    Arguments
    ---------
    prev_comp: time previous aircraft finished service
    min_sep: minimum separation time between current and previous aircraft
    rel_time: time current aircraft is released into queue
    trav_time: travel time between from leaving pool to runway threshold

    Returns:
    --------
    t_out: time flight is finished being served
    straight_into_service: indicates whether flight enters service immediately on joining queue
    """
    t1 = rel_time + trav_time
    t2 = prev_comp + min_sep

    if t1 < t2:
        straight_into_service=0
        t_out = t2
    else:
        straight_into_service=1
        t_out = t1

    return t_out, straight_into_service
