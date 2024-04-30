from typing import List, Tuple
import random
import math
import numpy as np
import scipy.special as ss
#from scipy.special import gdtr, gdtrix

def sample_pretac_delay(alpha: float, beta: float, a: float, h: float, x_bar: float) -> float:
    """
    Sample pre-tactical delay for given flight.

    Pretactical delay is Y - (a - h)
    where a is scheduled arrival time, h is time tactical uncertainty begins
    and Y is a Gamma distribution with shape and rate parameters alpha and beta.
    If alpha or beta parameters are not valid i.e. <= 0, then pre-tactical
    delay is a given constant x_bar (mean delay).

    See section 4 of 
    "A New Simheuristic Approach fo Stochastic Runway Scheduling" (2022) by Shone et al
    for more details.
    """
    if alpha > 0 and beta > 0:
        # Note that numpy is parameterised by shape=alpha and scale = 1/beta
        pretac_delay = np.random.gamma(alpha, 1/beta) - (a - h)
    else: # in this case there isn't a well defined Gamma distribution
        # In this case the pre-tactical delay is set equal to the average lateness rather than being sampled randomly.
        pretac_delay = x_bar
    return pretac_delay

def sample_cond_gamma(t: float, k: int) -> float:
    """
    Sample from a Gamma distribution with shape k and rate 1, conditional on being
    above the value t.
    """
    q = ss.gdtr(1, k, t) # gdtr is a fast function in scipy for evaluating cdf of gamma dist
    z = q + (1-q) * random.random() # Unform RN between q and 1
    return ss.gdtrix(1, k, z) # gdtrix is fast function for quantile


def Gamma_GetServ(k: int, Time_Sep: List[List[int]], rel_time: float, trav_time: float, prev_class: int, cur_class: int, tm: float, weather_state: int, w_rho: float) -> Tuple[float, int]:
    """
    Simulate time that a flight finished its service.

    This is used when flight is already in the queue.

    Arguments:
    ---------
    k: Erlang shape parameter for service time
    Time_Sep: separation times between different classes
    rel_time: time flight was released to queue
    trav_time: travel time between entering pool to runway
    prev_class: weight class of previous flight
    cur_class: weight class of current flight
    tm: current time
    weather_state: code for current state of weather (0, 1 or 2)
    w_rho: multiplier for service time in case of bad weather

    Returns:
    --------
    t_out: time flight is finished being served
    straight_into_service: indicates whether flight enters service immediately on joining queue
    """

    # JF Question: I think I don't understand outputs here - check with Rob
    # JF Question: shouldn't current time be updated here to time of previous service completion?
    # This is for ACs that are already in the queue but not yet in service

    t1 = rel_time + trav_time # time at which flight reaches runway - JF - ask Rob about this

    rate = k / (Time_Sep[prev_class][cur_class]/60)
    if weather_state == 1:
        rate *= 1/w_rho

    getserv = np.random.gamma(k, 1/rate) # service time

    t2 = tm + getserv

    # Case 1: time to arrive at runway before current time plus time to be serviced
    if t1 < t2:
        straight_into_service = 0 # JF Question: is this the wrong way round? Rob says yes - if aircraft arrives at runway before previous service is finished it has to wait
        t_out = t2
    else:
        straight_into_service = 1
        t_out = t1

    return t_out, straight_into_service

def Gamma_GetServ_Future(k: int, Time_Sep: List[List[int]], rel_time, serv_time, trav_time, prev_class, cur_class, tm, weather_state, w_rho: float):

    # This is for ACs that have not yet been added to the queue

    t1 = rel_time + trav_time

    rate=k/(Time_Sep[prev_class][cur_class]/60)
    if weather_state==1:
        rate*=1/w_rho #=0.5

    getserv = serv_time # Look up the stored Gamma(k,1) value, serv_time
    getserv *= 1/rate # Convert it to the correct scale

    t2 = tm + getserv # JF Question: how can this be

    if t1 < t2:
        straight_into_service = 0
        t_out = t2
    else:
        straight_into_service = 1
        t_out = t1

    return t_out, straight_into_service

def Gamma_Conditional_GetServ(k: int, Time_Sep: List[List[int]], trav_time, rel_time, sv_time, prev_class, cur_class, tm, weather_state, w_rho: float):

    # This is for the AC currently in service

    t1=rel_time+trav_time

    rate=k/(Time_Sep[prev_class][cur_class]/60)
    if weather_state==1:
        rate*=1/w_rho #=0.5

    cond_serv = sample_cond_gamma(rate*(tm-sv_time), k)
    #cond_serv = sample_cond_gamma(rate*(tm-sv_time), gamma_cdf) # old code
    cond_serv *= 1/rate #Convert back to the correct scale

    t2=sv_time+cond_serv

    #max_t=max(t1,t2)

    if t1<t2:
        straight_into_service=0
        t_out=t2
    else:
        straight_into_service=1
        t_out=t1

    return t_out,straight_into_service

def gamma_cond_exp(t, alpha, beta):

    # Site: https://stats.stackexchange.com/questions/338378/closed-form-conditional-expectation-of-gamma-distributed-variable

    # This is assuming mean is alpha/beta, variance is alpha/(beta^2)
    # t is the amount of time that the service has already been in progress

    denom=0
    for j in range(alpha):
        denom+=((beta*t)**j)/math.factorial(j)
    num=((beta*t)**alpha)/math.factorial(alpha)
    x=(alpha/beta)*(1+num/denom)

    return x
