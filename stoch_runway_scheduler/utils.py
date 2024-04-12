from typing import List
from enum import Enum
from dataclasses import dataclass
import csv
import math
import random

#0: not ready yet (arrival), 1: in arrival pool, 2: added to arrival queue, 3: not ready yet (departure), 4: in departure pool, 5: added to departure queue, 6: finished.
class FlightStatus(Enum):
    NOT_READY = 0
    IN_POOL = 1
    IN_QUEUE = 2
    DEP_NOT_READY = 3 # JF Question: should these final two be combined?
    FINISHED = 6

@dataclass
class FlightInfo:
    status: FlightStatus
    ac_class: int 
    ps_time: float # 2 pre-scheduled time (plus pre-tactical delay)
    eta: float # 3 latest ETA
    release_time: float # 4 time at which aircraft is released from pool
    enters_service: float # 5 the time at which aircraft enters service
    travel_time: float # 6 travel time (generated in advance), # JF Question note from entering pool to runway? Update_ETAs suggests this is the case
    service_rns: List[float] # 7 list of random numbers used to calculate service times - JF NOTE: NOT USED!
    service_time: float # 8 actual service time s1+Z2 (worked out after class information is known)
    pool_time: float # 9 actual time that they join the pool (generated in advance)
    passenger_weight: float # g_i
    travel_time_indicator: bool # 11 indicator to show whether or not the AC's travel time has already been completed - used in Update_ETAs
    weather_state: int # 12 weather state at the time of release # JF Note - could possibly be removed with refactoring
    counter: int # counter
    qp: float # Perm only
    pred_cost: float # 15 predicted total cost at time of release
    service_completion_time: float # actual service completion time
    sched_dep_time: float # scheduled departure time - used?
    orig_sched_time: float # original pre-scheduled arrival time before adding pre-tactical delay
    flight_id: int

@dataclass
class SequenceInfo:
    sequence: List[int]
    n_traj: int # number of trajectories sampled so far
    v: float # V_s^n in paper  (eqn 14) (performance measure)
    queue_probs: List[float] # stores "probability" that aircraft will immediately go into service at end of remaining travel time - upsilon in paper (page 18)
    # Rob thinks this upsilon is redundant as Rob set lambda to zero in paper which means that if one random trajectory where aircraft enters
    # service immediately, this triggers release from pool.
    w: float # W_s^n in paper (eqn 15)

def read_flight_data(data_fn: str, min_time: int, max_time: int, wiener_sig: float):
    Ac_class = [] # this will store the weight class for each aircraft
    Orig_Ps = [] # original pre-scheduled times of aircraft, before applying the pre-tactical delay
    Dep_Ps = []  # h (i.e. time at which tactical uncertainty begins)
    flight_id = []
    Alpha_Ps = [] # Alpha parameters for Gamma dist for sampling pretactical delays
    Beta_Ps = []  # Beta parameters
    late_means = []
    with open(data_fn, 'r') as csvfile: # data file includes the on-time performance data such as means, variance of lateness based on 1 year of historical info from FlightRadar [DON'T CHANGE THIS FILE]
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata = list(datareader)
        n_lines = len(inputdata)
        for i in range(1, n_lines): # start from 1 because there's a title row
            ps_time = float(inputdata[i][1]) # pre-scheduled time
            # When reading in the flight data from the data file we only want to include flights with a pre-scheduled time between 6AM (360 mins) and 2PM (840 mins), including 6AM but not including 2PM
            if ps_time >= min_time and ps_time < max_time:
                ## Read Data ##
                ac_class = int(inputdata[i][2]) # weight class
                flight_name = str(inputdata[i][3]) # flight number
                dep_time = int(inputdata[i][4]) # departure time from the origin airport
                # scheduled flight dur, i.e. difference between scheduled departure and arrival time
                # JF Question: This doesn't seem to be used at the moment?
                sched_dur = int(inputdata[i][5]) 
                lateness_mn = float(inputdata[i][6]) # mean lateness based on historical data
                lateness_var = float(inputdata[i][7]) # variance of lateness based on historical data

                h = dep_time - 15 # JF Question: this 15 probably shouldn't be hard-coded - I think this is q_i in the paper? Probably yes

                ## Transform times ##
                # HERE WE RE-SCALE TIME SO THAT TIME '6AM' IS COUNTED AS TIME (ZERO+60).
                # WE START SIMULATING FROM TIME ZERO, I.E. AN HOUR BEFORE 6AM.
                h += -min_time+60
                ps_time += -min_time+60
                
                # The equations for xi_bar, si2, h_i, alpha and beta below are for calculating the parameters of the gamma distribution 
                # used for the pre-tactical delay. Details of this method are in Section 4 of the paper.
                xibar = ps_time + lateness_mn
                si2 = lateness_var 

                alpha = ((xibar-h)**2)/(si2-(wiener_sig**2)*(xibar-h)) # Can be negative via denominator
                beta = (xibar-h)/(si2-(wiener_sig**2)*(xibar-h))

                Ac_class.append(ac_class)
                Orig_Ps.append(ps_time) # original pre-scheduled time (as opposed to scheduled time following pre-tactical delay)
                flight_id.append(flight_name)
                Dep_Ps.append(h)
                Alpha_Ps.append(alpha)
                Beta_Ps.append(beta)
                late_means.append(lateness_mn)
            elif ps_time >= max_time: # Indicates that we have got to the end of the set of flights scheduled to arrive by 2PM
                break
    return flight_id, Ac_class, Orig_Ps, Dep_Ps, Alpha_Ps, Beta_Ps, late_means


def weather(tm, wlb, wub): # wlb is starting time for bad weather period, wub is ending time

    if wlb < wub:
        if tm < wlb:
            get_weather_state=0
        elif tm < wub:
            get_weather_state=1
        else:
            get_weather_state=2
    else:
        get_weather_state=0

    return get_weather_state

def getcost(ps_time, pool_time, trav_time, landing_time, pax_weight, thres1, thres2, lam1: float, lam2: float):

    cost=0
    #lam1=0.5 #weight for punctuality
    #lam2=0.5 #weight for queueing HMMM

    if landing_time > ps_time+thres1:
        cost+=lam1*pax_weight*(landing_time-(ps_time+thres1))**2

    if landing_time > pool_time+trav_time+thres2:
        cost+=lam2*pax_weight*(landing_time-(pool_time+trav_time+thres2))**2

    return cost
