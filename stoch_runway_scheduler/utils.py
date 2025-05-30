from typing import List, Tuple
from enum import Enum
from dataclasses import dataclass
import csv
import math
import random

# 0: not ready yet (arrival), 1: in arrival pool, 2: added to arrival queue, 3: not ready yet (departure), 4: in departure pool, 5: added to departure queue, 6: finished.
class FlightStatus(Enum):
    NOT_READY = 0
    IN_POOL = 1
    IN_QUEUE = 2
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
    service_rns: List[float] # 7 list of random numbers used to calculate service times
    service_time: float # 8 actual service time s1+Z2 (worked out after class information is known)
    pool_time: float # 9 actual time that they join the pool (generated in advance)
    passenger_weight: float # g_i
    travel_time_indicator: bool # 11 indicator to show whether or not the AC's travel time has already been completed - used in Update_ETAs
    weather_state: int # 12 weather state at the time of release # JF Note - could possibly be removed with refactoring
    counter: int # counter
    pred_cost: float # 15 predicted total cost at time of release
    service_completion_time: float # actual service completion time
    sched_dep_time: float # scheduled departure time - used?
    orig_sched_time: float # original pre-scheduled arrival time before adding pre-tactical delay
    flight_id: int


def read_flight_data(data_fn: str, min_time: int, max_time: int, wiener_sig: float) -> Tuple[List[str], List[int], List[int], List[int], List[float], List[float], List[float]]:
    """
    Reads flight data from CSV file.

    Each flight corresponds to a row of the table and has the following attributes:

    Arr Time (int): Pre-scheduled arrival time at airport. 
    Class (int): Weight class of aircraft (0 = H (Heavy), 1 = U, 2 = M, 3 = L)
    Flight number (str): IATA flight designator (e.g. BA032)
    Dep time (int): Departure time from origin aiport (< Arr Time)
    Flight time (int): duration of flight
    Pretac_mean (float): mean of pretactical data (i.e. E[A] - a)
    Pretac_var (float): variance of arrival time (or pretactical delay) (i.e. Var[A])

    Times are given in minutes after midnight (e.g. 360 is 6AM).
    Durations are also given in minutes.

    Function processes these inputs and outputs flight attributes as lists.
    In particular, only flights which occur within given time interval
    are included, and for each of these flights, the times are adjusted
    so that the beginning of the time horizon is 0, and parameters
    for a Gamma random variable for each flight, Y, which is used to calculate
    pre-tactical delays, are calculated. 
    Pre-tactical delays in particular are given by Y - (a - h).
    These parameters are set so that
    the overall mean and variance of the actual arrival time (calculated empirically)
    match those simulated by pre-tactical and tactical delays.
    See section 4 of 
    "A New Simheuristic Approach fo Stochastic Runway Scheduling" (2022) by Shone et al
    for more details.

    Arguments:
    ---------
    data_fn: path of CSV to load
    min_time: Start time of time horizon under consideration
    max_time: End time of time horizon under consideration
    wiener_sig: standard deviation of Brownian motion used in tactical delays

    Returns:
    ---------
    flight_id: IATA flight designators for each flight
    Ac_class: Weight classes
    Orig_Ps: Adjusted scheduled arrival times of flights at destination airport
    Dep_Ps: Adjusted times at which tactical uncertainty begins (i.e. h)
    Alpha_Ps: Shape parameters of Gamma distributions used to calculate pre-tactical delay
    Beta_Ps: Rate parameters of Gamma distributions used to calculate pre-tactical delay
    late_means: Average lateness values for arriving
    """
    # JF Question: ask Rob to check accuracy of above docstring
    # JF Question: is a column fo Flight time necessary? This could be inferred from Dep and Arr times

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
                sched_dur = int(inputdata[i][5]) # JF Question: This doesn't seem to be used at the moment?
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
                xibar = ps_time + lateness_mn # I.e. mean arrival time
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



@dataclass
class Cost:
    thres1: float # gamma^S in paper - tolerance for schedule delay i.e. if this is 15 minutes, then delays less than 15 minutes don't induce any penalty
    thres2: float # gamma^W in paper - tolerance for extra airborne delay
    lam1: float # theta^S in paper - relative weight on schedule delay
    lam2: float # theta^W in paper - relative weight on extra airborne delay

    def __call__(self, ps_time: float, pool_time: float, trav_time: float, landing_time: float, pax_weight: float) -> float:
        """Calculates cost for a given flight which is a weighted sum of penalties for schedule delay and extra airborne delay.

        Arguments:
        ----------
        ps_time: scheduled arrival time
        pool_time: time at which aircraft joins pool
        trav_time: travel time between joining pool and reaching runway
        landing_time: actual time aircraft lands (after serviced in queue)
        pax_weight: weighting dependent on importance of flight (e.g. may depend on passengers)
        """
        C_S = max(landing_time - (ps_time + self.thres1), 0)**2
        # JF Question: this needs explaining to me - formula seems a bit different to what used in paper
        C_W = max(landing_time - (pool_time + trav_time + self.thres2), 0)**2
        return pax_weight * (self.lam1 * C_S + self.lam2 * C_W)