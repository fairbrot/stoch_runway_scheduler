from typing import List
import logging
from .utils import FlightInfo, Cost, FlightStatus
from .weather import WeatherProcess, WeatherStatus
from .sequence import SequenceInfo
from .separation import StochasticSeparation, landing_time
from .gamma import gamma_cond_exp

# JF: This is the main deterministic heuristic
# Name may not be best choice
def Genetic_determ(Ac_Info: List[FlightInfo], Arr_Pool: List[int], Arr_NotReady: List[int], 
                    Ac_queue: List[int], tm: float, NoA: int, sep: StochasticSeparation,
                    prev_class: int, GA_Info: List[SequenceInfo], weather: WeatherProcess, tau: int, Max_LookAhead: int, cost_fn: Cost, 
                    basecost: float):

    stepthrough_logger = logging.getLogger("stepthrough")
    step_summ_logger = logging.getLogger("step_summ")

    ArrTime = [0] * NoA

    # Generate arrival and service time percentiles for AC not yet in queue
    ArrTime_Sorted = []

    for AC in Arr_Pool:
        ArrTime[AC] = Ac_Info[AC].pool_time
        ArrTime_Sorted.append([ArrTime[AC], AC])

    for AC in Arr_NotReady:
        ArrTime[AC]=max(0, Ac_Info[AC].eta - tau)
        ArrTime_Sorted.append([ArrTime[AC], AC])

    ArrTime_Sorted.sort(key=lambda x: x[0])

    stepthrough_logger.info('basecost is %.2f', basecost)
    stepthrough_logger.info('Generated results for ACs already in the queue are as follows:')
    stepthrough_logger.info('AC, Class, Time Sep, Release time, Travel time, Enters serv, Actual serv, Finish time, Pax weight, Cost')

     # The case where queue is empty needs fixing
     # This function needs to know time of last landing in order to set prev_landing here
    if len(Ac_queue) > 0:

        # Need to generate service times for AC already in the queue; first consider the customer in position 0
        AC = Ac_queue[0]
        Ac_Infoi = Ac_Info[AC]
        rel_time = Ac_Infoi.release_time
        sv_time = Ac_Infoi.enters_service
        cur_class = Ac_Infoi.ac_class

        t1 = Ac_Infoi.eta

        # Get the conditional expectation of service time based on service time elapsed so far
        cond_tm = sep.expected_conditional_seperation(tm - sv_time, prev_class, cur_class,
                                                      Ac_Infoi.weather_state)
        t2 = sv_time + cond_tm

        queue_complete = max(t1,t2)
        basecost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, tau, queue_complete, Ac_Infoi.passenger_weight)
        perm_prev_class = cur_class

        # Now consider the rest of the customers in the queue
        for AC in Ac_queue:
            Ac_Infoi = Ac_Info[AC]
            rel_time = Ac_Infoi.release_time
            cur_class = Ac_Infoi.ac_class

            t1 = Ac_Infoi.eta
            weather_state = weather(rel_time) # weather(queue_complete)
            t2 = queue_complete + sep.expected_separation(prev_class, cur_class, weather_state)
            queue_complete = max(t1, t2)

            perm_prev_class = cur_class
            basecost += cost_fn(Ac_Infoi.orig_sched_time,Ac_Infoi.pool_time,tau,queue_complete,Ac_Infoi.passenger_weight)
    else: 
        queue_complete = tm # not right
        perm_prev_class = prev_class

    stored_prev_class = perm_prev_class

    # Try all the sequences in the population
    for info in GA_Info:
        permcost = basecost
        latest_tm = tm
        perm_prev_class = stored_prev_class
        perm_queue_complete = queue_complete

        perm = info.sequence

        # JF Question: why would this not be len(perm)?
        xi_list = [] # for each flight indicates whether or not it goes straight to service - called xi in paper
        for AC in perm:
            Ac_Infoi = Ac_Info[AC]
            perm_class = Ac_Infoi.ac_class
            reltime = max(latest_tm, ArrTime[AC])

            t1 = reltime + tau
            # JF Question - should we use time start servicing here?
            # Like begin_serv in Genetic?
            exp_serv = sep.expected_separation(perm_prev_class, perm_class, weather(rel_time))
            t2 = perm_queue_complete + exp_serv

            AC_FinishTime, straight_into_service = landing_time(perm_queue_complete, exp_serv,
                                                                rel_time, tau)
            xi_list.append(straight_into_service)
            permcost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC], tau, AC_FinishTime, Ac_Infoi.passenger_weight)
            
            latest_tm = reltime
            perm_queue_complete = AC_FinishTime
            perm_prev_class = perm_class
        

        info.add_observation(permcost, xi_list)

    GA_Info.sort(key=lambda x: x.v)

    Ac_added = []
    # JF Note: ask Rob why no additional checks on flights
    assert len(GA_Info) > 0
    perm = GA_Info[0]
    for AC in perm:
        if Ac_Info[AC].status != FlightStatus.IN_POOL:
            break
        Ac_added.append(AC)

    return Ac_added
