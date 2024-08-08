from typing import List, TextIO, Tuple
import logging
from collections import deque
from .utils import FlightStatus, FlightInfo, Cost
from .sequence import SequenceInfo
from .weather import WeatherProcess, StochasticWeatherProcess
from .trajectory import StochasticTrajectory
from .separation import StochasticSeparation, landing_time
from .annealing_cost import Annealing_Cost


def simulate_sequences(GA_Info: List[SequenceInfo], tm: float, Ac_Info: tuple[FlightInfo], Ac_queue: deque[int], trajecs: tuple[StochasticTrajectory], sep: StochasticSeparation,
                        weather: StochasticWeatherProcess, prev_ld: float, prev_class: int, cost_fn: Cost) -> Tuple[List[float], List[List[int]]]:
    """
    Simulates landing times and associated costs of a set of sequences using common random numbers.

    Arguments:
    ----------
    GA_Info: list of current sequences
    tm: current time
    Ac_Info: information for all flights at current time
    Ac_queue: list of aircraft currently in the queue
    tau: threshold for reaching pool (flight in pool when ETA - tm <= tau)
    sep: object for sampling normalized separation times
    weather: object for generating random weather process conditional on state at current time
    wiener_sig: standard deviation of Brownian motion
    prev_ld: time of previous landing (service completion)
    prev_class: class of previous aircraft to land
    cost_fn: function for calculating costs of flight

    Returns:
    -------
    costs: costs of each sequence
    xi_lists: for each sequence, indicators of whether each flight went straight into service
    """
    weather_sample = weather.sample_process(tm)
    new_cost = 0
    perm_prev_class = prev_class

     # The case where queue is empty needs fixing - the same issue existed before refactoring Genetic
     # This function needs to know time of last landing in order to set prev_landing here

    # Simulate current queue
    for i, AC in enumerate(Ac_queue):
        # Need to generate service times for AC already in the queue; first consider the customer in position 0
        AC = Ac_queue[0]
        Ac_Infoi = Ac_Info[AC]
        # JF: this code was used previously for case where tm < eta
        # z=int(random.randrange(1,999))
        # sched=int(10*round(Ac_Infoi.eta-tm,1))
        # trav_time=wiener_cdf[sched][z]
        # sched = int(round(info.eta - tm, 1)) # similar to NOT_READY case
        trav_time = trajecs[AC].simulate_travel_time(tm, Ac_Infoi)
        if i == 0:
            min_sep = sep.sample_conditional_separation(tm - prev_ld, prev_class, Ac_Infoi.ac_class, Ac_Infoi.weather_state)
        else:
            # JF Question: previous code was sampling whether at release (which is in past) which didn't seem right
            min_sep = sep.sample_separation(perm_prev_class, Ac_Infoi.ac_class, Ac_Infoi.weather_state)

        landing_complete, straight_into_service = landing_time(prev_ld, min_sep, Ac_Infoi.release_time, trav_time)
        new_cost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, trav_time, landing_complete, Ac_Infoi.passenger_weight)

        prev_ld = landing_complete
        perm_prev_class = Ac_Infoi.ac_class

    stored_prev_class = perm_prev_class
    stored_prev_ld = prev_ld

    seq_ACs = set(ac for seq in GA_Info for ac in seq.sequence)
    ac_attrs = dict()

    # Generate arrival, travel and normalized service times for flights in sequences
    for AC in seq_ACs:
        info = Ac_Info[AC]
        pool_time = trajecs[AC].simulate_pool_time(tm, info)
        trav_time = trajecs[AC].simulate_travel_time(tm, info)
        serv_time = sep.sample_normalized_separation()
        ac_attrs[AC] = (pool_time, trav_time, serv_time)

    costs = []
    xi_lists = []
    for info in GA_Info:
        # stepthrough_logger.info('Now trying sequence %s', info.sequence)
        # stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost = new_cost
        latest_tm = tm
        perm_prev_class = stored_prev_class
        prev_ld = stored_prev_ld

        perm = info.sequence
        # no_ACs = min(Max_LookAhead, len(perm)) # JF Question: can this not just be len(perm)? Also, don't we want to evaluate full length of sequence?
        xi_list = [] # for each flight indicates whether or not it goes straight to service - called xi in paper
        for AC in perm:
            pool_time, trav_time, serv_time = ac_attrs[AC]
            Ac_Infoi = Ac_Info[AC]
            perm_class = Ac_Infoi.ac_class
            reltime = max(latest_tm, pool_time)
            weather_state = weather_sample(reltime)
            
            min_sep = sep.sample_separation(perm_prev_class, perm_class, weather_state, norm_service_time = serv_time)
            AC_FinishTime, straight_into_service = landing_time(prev_ld, min_sep, reltime, trav_time)
            xi_list.append(straight_into_service)
            permcost += cost_fn(Ac_Infoi.orig_sched_time, pool_time, trav_time, AC_FinishTime, Ac_Infoi.passenger_weight)
            #latest_tm = reltime # JF Question - I don't think this needs updating - all flights should join queue now or when they enter pool

            prev_ld = AC_FinishTime
            perm_prev_class = perm_class
        costs.append(permcost)
        xi_lists.append(xi_list)

    return costs, xi_lists

def Calculate_FCFS(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process: WeatherProcess, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], cost_fn: Cost):

    tm=0

    totserv=0
    totcost=0
    prev_class=4
    queue_complete=0
    weather_state=0

    FCFS_Seq = [0]*NoA
    for i in range(NoA):
        FCFS_Seq[i]=ArrTime_Sorted[i][1]

    FCFS_cost = Annealing_Cost(FCFS_Seq, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process, 0, NoA, w_rho, k, Time_Sep, cost_fn)

    return FCFS_cost

# JF Question: what does this do?
def Posthoc_Check(seq: list[int], Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process: WeatherProcess, output, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], cost_fn: Cost):
    # seq: order in which flights were served
    perm=seq
    perm_cost=0
    latest_tm=0
    perm_prev_class=4
    perm_queue_complete=0
    perm_weather_state=0
    j=0

    while j<NoA:

        AC = perm[j]
        Ac_Infoi = Ac_Info[AC]
        release_time = Ac_Infoi.release_time # max(latest_tm, ArrTime[AC][0])
        #print('j: '+str(j)+' AC: '+str(AC)+' release_time: '+str(release_time))
        trav_time=Ac_Infoi.travel_time
        perm_class=Ac_Infoi.ac_class
        begin_serv=max(release_time,perm_queue_complete)
        perm_weather_state=weather_process(release_time) # JF Question: why not weather_process(begin_serv)?

        if perm_weather_state==1:
            ws=1/w_rho
        else:
            ws=1
        rate=ws*k/(Time_Sep[perm_prev_class][perm_class]/60)
        serv=ServTime[AC]/rate


        t1=release_time+trav_time
        t2=perm_queue_complete+serv
        perm_queue_complete=max(t1,t2)

        perm_cost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC][0], trav_time, perm_queue_complete, Ac_Infoi.passenger_weight)

        if output==1:
            print('AC: '+str(AC)+' class: '+str(perm_class)+' release_time: '+str(release_time)+' trav_time: '+str(trav_time)+' begin_serv: '+str(begin_serv)+' t1: '+str(t1)+' t2: '+str(t2)+' finish time: '+str(perm_queue_complete)+' weather state: '+str(perm_weather_state)+' pax weight: '+str(Ac_Infoi.passenger_weight)+' cost: '+str(cost_fn(Ac_Infoi.ps_time,ArrTime[AC][0],trav_time,perm_queue_complete,Ac_Infoi.passenger_weight)))

        latest_tm=release_time
        perm_prev_class=perm_class

        j+=1

    return perm_cost