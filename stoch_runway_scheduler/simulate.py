from typing import List, TextIO, Tuple
import logging
import math
import random
import numpy as np
from .utils import FlightStatus, FlightInfo, Cost
from .sequence import SequenceInfo
from .weather import WeatherProcess, StochasticWeatherProcess
from .separation import StochasticSeparation, landing_time
from .annealing_cost import Annealing_Cost


def simulate_sequences(GA_Info: List[SequenceInfo], tm: float, Ac_Info: List[FlightInfo], Ac_queue: List[int], tau: float, sep: StochasticSeparation,
                        weather: StochasticWeatherProcess, wiener_sig: float, prev_class: int, cost_fn: Cost) -> Tuple[List[float], List[List[int]]]:
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
    prev_ld = tm if not Ac_queue else Ac_Info[Ac_queue[0]].enters_service

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
        trav_time = Ac_Infoi.travel_time if tm >= Ac_Infoi.eta else np.random.wald(tau, (tau/wiener_sig)**2)
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
        match info.status:
            case FlightStatus.NOT_READY:
                sched = int(round(Ac_Info[AC].eta-(tm+tau),1)) # JF Question: what is this?
                if sched<=0:
                    arr_time = tm
                else:
                    arr_time = np.random.wald(sched,(sched/wiener_sig)**2) + tm
            case FlightStatus.IN_POOL:
                arr_time = max(0, Ac_Info[AC].eta - tau) # JF Question: Time aircraft arrives in Pool - could this not be replaced with Ac_Info[AC].pool_time?
            case _:
                raise RuntimeError(f'Aircraft {AC} is in a sequence but has status {info.status}')
        trav_time = np.random.wald(tau, (tau / wiener_sig)**2)
        serv_time=sep.sample_normalized_separation()
        ac_attrs[AC] = (arr_time, trav_time, serv_time)

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
            arr_time, trav_time, serv_time = ac_attrs[AC]
            Ac_Infoi = Ac_Info[AC]
            perm_class = Ac_Infoi.ac_class
            reltime = max(latest_tm, arr_time)
            weather_state = weather_sample(reltime)
            
            min_sep = sep.sample_separation(perm_prev_class, perm_class, weather_state, norm_service_time = serv_time)
            AC_FinishTime, straight_into_service = landing_time(prev_ld, min_sep, reltime, trav_time)
            xi_list.append(straight_into_service)
            permcost += cost_fn(Ac_Infoi.orig_sched_time, arr_time, trav_time, AC_FinishTime, Ac_Infoi.passenger_weight)
            #latest_tm = reltime # JF Question - I don't think this needs updating - all flights should join queue now or when they enter pool

            prev_ld = AC_FinishTime
            perm_prev_class = perm_class
        costs.append(permcost)
        xi_lists.append(xi_list)

    return costs, xi_lists

def generate_trajectory(Dep_time: float, Ps_time: float, tau: int, wiener_sig: float, freq: int):
    """
    Generates a trajectory for an aircraft.

    Parameters
    ----------
    Dep_time: time at which Brownian motion starts being used to predict ETA (called h in paper)
    Ps_time: pre-scheduled arrival time at destination airport plus pre-tactical delay
    tau: flight is considered to have joined pool when current time >= ETA - tau
    wiener_sig: standard deviation of Brownian motion
    freq: frequency (per minute) at which trajectories are updated

    Returns
    -------
    pool_arr_time: float
        time at which flight arrives in pool
    travel_time: float
        time to travel from pool threshold to runway
    brown_motion: List[float]
        Brownian motion. Element j is the ETA at time j/freq.
    """

    brown_motion = []
    # For flights with scheduled departure less than 0 we need to initially simulate where BM would be at time 0
    if Dep_time < 0:
        ETA = random.gauss(Ps_time, math.sqrt(0-Dep_time)*wiener_sig) # Update the latest ETAs for ACs that already had their dep time before time zero
    else:
        ETA = Ps_time # ETA = pre-scheduled time
    brown_motion.append(ETA)

    # i.e. when ETA is within 30 minutes of current time (0)
    # JF note: perhaps add a check for aircraft already arriving at runway
    if 0 >= ETA-tau:
        pool_arr_time = 0
        chk = 1
        ETA = tau # JF Question - is this right?
    else:
        chk = 0

    j = 0
    while True:
        j += 1 # step forward in increments of 1/freq minutes
        if j > Dep_time*freq: # only update ETA if we've gone beyond the AC's departure time
            ETA = random.gauss(ETA, 0.1*wiener_sig)
        brown_motion.append(ETA)
        if j/freq >= ETA-tau and chk == 0:
            pool_arr_time = round(j/freq,2) # pool arrival time - j/freq is the 'current time'
            chk = 1
        elif j/freq >= ETA and chk == 1:
            runway_time = round(j/freq,2)
            # Plane arrives at runway
            travel_time = runway_time - pool_arr_time # travel time between entering pool and arriving at runway
            chk = 2
            break

    return pool_arr_time, travel_time, brown_motion

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
def Posthoc_Check(seq, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process: WeatherProcess, output, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], cost_fn: Cost):

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


def round_down(tm: float, freq: int) -> float:
    """Rounds tm down to beginning of interval defined by freq.

    For example, if freq is 100, then we round down for step size of 0.01,
    that is 0.32453 would become 0.32.
    """
    int_size = 1/freq
    i = tm // int_size
    return i*int_size


def Update_ETAs(Ac_Info: List[FlightInfo], Arr_NotReady: List[int], Ac_queue: List[int], 
                tm: float, Brown_Motion: List[List[float]], Arr_Pool: List[int], tau: float,
                freq: int):
    
    # JF Question: this only updates ETAs for aircraft not ready, which join the pool, or are in the queue.
    # It does not seemingly update the ETAs for flights aleady in the pool.
    
    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    step_new_logger = logging.getLogger('step_new')

    i = 0
    # Updates flights not aleady in Pool
    to_remove = []
    for i, AC in enumerate(Arr_NotReady):
        Ac_Infoi = Ac_Info[AC]
        if tm >= Ac_Infoi.pool_time:
            Arr_Pool.append(AC)
            to_remove.append(i)
            Ac_Infoi.status = FlightStatus.IN_POOL
            Ac_Infoi.eta = Ac_Infoi.pool_time + tau # JF Question: is this necessary?

            msg = '* Added aircraft '+str(AC)+' to the arrival pool at time '+str(tm)+' (new readiness time is '+str(Ac_Infoi.pool_time+tau)+')'+'\n'+'\n'
            stepthrough_logger.info(msg)
            step_summ_logger.info(msg)
            step_new_logger.info(msg)
        else:
            Ac_Infoi.eta = Brown_Motion[AC][int(tm*freq)]
    for i in reversed(to_remove):
        Arr_NotReady.pop(i)

    # Updates flights which are in the queue
    for i, AC in enumerate(Ac_queue):
        Ac_Infoi = Ac_Info[AC]
        if not Ac_Infoi.travel_time_indicator:
            rel_time = Ac_Infoi.release_time
            trav_so_far = tm - rel_time # amount of time spent travelling to the runway so far
            # JF Question: why round? Should this be round down? Is this related to freq?
            # I think we need rounding down to the beginning of interval defined by freq
            rounded_trav_so_far = round_down(trav_so_far, freq)
            if rounded_trav_so_far >= Ac_Infoi.travel_time:
                msg = '* Aircraft '+str(AC)+' has finished its travel time at '+str(tm)+'; still waiting for service time.'+'\n'+'\n'
                stepthrough_logger.info(msg)
                step_summ_logger.info(msg)
                Ac_Infoi.travel_time_indicator = True
            else:
                Ac_Infoi.eta = Brown_Motion[AC][int((Ac_Infoi.pool_time + rounded_trav_so_far) * freq)]

def Update_Stats(tm: float, AC: int, Ac_Info: List[FlightInfo], Ac_queue: List[int], real_queue_complete: float, weather_process: WeatherProcess, latest_class, next_completion_time, sep: StochasticSeparation, SubPolicy: str, counter: int):
    """
    Updates various states when after flight is released into queue.

    Arguments
    ---------
    real_queue_complete: time last aircraft to join queue was or will be serviced
    weather_process: object which can be queried for weather state at given time
    latest_class: class of previous aircraft to join queue
    next_completion_time: time of next service in queue
    k: Erlang service parameter
    Time_Sep: time separation array
    w_rho: bad weather multiplier
    Subpolicy:
    counter:

    Returns
    -------
    real_queue_complete: time aircraft currently joining queue will be finished being served (function argument is updated)
    next_completion_time: time of next service in queue (function argument is updated)
    latest_class: class of flight just added (function argument is updated)
    """
    # Function which produces statistics about a flight which has just been released - also calculates
    # when this flight will be finished being serviced

    Ac_Infoi = Ac_Info[AC]

    Ac_Infoi.status = FlightStatus.IN_QUEUE
    release_time = tm # release time
    #begin_serv = max(release_time, real_queue_complete) # JF Question: is this a mistake - I think we should use below
    begin_serv = real_queue_complete # JF: i.e. time previous aircraft lands
    cur_class = Ac_Infoi.ac_class

    # Weather state at time flight joins queue
    weather_state = weather_process(release_time)

    # Separation/Service time for flight to land (once it is being processed)
    actual_serv = sep.sample_separation(latest_class, cur_class, weather_state, norm_service_time=Ac_Infoi.service_rns)

    trav_time = Ac_Infoi.travel_time
    finish_time, straight_into_service = landing_time(begin_serv, actual_serv, release_time, trav_time)

    real_queue_complete = finish_time

    Ac_Infoi.release_time = release_time
    Ac_Infoi.enters_service = begin_serv
    Ac_Infoi.weather_state = weather_state
    Ac_Infoi.service_time = actual_serv
    Ac_Infoi.service_completion_time = finish_time

    latest_class = cur_class

    # I.e. if flight just added to the queue is only flight in queue
    if len(Ac_queue) == 1:
        next_completion_time = finish_time

    return real_queue_complete, next_completion_time, latest_class

def Serv_Completions(Ac_Info, Ac_queue, prev_class, totserv, Ac_finished, tm, next_completion_time, cost_fn: Cost, f: TextIO, SubPolicy, rep, Left_queue):

    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    # print('Entered Serv_Completions')
    # print('ESC tm: '+str(tm)+' next_completion_time: '+str(next_completion_time))

    arr_cost=0
    dep_cost=0

    j=0

    while len(Ac_queue) > 0:

        AC = Ac_queue[0]
        Ac_Infoi = Ac_Info[AC]

        finish_time = Ac_Infoi.service_completion_time
        current_class = Ac_Infoi.ac_class

        #print('finish_time: '+str(finish_time))

        if tm >= finish_time: #release_time+trav_time and phase==k:
            #print('* Service phase '+str(phase)+' completed for aircraft '+str(Ac_queue[0])+' at time '+str(tm+delta))
            Ac_finished[AC] = finish_time
            #print('* Service completion finished for aircraft '+str(AC))
            stepthrough_logger.info('* Service completion finished for aircraft '+str(AC)+'\n'+'\n')
            step_summ_logger.info('* Service completion finished for aircraft '+str(AC)+'\n'+'\n')
            if Ac_Infoi.status == FlightStatus.IN_QUEUE:
                Ac_Infoi.status = FlightStatus.FINISHED
                arr_cost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight)
                totserv+=1
            else: # JF Question: is this clause needed?
                Ac_Infoi.status = FlightStatus.FINISHED
                arr_cost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight)

            #f.write(str(SubPolicy)+','+str(rep)+','+str(AC)+','+str(Ac_Infoi.flight_id)+','+str(prev_class)+','+str(current_class)+','+str(Time_Sep[prev_class][current_class]/60)+','+str(Ac_Infoi.orig_sched_time)+','+str(Ac_Infoi.ps_time)+','+str(Ac_Infoi.pool_time)+','+str(Ac_Infoi.release_time)+','+str(Ac_Infoi.travel_time)+','+str(Ac_Infoi.weather_state)+','+str(Ac_Infoi.enters_service)+','+str(Ac_Infoi.service_time)+','+str(Ac_Infoi.service_completion_time)+','+str(max(0,finish_time-(Ac_Infoi.ps_time+cost_fn.thres1)))+','+str(finish_time-(Ac_Infoi.pool_time+Ac_Infoi.travel_time))+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight))+',')

            f.write(str(Ac_Infoi.pred_cost)+',')
            f.write('\n')

            prev_class = current_class

            Left_queue.append(AC)
            Ac_queue.remove(AC)

            print(str(SubPolicy)+' totserv: '+str(totserv))

            if len(Ac_queue)>0:
                New_AC=Ac_queue[0]
                next_completion_time=Ac_Info[New_AC].service_completion_time

        else:
            break

    return arr_cost, dep_cost, totserv, prev_class, Ac_finished, next_completion_time
