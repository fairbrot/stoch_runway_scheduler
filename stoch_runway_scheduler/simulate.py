from typing import List, TextIO, Tuple
import logging
import math
import random
import numpy as np
from .utils import FlightStatus, FlightInfo, Cost
from .weather import WeatherProcess
from .annealing_cost import Annealing_Cost


def simulate_flight_times(tm: float, Ac_Info: List[FlightInfo], tau: float, k: int, wiener_sig: float) -> Tuple[List[float], List[float], List[float]]:
    """
    Simulates arrival in pool times and travel times conditional on situation at time tm.

    Arguments:
    ----------
    tm: current time
    Ac_Info: information for all flights at current time
    tau: threshold for reaching pool (flight in pool when ETA - tm <= tau)
    k: Erlang shape parameter for service times
    wiener_sig: standard deviation of Brownian motion

    Returns:
    -------
    ArrTime: List of arrival times in Pool (0 for any flights already served)
    Trav_Time: List of travel times between entering pool and reaching runway (0 for any flights already served)
    ServTimes: List of standardized service times (0 for any flights in queue or already served)
    """
    NoA = len(Ac_Info)

    ArrTime = [0]*NoA # pool time
    Trav_Time = [0]*NoA # travel time between pool and runway
    ServTime = [0]*NoA # standardized service time

    for AC, info in enumerate(Ac_Info):
        match info.status:
            case FlightStatus.NOT_READY:
                sched = int(round(Ac_Info[AC].eta-(tm+tau),1)) # JF Question: what is this?
                if sched<=0:
                    ArrTime[AC]=tm
                else:
                    ArrTime[AC]=np.random.wald(sched,(sched/wiener_sig)**2) + tm
                Trav_Time[AC]=np.random.wald(tau,(tau/wiener_sig)**2)
                ServTime[AC]=np.random.gamma(k,1)

            case FlightStatus.IN_POOL:
                ArrTime[AC] = max(0, Ac_Info[AC].eta - tau) # JF Question: Time aircraft arrives in Pool - could this not be replaced with Ac_Info[AC].pool_time?
                Trav_Time[AC] = np.random.wald(tau, (tau / wiener_sig)**2)
                ServTime[AC] = np.random.gamma(k, 1)

            case FlightStatus.IN_QUEUE:
                # RS: this block of code is probably needed but wasn't included in the 5000 experiments for the paper
                # JF Question : current code seemed incorrect - I've tried to adapt below to use np.wald. Please check this.
                if tm >= info.eta:
                	Trav_Time[AC] = info.travel_time # travel time has already finished
                else:
                	# z=int(random.randrange(1,999))
                	# sched=int(10*round(Ac_Infoi.eta-tm,1))
                	# trav_time=wiener_cdf[sched][z]
                    # sched = int(round(info.eta - tm, 1)) # similar to NOT_READY case
                    # trav_time = np.random.wald(sched,(sched/wiener_sig)**2) # JF Question: is this right?
                    Trav_Time[AC] = np.random.wald(tau, (tau/wiener_sig)**2) # JF: this was previously used
                    ArrTime[AC] = info.pool_time
                
    return ArrTime, Trav_Time, ServTime

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

def Get_Actual_Serv(AC: int, prev_class: int, cur_class: int, weather_state: int, k: int, Time_Sep: List[List[int]], Ac_Info: List[FlightInfo], w_rho: float) -> float:
    """
    Calculates service time of landing aircraft given aircraft's class and that of previous aircraft, and state of weather.
    
    Arguments
    ---------
    AC: aircraft index
    prev_class: class of previous aircraft to be added to the queue
    cur_class: class of aircraft currently joining the queue
    weather_state: type of weather (1 bad)
    k: Erlang service parameterr
    Time_Sep: time separation array
    Ac_Info: list of flight information
    w_rho: bad weather multiplier

    Returns:
    --------
    servtime: random service time for aircraft which is joining queue
    """
    # Samping queue service time
    # Service time also depends on state of weather
    # Ac_Info[AC].service_rns is pre-generated random number which is then scaled appropriately

    serv_percs=Ac_Info[AC].service_rns
    if weather_state==1:
        ws = 1/w_rho
    else:
        ws = 1

    # JF Question: does this depend on freq or conv_factor?
    rate = ws*k/(Time_Sep[prev_class][cur_class]/60)

    # Transformation causes serv_percs to go from [mean k, var k] to [mean e_{ij}, var e_{ij}^2/k]
    # See page 8 of paper (Section 3.1)
    servtime = serv_percs/rate

    return servtime

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

def Update_Stats(tm: float, AC: int, Ac_Info: List[FlightInfo], Ac_queue: List[int], real_queue_complete: float, weather_process: WeatherProcess, latest_class, Ov_GA_counter, next_completion_time, k: int, Time_Sep: List[List[int]], w_rho: float, SubPolicy: str, counter: int):
    """
    Updates various states when after flight is released into queue.

    Arguments
    ---------
    real_queue_complete: time last aircraft to join queue was or will be serviced
    weather_process: object which can be queried for weather state at given time
    latest_class: class of previous aircraft to join queue
    Ov_GA_counter:
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
    Ov_GA_counter:
    """
    # Function which produces statistics about a flight which has just been released - also calculates
    # when this flight will be finished being serviced

    Ac_Infoi = Ac_Info[AC]

    Ac_Infoi.status = FlightStatus.IN_QUEUE
    release_time = tm # release time
    begin_serv = max(release_time, real_queue_complete) # time that service begins
    cur_class = Ac_Infoi.ac_class

    # Weather state at time flight joins queue
    get_weather_state = weather_process(release_time) # JF Question: why not begin_serv here?

    # Service time for flight to land (once it is being processed)
    # JF Question: weather state might be wrong here - should this be weather at time flight begins service?
    actual_serv = Get_Actual_Serv(AC, latest_class, cur_class, get_weather_state, k, Time_Sep, Ac_Info, w_rho)
    trav_time = Ac_Infoi.travel_time

    t1 = release_time + trav_time
    t2 = begin_serv + actual_serv

    # time service finishes
    # JF Question: why is max needed here?
    finish_time = max(t1, t2)
    real_queue_complete = finish_time

    Ac_Infoi.release_time = release_time
    Ac_Infoi.enters_service = begin_serv
    Ac_Infoi.weather_state = get_weather_state
    Ac_Infoi.service_time = actual_serv
    Ac_Infoi.service_completion_time = finish_time

    Ac_Infoi.counter = Ov_GA_counter
    Ov_GA_counter = 0

    latest_class = cur_class

    # I.e. if flight just added to the queue is only flight in queue
    if len(Ac_queue) == 1:
        next_completion_time = finish_time

    return real_queue_complete, next_completion_time, latest_class, Ov_GA_counter

def Serv_Completions(Ac_Info, Ac_queue, prev_class, totserv, Ac_finished, tm, next_completion_time, cost_fn: Cost, f: TextIO, SubPolicy, rep, Time_Sep: List[List[int]], Left_queue):

    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    # print('Entered Serv_Completions')
    # print('ESC tm: '+str(tm)+' next_completion_time: '+str(next_completion_time))

    arr_cost=0
    dep_cost=0

    j=0

    while len(Ac_queue)>0:

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
                Ac_Infoi.status = FlightStatus.FINISHED # JF Note: check with Rob
                arr_cost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight)
                totserv+=1
            else: # JF Question: is this clause needed?
                Ac_Infoi.status = FlightStatus.FINISHED
                arr_cost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight)

            f.write(str(SubPolicy)+','+str(rep)+','+str(AC)+','+str(Ac_Infoi.flight_id)+','+str(prev_class)+','+str(current_class)+','+str(Time_Sep[prev_class][current_class]/60)+','+str(Ac_Infoi.orig_sched_time)+','+str(Ac_Infoi.ps_time)+','+str(Ac_Infoi.pool_time)+','+str(Ac_Infoi.release_time)+','+str(Ac_Infoi.travel_time)+','+str(Ac_Infoi.weather_state)+','+str(Ac_Infoi.enters_service)+','+str(Ac_Infoi.service_time)+','+str(Ac_Infoi.service_completion_time)+','+str(max(0,finish_time-(Ac_Infoi.ps_time+cost_fn.thres1)))+','+str(finish_time-(Ac_Infoi.pool_time+Ac_Infoi.travel_time))+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight))+',')
            f.write(str(Ac_Infoi.counter)+',')

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

    return arr_cost,dep_cost,totserv,prev_class,Ac_finished,next_completion_time
