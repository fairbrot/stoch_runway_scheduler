from typing import List, TextIO
import logging
import math
import random
from .utils import weather, getcost, FlightStatus, FlightInfo
from .gamma import sample_cond_gamma
from .annealing_cost import Annealing_Cost

def generate_weather(wlb: int, wub: int, T: int, weather_sig: float, freq: int):
    """
    Generates brownian motion for bad weather forecasts, as well as actual times of bad weather.

    Parameters
    ----------
    wlb: prediction at time 0 for start of bad weather
    wub: prediction at time 0 for end of bad weather
    T: length of time horizon in minutes (used to pre-allocate output arrays)
    weather_sig: standard deviation of Brownian motion
    freq: number of updates for each minute

    Returns
    -------
    wlb_tm: int
        actual start time of bad weather
    wub_tm:
        actual end time of bad weather
    weather_lb: List[float]
        forecast for start of bad weather over time. Each element i corresponds to forecast at time j*freq.
    weather_ub: List[float]
        forecast for end of bad weather over time.
    """
    N = T*freq # Size of weather array
    weather_lb = [0]*N # Dynamically forcast for start of bad weather T_0(t)
    weather_ub = [0]*N # Dynamically forcast for end of bad weather T_1(t)

    # called U_0 and U_1 in the paper
    wlb_tm = 0 # Actual (randomly generated) time at which bad weather starts; leave this as zero
    wub_tm = 0 # Actual (randomly generated) time at which bad weather ends; leave this as zero

    weather_lb[0] = wlb # Set initial forecasted values
    weather_ub[0] = wub
    old_lb = wlb
    old_ub = wub

    j = 0
    # Brownian motion again used to get dynamic forecast prediction
    for j in range(N):
        new_lb = random.gauss(old_lb,0.1*weather_sig) # random.gauss(old_lb,0.01) JF Question: why is 0.1 inside here - shouldn't this be incorporated into weather_sig
        new_ub = random.gauss(old_ub,0.1*weather_sig)
        if new_lb > new_ub: # JF Question: Why is this needed?
            new_ub = new_lb
        if j/freq >= new_lb and wlb_tm == 0:
            wlb_tm = j/freq
        if j/freq >= new_ub and wub_tm == 0:
            wub_tm = j/freq
        weather_lb[j] = new_lb
        weather_ub[j] = new_ub
        old_lb = new_lb
        old_ub = new_ub
    return wlb_tm, wub_tm, weather_lb, weather_ub

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
        ETA = tau
    else:
        chk = 0

    j = 0
    while True:
        j += 1 # step forward in increments of 1/freq minutes
        if j > Dep_time*freq: # only update ETA if we've gone beyond the AC's departure time ##- JF: I think logic is wrong here - j needs scaling (like below)
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

def Calculate_FCFS(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, wlb_tm, wub_tm, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float):

    tm=0

    totserv=0
    totcost=0
    prev_class=4
    queue_complete=0
    weather_state=0

    FCFS_Seq=[0]*NoA
    for i in range(NoA):
        FCFS_Seq[i]=ArrTime_Sorted[i][1]

    FCFS_cost = Annealing_Cost(FCFS_Seq, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, wlb_tm, wub_tm, 0, NoA, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)

    return FCFS_cost


def Posthoc_Check(seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,output, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float):

    perm=seq
    perm_cost=0
    latest_tm=0
    perm_prev_class=4
    perm_queue_complete=0
    perm_weather_state=0
    j=0

    while j<NoA:

        AC=perm[j]
        Ac_Infoi=Ac_Info[AC]
        release_time=Ac_Infoi.release_time #max(latest_tm,ArrTime[AC][0])
        #print('j: '+str(j)+' AC: '+str(AC)+' release_time: '+str(release_time))
        trav_time=Ac_Infoi.travel_time
        perm_class=Ac_Infoi.ac_class
        begin_serv=max(release_time,perm_queue_complete)
        perm_weather_state=weather(release_time,wlb_tm,wub_tm) #weather(begin_serv,wlb_tm,wub_tm)

        if perm_weather_state==1:
            ws=1/w_rho
        else:
            ws=1
        rate=ws*k/(Time_Sep[perm_prev_class][perm_class]/60)
        serv=ServTime[AC]/rate


        t1=release_time+trav_time
        t2=perm_queue_complete+serv
        perm_queue_complete=max(t1,t2)

        perm_cost+=getcost(Ac_Infoi.orig_sched_time,ArrTime[AC][0],trav_time,perm_queue_complete,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2)

        if output==1:
            print('AC: '+str(AC)+' class: '+str(perm_class)+' release_time: '+str(release_time)+' trav_time: '+str(trav_time)+' begin_serv: '+str(begin_serv)+' t1: '+str(t1)+' t2: '+str(t2)+' finish time: '+str(perm_queue_complete)+' weather state: '+str(perm_weather_state)+' pax weight: '+str(Ac_Infoi.passenger_weight)+' cost: '+str(getcost(Ac_Infoi.ps_time,ArrTime[AC][0],trav_time,perm_queue_complete,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2)))

        latest_tm=release_time
        perm_prev_class=perm_class

        j+=1

    return perm_cost

def Get_Actual_Serv(AC,prev_class,cur_class,weather_state,k, Time_Sep: List[List[int]], Ac_Info, w_rho: float):
    # Samping queue service time
    # Service time also depends on state of weather
    # Ac_Info[AC].service_rns is pre-generated random number which is then scaled appropriately

    serv_percs=Ac_Info[AC].service_rns
    if weather_state==1:
        ws=1/w_rho
    else:
        ws=1

    rate = ws*k/(Time_Sep[prev_class][cur_class]/60)

    servtime = serv_percs/rate #Transformation causes serv_percs to go from [mean k, var k] to [mean e_{ij}, var e_{ij}^2/k]

    return servtime

# JF: This isn't used anywhere currently - can it be deleted?
def GetServTime(trav_time,rel_time,prev_class,cur_class,tm,sv_time,ee,weather_state,gamma_cdf, k: int, Time_Sep: List[List[int]], w_rho: float):

    stepthrough_logger = logging.getLogger('stepthrough')
    #print('GetServTime')

    #This only for ACs that are already in the queue (either in service or not yet in service)

    t1=rel_time+trav_time

    rate=k/(Time_Sep[prev_class][cur_class]/60)
    if weather_state==1:
        rate*=1/w_rho #=0.5

    #t2=tm

    cond_serv=sample_cond_gamma(rate*(tm-sv_time),gamma_cdf) #rate*(tm_sv_time) gives the time that service has already been in progress after converting to the Gamma(k,1) scale
    cond_serv*=1/rate #Convert back to the correct scale

    t2=sv_time+cond_serv #t2+=rem_serv

    if ee==1:
        stepthrough_logger.info(str(t2-tm)+',')

    # if len(Ac_queue)==2 and ee==1:
    # 	print('Trav time completed at '+str(t1)+', serv time completed at '+str(t2))

    if t1<t2:
        straight_into_service=0
        t_out=t2
    else:
        straight_into_service=1
        t_out=t1

    #print('Out of GetServTime')

    return t_out,straight_into_service

# JF: again this isn't used anywhere
def GetServTime_Future(trav_time,serv_time,rel_time,prev_class,cur_class,tm,ee,weather_state, k: int, Time_Sep: List[List[int]], w_rho: float): 

    stepthrough_logger = logging.getLogger('stepthrough')
    #This is for ACs that have not yet been added to the queue

    #print('tm: '+str(tm))

    t1=rel_time+trav_time

    #Next, get s1+Z2
    # if type(prev_class)==list:
    # 	print('prev_class: '+str(prev_class))
    # if type(cur_class)==list:
    # 	print('cur_class: '+str(cur_class))

    rate=k/(Time_Sep[prev_class][cur_class]/60)
    if weather_state==1:
        rate*=1/w_rho #=0.5

    t2=tm

    sv_time=serv_time/rate #Convert from the Gamma(k,1) scale to the correct scale

    t2+=sv_time

    # for m in range(ph_B):
    # 	# if m>len(serv_time)-1:
    # 	# 	print('ph_B: '+str(ph_B))
    # 	t2+=(-1/rate)*math.log(serv_time[m])

    if ee==1:
        stepthrough_logger.info(str(t2-tm)+',')

    # if len(Ac_queue)==2 and ee==1:
    # 	print('Trav time completed at '+str(t1)+', serv time completed at '+str(t2))

    if t1<t2:
        straight_into_service=0
        t_out=t2
    else:
        straight_into_service=1
        t_out=t1

    return t_out,straight_into_service

def round_down(tm: float, freq: int) -> float:
    """Rounds tm down to beginning of interval defined by freq.

    For example, if freq is 100, then we round down for step size of 0.01,
    that is 0.32453 would become 0.32.
    """
    int_size = 1/freq
    i = tm // int_size
    return i*int_size


def Update_ETAs(Ac_Info: List[FlightInfo], Arr_NotReady: List[int], Dep_NotReady: List[int], 
                Ac_queue: List[int], tm: float, Brown_Motion: List[List[float]], Arr_Pool: List[int], tau: float,
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
        if Ac_Infoi.travel_time_indicator == 0:
            rel_time = Ac_Infoi.release_time
            trav_so_far = tm - rel_time # amount of time spent travelling to the runway so far
            # JF Question: why round? Should this be round down? Is this related to freq?
            # I think we need rounding down to the beginning of interval defined by freq
            rounded_trav_so_far = round_down(trav_so_far, freq)
            if rounded_trav_so_far >= Ac_Infoi.travel_time:
                msg = '* Aircraft '+str(AC)+' has finished its travel time at '+str(tm)+'; still waiting for service time.'+'\n'+'\n'
                stepthrough_logger.info(msg)
                step_summ_logger.info(msg)
                Ac_Infoi.travel_time_indicator = 1
            else:
                Ac_Infoi.eta = Brown_Motion[AC][int((Ac_Infoi.pool_time + rounded_trav_so_far) * freq)]

def Update_Stats(tm,AC,Ac_Info,Ac_queue,real_queue_complete,wlb_tm,wub_tm,latest_class,Ov_GA_counter,next_completion_time, k: int, Time_Sep: List[List[int]], w_rho: float, SubPolicy, counter, qp):
    # Function which produces statistics about a flight which has just been released - also calculates
    # when this flight will be finished being serviced

    Ac_Infoi=Ac_Info[AC]

    #Ac_Infoi.status += 1 # update status
    Ac_Infoi.status = FlightStatus.IN_QUEUE
    release_time = tm # release time
    begin_serv = max(release_time,real_queue_complete) # time that service begins
    cur_class = Ac_Infoi.ac_class

    get_weather_state = weather(release_time, wlb_tm, wub_tm) #weather(begin_serv,wlb_tm,wub_tm)

    actual_serv=Get_Actual_Serv(AC,latest_class,cur_class,get_weather_state,k,Time_Sep, Ac_Info, w_rho)
    trav_time=Ac_Infoi.travel_time

    t1=release_time+trav_time
    t2=begin_serv+actual_serv

    finish_time=max(t1,t2)
    real_queue_complete=finish_time

    Ac_Infoi.release_time=release_time
    Ac_Infoi.enters_service=begin_serv
    Ac_Infoi.weather_state=get_weather_state
    Ac_Infoi.service_time=actual_serv
    Ac_Infoi.service_completion_time=finish_time

    if SubPolicy in ('GA','GAD','VNS','VNSD'):
        Ac_Infoi.counter=Ov_GA_counter
        Ov_GA_counter=0
    else:
        Ac_Infoi.counter=counter
    Ac_Infoi.qp=qp

    latest_class=cur_class

    if len(Ac_queue)==1:
        next_completion_time=finish_time

    return real_queue_complete,next_completion_time,latest_class,Ov_GA_counter

def Serv_Completions(Ac_Info, Ac_queue, prev_class, totserv, Ac_finished, tm, next_completion_time, thres1: int, thres2: int, lam1: float, lam2: float, f: TextIO, SubPolicy, rep, Time_Sep: List[List[int]], Left_queue):

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
                Ac_Infoi.status = FlightStatus.DEP_NOT_READY # JF Note: check with Rob
                arr_cost += getcost(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2)
                #print('* Cost incurred is '+str(arr_cost))
                totserv+=1
            else: # JF Question: is this clause needed?
                Ac_Infoi.status = FlightStatus.FINISHED
                arr_cost += getcost(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight, thres1, thres2, lam1, lam2)

            f.write(str(SubPolicy)+','+str(rep)+','+str(AC)+','+str(Ac_Infoi.flight_id)+','+str(prev_class)+','+str(current_class)+','+str(Time_Sep[prev_class][current_class]/60)+','+str(Ac_Infoi.orig_sched_time)+','+str(Ac_Infoi.ps_time)+','+str(Ac_Infoi.pool_time)+','+str(Ac_Infoi.release_time)+','+str(Ac_Infoi.travel_time)+','+str(Ac_Infoi.weather_state)+','+str(Ac_Infoi.enters_service)+','+str(Ac_Infoi.service_time)+','+str(Ac_Infoi.service_completion_time)+','+str(max(0,finish_time-(Ac_Infoi.ps_time+thres1)))+','+str(finish_time-(Ac_Infoi.pool_time+Ac_Infoi.travel_time))+','+str(Ac_Infoi.passenger_weight)+','+str(getcost(Ac_Infoi.ps_time,Ac_Infoi.pool_time,Ac_Infoi.travel_time,finish_time,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2))+',')
            f.write(str(Ac_Infoi.counter)+','+str(Ac_Infoi.qp)+',')

            if SubPolicy in ('GA','GAD','VNS','VNSD'):
                f.write(str(Ac_Infoi.pred_cost)+',')
            # if SubPolicy in ('VNS'):
            # 	f.write(str(Loop_Evals/Loop_Nums)+',')
            # elif SubPolicy in ('VNSD'):
            # 	f.write(',')
            # if SubPolicy in ('GA','GAD','VNS','VNSD'):
            # 	f.write(str(elap_tot/elap_num)+','+str(Repop_elap_tot/Repop_elap_num)+','+str(Pop_elap_tot/Pop_elap_num)+',')
                #f.write(str(elap)+','+str(Repop_elap)+','+str(Pop_elap)+',')
            f.write('\n')

            prev_class = current_class

            Left_queue.append(AC)
            Ac_queue.remove(AC)

            #print('Arr_NotReady is: '+str(Arr_NotReady))
            #print('Arr_Pool is: '+str(Arr_Pool))
            #print('Ac_queue is: '+str(Ac_queue))
            #print('Left_queue is: '+str(Left_queue))
            print(str(SubPolicy)+' totserv: '+str(totserv))

            if len(Ac_queue)>0:
                New_AC=Ac_queue[0]
                next_completion_time=Ac_Info[New_AC][16]

        else:
            break

    return arr_cost,dep_cost,totserv,prev_class,Ac_finished,next_completion_time
