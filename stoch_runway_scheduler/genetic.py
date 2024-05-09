from typing import List
import logging
import numpy as np 

from .utils import FlightInfo, Cost
from .weather import StochasticWeatherProcess
from .sequence import SequenceInfo
from .gamma import Gamma_GetServ, Gamma_Conditional_GetServ
from .simulate import simulate_flight_times

# JF: this is the main sim heuristic - it is doing too much
def Genetic(Ac_Info: List[FlightInfo], Arr_Pool, Ac_queue, tm, k, prev_class, GA_Info, GA_LoopSize, GA_CheckSize, GA_counter, basecost, weather: StochasticWeatherProcess, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], cost_fn: Cost, GA_Check_Increment: int, S_min: int, w_rho: float, wiener_sig: float):
    # JF Note: could maybe remove argument Max_LookAhead if no_ACs can be inferred from other arguments
    
    # JF Question: is it an issue that Arr_NotReady is not an argument here? Sequences could possibly contain flights in this set
    stepthrough_logger = logging.getLogger("stepthrough")
    step_summ_logger = logging.getLogger("step_summ")
    step_new_logger = logging.getLogger("step_new")

    stepthrough_logger.info('Now entering Genetic procedure')

    ArrTime, Trav_Time, ServTime = simulate_flight_times(tm, Ac_Info, tau, k, wiener_sig)
    
    # Before proceeding, randomly generate wlb_gen and wub_gen
    weather_sample = weather.sample_process(tm)

    stepthrough_logger.info("basecost is %s\n", basecost)
    stepthrough_logger.info('Generated results for ACs already in the queue are as follows:')
    stepthrough_logger.info('AC, Class, Time Sep , Release time, Travel time, Enters serv, Actual serv, Finish time, Pax weight, Cost')

    if len(Ac_queue) > 0:

        # Need to generate service times for AC already in the queue; first consider the customer in position 0
        AC = Ac_queue[0]
        Ac_Infoi = Ac_Info[AC]
        rel_time = Ac_Infoi.release_time
        sv_time = Ac_Infoi.enters_service
        cur_class = Ac_Infoi.ac_class
        weather_state = Ac_Infoi.weather_state

        queue_complete, straight_into_service = Gamma_Conditional_GetServ(k, Time_Sep, Trav_Time[AC], rel_time, sv_time, prev_class, cur_class, tm, weather_state, w_rho)
        basecost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC], Trav_Time[AC], queue_complete, Ac_Infoi.passenger_weight)
        perm_prev_class = cur_class


        # Now consider the rest of the customers in the queue
        for j in range(1,len(Ac_queue)):
            AC = Ac_queue[j]
            Ac_Infoi = Ac_Info[AC]
            rel_time = Ac_Infoi.release_time
            cur_class = Ac_Infoi.ac_class
            # ROB TO CHECK - should we use max(queue_complete, rel_time) in line below?
            weather_state = weather_sample(rel_time)
            trav_time = Trav_Time[AC]

            # JF Question: Not sure how this first case arises - perhaps it is a mistake? trav_time isn't even defined for flight j
            if trav_time <= 0:
                trav_time=0
            else:
                # Why is this sampled again? Don't we already have travel time for flights in queue?
                trav_time = np.random.wald(tau, (tau/wiener_sig)**2)

            # if tm>=Ac_Infoi.eta: #this block of code is probably needed but wasn't included in the 5000 experiments for the paper
            # 	trav_time=Ac_Infoi.travel_time #travel time has already finished
            # else:
            # 	z=int(random.randrange(1,999))
            # 	sched=int(10*round(Ac_Infoi.eta-tm,1))
            # 	trav_time=wiener_cdf[sched][z]

            # JF: does tm need updating here for other flights
            queue_complete, straight_into_service = Gamma_GetServ(k, Time_Sep, rel_time, trav_time, perm_prev_class, cur_class, queue_complete, weather_state, w_rho)
            perm_prev_class = cur_class
            basecost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, trav_time, queue_complete, Ac_Infoi.passenger_weight)

    else:
        queue_complete=tm
        perm_prev_class=prev_class

    stored_prev_class = perm_prev_class

    # Try all the sequences in the population
    for info in GA_Info:
        stepthrough_logger.info('Now trying sequence %s', info.sequence)
        stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost = basecost
        latest_tm = tm
        perm_prev_class = stored_prev_class
        perm_queue_complete = queue_complete

        perm = info.sequence
        no_ACs = min(Max_LookAhead, len(perm)) # JF Question: can this not just be len(perm)?
        xi_list = [] # for each flight indicates whether or not it goes straight to service - called xi in paper
        for index in range(no_ACs):
            AC = perm[index]
            Ac_Infoi = Ac_Info[AC]
            perm_class = Ac_Infoi.ac_class
            reltime = max(latest_tm, ArrTime[AC])
            # JF Question: should we use begin_serv rather than reltime below?
            # begin_serv = max(reltime, perm_queue_complete)
            weather_state = weather_sample(reltime)

            stepthrough_logger.info('%d, %s, %.2f, %.2f, %.2f, %.2f, %.2f,',
                                    AC, perm_class, Time_Sep[perm_prev_class][perm_class]/60, ArrTime[AC], reltime, Trav_Time[AC], perm_queue_complete)

            AC_FinishTime, straight_into_service = Gamma_GetServ(k, Time_Sep, reltime, Trav_Time[AC], perm_prev_class, perm_class, 
                                                                perm_queue_complete,weather_state, w_rho, serv_time = ServTime[AC])
            xi_list.append(straight_into_service)

            permcost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC], Trav_Time[AC], AC_FinishTime, Ac_Infoi.passenger_weight)
            latest_tm = reltime

            # JF Question: cost_fn below not consistent with that above. Why?
            stepthrough_logger.info(str(AC_FinishTime)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight))+'\n')

            perm_queue_complete = AC_FinishTime
            perm_prev_class = perm_class

        stepthrough_logger.info('Final cost: '+','+str(permcost)+','+'N traj: '+','+str(info.n_traj)+','+'Old cost: '+','+str(info.v)+',')
        info.add_observation(permcost, xi_list)
        stepthrough_logger.info('Total cost: '+','+str(info.v)+','+'Queue probs: '+','+str(info.queue_probs)+'\n'+'\n')

    GA_Info.sort(key=lambda x: x.sequence)
    for info in GA_Info:
        step_summ_logger.info('%.2f, ', info.v)
    step_summ_logger.info('\n')

    GA_Info.sort(key=lambda x: x.v)

    GA_counter += 1

    if GA_counter >= GA_CheckSize:

        step_new_logger.info('GA_counter: %s, ', GA_counter)
        step_new_logger.info('Arr_Pool: %s, ', Arr_Pool)
        step_new_logger.info('Ac_queue: %s, ', Ac_queue)
        step_new_logger.info('GA_Info:'+'\n')
        for info in GA_Info:
            step_new_logger.info('%s', info)

        to_remove = SequenceInfo.rank_and_select(GA_Info)
        to_remove.sort(reverse=True)
        
        for i in to_remove:
            info = GA_Info.pop(i)

        # When iterations reaches GA_LoopSize (500) we reset GA_CheckSize - otherwise
        # we increase so ranking 
        if GA_counter >= GA_LoopSize:
            GA_CheckSize = GA_Check_Increment
        else:
            GA_CheckSize += GA_Check_Increment
        step_new_logger.info('New GA_CheckSize: '+str(GA_CheckSize)+'\n'+'\n')

    Ac_added=[]

    # Mark some flights in best sequence for release? Yes
    if len(Arr_Pool)>0:
        assert len(GA_Info) > 0
        perm = GA_Info[0] # JF Question: is this assuming the list is in a particular order? Yes
        if perm.sequence[0] in Arr_Pool:
            counter = perm.n_traj
            if counter >= GA_Check_Increment or (len(GA_Info) <= S_min):
                for j, AC in enumerate(perm.sequence):
                    if perm.queue_probs[j] <= 0:
                        break
                    Ac_added.append(AC)
                    step_new_logger.info('Counter is '+','+str(counter)+', ss_prob is '+','+str(perm.queue_probs[j])+', Adding AC '+','+str(AC)+'\n')

        else:
            counter=0

    else:
        counter=0

    return Ac_added, counter, GA_CheckSize, GA_counter
