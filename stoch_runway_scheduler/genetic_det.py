from typing import List
import logging
import time
from .utils import FlightInfo, Cost
from .weather import WeatherProcess, WeatherStatus
from .sequence import SequenceInfo
from .gamma import gamma_cond_exp

# JF: This is the main deterministic heuristic
# Name may not be best choice
def Genetic_determ(Ac_Info: List[FlightInfo], Arr_Pool: List[int], Arr_NotReady: List[int], 
                    Ac_queue: List[int], tm: float, NoA: int, k:int, 
                    prev_class: int, GA_Info: List[SequenceInfo], weather: WeatherProcess, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], cost_fn: Cost, 
                    tot_arr_cost: float, tot_dep_cost: float, w_rho: float):

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

    basecost = tot_arr_cost + tot_dep_cost

    stepthrough_logger.info('basecost is %.2f', basecost)
    stepthrough_logger.info('Generated results for ACs already in the queue are as follows:')
    stepthrough_logger.info('AC, Class, Time Sep, Release time, Travel time, Enters serv, Actual serv, Finish time, Pax weight, Cost')

    if len(Ac_queue) > 0:

        # Need to generate service times for AC already in the queue; first consider the customer in position 0
        AC = Ac_queue[0]
        Ac_Infoi = Ac_Info[AC]
        rel_time = Ac_Infoi.release_time
        sv_time = Ac_Infoi.enters_service
        cur_class = Ac_Infoi.ac_class

        t1=Ac_Infoi.eta

        # Get the conditional expectation of service time based on service time elapsed so far

        if Ac_Infoi.weather_state==1:
            beta = k/(w_rho*Time_Sep[prev_class][cur_class]/60)
        else:
            beta = k/(Time_Sep[prev_class][cur_class]/60)
        alpha = k
        cond_tm = gamma_cond_exp(tm - sv_time, alpha, beta)

        t2 = sv_time + cond_tm
        #print('Time remaining based on mean value: '+str(sv_time+cond_tm-tm))

        queue_complete = max(t1,t2)

        basecost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, tau, queue_complete, Ac_Infoi.passenger_weight)
        stepthrough_logger.info('%.2f, %.2f, %.2f', queue_complete, Ac_Infoi.passenger_weight, 
                                cost_fn(Ac_Infoi.ps_time, Ac_Infoi.pool_time, tau, queue_complete, Ac_Infoi.passenger_weight))

        perm_prev_class=cur_class


        # Now consider the rest of the customers in the queue
        for AC in Ac_queue:
            Ac_Infoi=Ac_Info[AC]
            rel_time=Ac_Infoi.release_time
            cur_class=Ac_Infoi.ac_class

            stepthrough_logger.info(str(AC)+','+str(Ac_Infoi.ac_class)+','+str(Time_Sep[perm_prev_class][cur_class]/60)+','+str(Ac_Infoi.release_time)+','+str(Ac_Infoi.travel_time)+','+str(Ac_Infoi.enters_service)+',')

            t1 = Ac_Infoi.eta
            weather_state = weather(rel_time) # weather(queue_complete)
            if weather_state == 1:
                t2 = queue_complete + (w_rho * Time_Sep[perm_prev_class][cur_class]/60)
                stepthrough_logger.info(str(w_rho*Time_Sep[perm_prev_class][cur_class]/60)+',')
            else:
                t2=queue_complete+(Time_Sep[perm_prev_class][cur_class]/60)
                stepthrough_logger.info(str(Time_Sep[perm_prev_class][cur_class]/60)+',')
            queue_complete=max(t1,t2)

            perm_prev_class=cur_class
            basecost+=cost_fn(Ac_Infoi.orig_sched_time,Ac_Infoi.pool_time,tau,queue_complete,Ac_Infoi.passenger_weight)

            stepthrough_logger.info(str(queue_complete)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time,Ac_Infoi.pool_time,tau,queue_complete,Ac_Infoi.passenger_weight))+'\n')

    else:

        queue_complete=tm
        perm_prev_class=prev_class

    stored_prev_class=perm_prev_class
    stored_queue_complete=queue_complete

    #Try all the sequences in the population

    for info in GA_Info:

        stepthrough_logger.info('\n'+'Now trying sequence '+','+str(info.sequence)+'\n')
        stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost = basecost
        latest_tm = tm
        perm_prev_class = stored_prev_class
        perm_queue_complete = queue_complete

        perm = info.sequence
        info.n_traj += 1
        gam = 1/info.n_traj

        # JF Question: why would this not be len(perm)?
        no_ACs = min(Max_LookAhead, len(perm)) #min(AC_List_Length,len(Arr_Pool)+len(Arr_NotReady))
        for index in range(no_ACs):

            #index=perm[i]
            AC = perm[index]
            Ac_Infoi = Ac_Info[AC]
            perm_class = Ac_Infoi.ac_class
            reltime = max(latest_tm,ArrTime[AC])
            stepthrough_logger.info(str(AC)+','+str(Ac_Infoi.ac_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(tau)+','+str(perm_queue_complete)+',')

            t1=reltime+tau
            # JF Question - should we use time start servicing here?
            # Like begin_serv in Genetic?
            if weather(reltime) == WeatherStatus.BAD:
                exp_serv = w_rho*Time_Sep[perm_prev_class][perm_class]/60
            else:
                exp_serv = Time_Sep[perm_prev_class][perm_class]/60
            t2 = perm_queue_complete + exp_serv

            stepthrough_logger.info(str(exp_serv)+',')

            AC_FinishTime = max(t1,t2)
            if t1>=t2:
                straight_into_service=1
            else:
                straight_into_service=0

            info.queue_probs[index] = (1-gam)*info.queue_probs[index] + gam*straight_into_service

            permcost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC], tau, AC_FinishTime, Ac_Infoi.passenger_weight)
            latest_tm = reltime

            stepthrough_logger.info(str(AC_FinishTime) + ',' + str(Ac_Infoi.passenger_weight) + ',' + str(cost_fn(Ac_Infoi.ps_time, ArrTime[AC], tau, AC_FinishTime, Ac_Infoi.passenger_weight))+'\n')

            perm_queue_complete = AC_FinishTime
            perm_prev_class = perm_class


        stepthrough_logger.info('Final cost: ' + ',' +str(permcost) + ',' + 'Gamma: ' + ',' + str(gam) + ',' + 'Old cost: ' + ','+str(info.v) + ',')

        info.v = (1-gam)*info.v + gam*permcost

        stepthrough_logger.info('Total cost: '+','+str(info.v)+','+'Queue probs: '+','+str(info.queue_probs)+'\n'+'\n')

    GA_Info.sort(key=lambda x: x.sequence)
    for j in range(len(GA_Info)):
        step_summ_logger.info(str(info.v)+',')

    GA_Info.sort(key=lambda x: x.v)

    Ac_added=[]
    counter=0

    if len(Arr_Pool)>0:
        assert len(GA_Info) > 0
        perm = GA_Info[0]

        if perm.sequence[0] in Arr_Pool:

            if 1==1: #perm.n_traj>=100 and 1==1: #perm.queue_probs[0]>0: #0.05:
                j=0
                counter=perm.n_traj
                while 1==1: #perm[3][j]>0: #0.05:
                    AC=perm.sequence[j]
                    Ac_added.append(AC)
                    j+=1
                    if j==len(perm.sequence):
                        break

    # end_time=time.time()
    # elap=(end_time-start_time)/conv_factor
    # #elap=0.01

    #elap=fixed_elap_vnsd/conv_factor

    return Ac_added, counter, stored_queue_complete
