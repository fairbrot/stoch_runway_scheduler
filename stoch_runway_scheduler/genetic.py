from typing import List, Tuple
import logging
import random
import math
import time
import numpy as np 

from .utils import weather, SequenceInfo, FlightInfo, Cost
from .gamma import Gamma_GetServ, Gamma_GetServ_Future, Gamma_Conditional_GetServ
from .simulate import simulate_weather, simulate_flight_times

# JF: this is the main sim heuristic
def Genetic(Ac_Info: List[FlightInfo], Arr_Pool, Arr_NotReady, Ac_queue, tm, k, prev_class, GA_PopList, GA_Info, GA_LoopSize, GA_CheckSize, GA_counter, basecost, wlb, wub, Opt_List, soln_evals_tot, soln_evals_num, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], cost_fn: Cost, GA_Check_Increment: int, Opt_Size: int, w_rho: float, wiener_sig: float, weather_sig: float):
    # JF Note: could maybe remove argument Max_LookAhead if no_ACs can be inferred from other arguments
    stepthrough_logger = logging.getLogger("stepthrough")
    step_summ_logger = logging.getLogger("step_summ")
    step_new_logger = logging.getLogger("step_new")

    output = 0 # output=1 means we're printing results as we go along; output=2 means we're outputting results to "Detailed" csv file
    pruned = 0 # indicator of whether or not the number of sequences has gone below the minimum number

    stepthrough_logger.info('Now entering Genetic procedure')

    NoA = len(Ac_Info)

    ArrTime, Trav_Time, ServTime = simulate_flight_times(tm, Ac_Info, tau, k, wiener_sig)
    
    # Before proceeding, randomly generate wlb_gen and wub_gen
    wlb_gen, wub_gen = simulate_weather(tm, wlb, wub, weather_sig)

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
            weather_state = weather(rel_time, wlb_gen, wub_gen) # weather(queue_complete, wlb_gen, wub_gen)
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
            queue_complete, straight_into_service = Gamma_GetServ(k, Time_Sep, rel_time, trav_time, perm_prev_class, cur_class, tm, weather_state, w_rho)
            perm_prev_class = cur_class
            basecost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, trav_time, queue_complete, Ac_Infoi.passenger_weight)


    else:

        queue_complete=tm
        perm_prev_class=prev_class

    stored_prev_class = perm_prev_class
    stored_queue_complete = queue_complete

    # Try all the sequences in the population

    for info in GA_Info:

        stepthrough_logger.info('Now trying sequence %s', info.sequence)
        stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost = basecost
        latest_tm = tm
        perm_prev_class = stored_prev_class
        perm_queue_complete = queue_complete

        perm = info.sequence
        info.n_traj += 1
        gam=1/info.n_traj #

        no_ACs = min(Max_LookAhead, len(perm))
        for index in range(no_ACs):

            AC = perm[index]
            Ac_Infoi = Ac_Info[AC]
            perm_class = Ac_Infoi.ac_class
            reltime = max(latest_tm,ArrTime[AC])
            begin_serv = max(reltime,perm_queue_complete)
            weather_state = weather(reltime,wlb_gen,wub_gen)


            stepthrough_logger.info(str(AC) + ',' + str(perm_class) + ',' + str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC]) + ',' + str(reltime) + ',' + str(Trav_Time[AC]) + ',' + str(perm_queue_complete)+',')

            AC_FinishTime, straight_into_service = Gamma_GetServ_Future(k, Time_Sep, reltime, ServTime[AC], Trav_Time[AC], perm_prev_class, perm_class, perm_queue_complete,weather_state, w_rho)

            info.queue_probs[index] = (1 - gam) * info.queue_probs[index] + gam*straight_into_service

            permcost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC], Trav_Time[AC], AC_FinishTime, Ac_Infoi.passenger_weight)
            latest_tm = reltime


            stepthrough_logger.info(str(AC_FinishTime)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight))+'\n')

            perm_queue_complete=AC_FinishTime
            perm_prev_class=perm_class

        stepthrough_logger.info('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(info.v)+',')

        info.v=(1-gam)*info.v+gam*permcost
        info.w += permcost**2
        stepthrough_logger.info('Total cost: '+','+str(info.v)+','+'Queue probs: '+','+str(info.queue_probs)+'\n'+'\n')

    GA_Info.sort(key=lambda x: x.sequence)
    for j in range(len(GA_Info)):
        step_summ_logger.info(str(GA_Info[j].v)+',')
    step_summ_logger.info('\n')

    ######### NOW UPDATE OPT LIST SEQS ###################

    for info in Opt_List:

        stepthrough_logger.info('\n'+'Now trying sequence '+','+str(info.sequence)+'\n')
        stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost=basecost
        latest_tm=tm
        perm_prev_class=stored_prev_class
        perm_queue_complete=queue_complete
        #perm_weather_state=weather_state

        info=info

        #print('info.sequence: '+str(info.sequence))
        perm=info.sequence
        info.n_traj+=1
        #gam=0.01
        #gam=2/(GA_Info[j].n_traj+1)
        gam=1/info.n_traj
        #gam=0.1

        #print('GA_Info: '+str(GA_Info))

        #print('Opt_List: '+str(Opt_List))
        no_ACs = min(Max_LookAhead,len(perm))
        for index in range(no_ACs):
            # index=perm[i]
            AC=perm[index]
            Ac_Infoi=Ac_Info[AC]
            perm_class=Ac_Infoi.ac_class
            reltime=max(latest_tm,ArrTime[AC])
            begin_serv=max(reltime,perm_queue_complete)
            weather_state=weather(reltime,wlb_gen,wub_gen) #weather(begin_serv,wlb_gen,wub_gen)

            stepthrough_logger.info(str(AC)+','+str(perm_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(Trav_Time[AC])+','+str(perm_queue_complete)+',')


            #perm_weather_state=weather(reltime,wlb,wub)
            AC_FinishTime, straight_into_service=Gamma_GetServ_Future(k, Time_Sep, reltime, ServTime[AC], Trav_Time[AC], perm_prev_class, perm_class, perm_queue_complete, weather_state, w_rho)

            info.queue_probs[index]=(1-gam)*info.queue_probs[index]+gam*straight_into_service

            permcost+=cost_fn(Ac_Infoi.orig_sched_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight) #Ac_Infoi.passenger_weight*(AC_FinishTime-(Ac_Infoi.ps_time+thres))**2
            latest_tm=reltime

            stepthrough_logger.info(str(AC_FinishTime)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight))+'\n')

            perm_queue_complete=AC_FinishTime
            perm_prev_class=perm_class

        stepthrough_logger.info('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j].v)+',')



        info.v=(1-gam)*info.v+gam*permcost
        info.w+=permcost**2 #Sum of squares

        stepthrough_logger.info('Total cost: '+','+str(info.v)+','+'Queue probs: '+','+str(info.queue_probs)+'\n'+'\n')

    Opt_List.sort(key=lambda x: x.sequence)
    for j in range(len(Opt_List)):
        step_summ_logger.info(str(info.v)+',')
    step_summ_logger.info('\n')

    # minperm=0
    # mincost=0
    # for j in range(len(GA_Info)):
    # 	if j==0 or GA_Info[j].v<mincost:
    # 		mincost=GA_Info[j].v
    # 		minperm=j

    GA_Info.sort(key=lambda x: x.v)
    Opt_List.sort(key=lambda x: x.v)

    GA_counter+=1

    if GA_counter>=GA_CheckSize:

        step_new_logger.info('GA_counter: '+','+str(GA_counter)+'\n')
        step_new_logger.info('Arr_Pool: '+','+str(Arr_Pool)+'\n')
        step_new_logger.info('Ac_queue: '+','+str(Ac_queue)+'\n'+'\n')
        step_new_logger.info('GA_PopList:'+'\n')
        for i in range(len(GA_Info)):
            step_new_logger.info(str(GA_Info[i])+'\n')
        step_new_logger.info('\n')
        step_new_logger.info('Opt_List:'+'\n')
        for i in range(len(Opt_List)):
            step_new_logger.info(str(Opt_List[i])+'\n')
        step_new_logger.info('\n')

        # if totserv>315:
        # 	print('tm: '+str(tm)+' GA_counter: '+str(GA_counter)+' len(GA_Info): '+str(len(GA_Info))+' len(Opt_List): '+str(Opt_List))

        t_val=1.96 #97.5th percentile of normal dist

        for info in GA_Info:

            info=GA_Info[j]

            if info.v>0:

                mn1=info.v
                n1=info.n_traj
                var1=(info.w-(n1*mn1**2))/(n1-1)

                for m in range(len(GA_Info)):

                    GA_Infom=GA_Info[m]

                    if GA_Infom.v>0:

                        mn2=GA_Infom.v
                        n2=GA_Infom.n_traj
                        var2=(GA_Infom.w-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:


                        # 	step_new_logger.info(str(info.sequence)+','+' loses to '+','+str(GA_Infom.sequence)+'\n')

                            info.v=-1
                            break

                        elif mn2>mn1+w_val:


                        # 	step_new_logger.info(str(GA_Infom.sequence)+','+' loses to '+','+str(info.sequence)+'\n')

                            GA_Infom.v=-1

            if info.v>0:

                for m in range(len(Opt_List)):

                    Opt_Listm=Opt_List[m]

                    if Opt_Listm.v>0:

                        mn2=Opt_Listm.v
                        n2=Opt_Listm.n_traj
                        var2=(Opt_Listm.w-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:

                        # 	step_new_logger.info(str(info.sequence)+','+' loses to '+','+str(Opt_Listm.sequence)+'\n')

                            info.v=-1
                            break

                        elif mn2>mn1+w_val:


                        # 	step_new_logger.info(str(Opt_Listm.sequence)+','+' loses to '+','+str(info.sequence)+'\n')

                            Opt_Listm.v=-1

        for j in range(len(Opt_List)):

            Opt_Listj=Opt_List[j]

            if Opt_Listj.v>0:

                mn1=Opt_Listj.v
                n1=Opt_Listj.n_traj
                var1=(Opt_Listj.w-(n1*mn1**2))/(n1-1)

                for m in range(len(GA_Info)):

                    GA_Infom=GA_Info[m]

                    if GA_Infom.v>0:

                        mn2=GA_Infom.v
                        n2=GA_Infom.n_traj
                        var2=(GA_Infom.w-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:

                        # 	step_new_logger.info(str(Opt_Listj.sequence)+','+' loses to '+','+str(GA_Infom.sequence)+'\n')

                            Opt_Listj.v=-1
                            break

                        elif mn2>mn1+w_val:

                        # 	step_new_logger.info(str(GA_Infom.sequence)+','+' loses to '+','+str(Opt_Listj.sequence)+'\n')

                            GA_Infom.v=-1

            if Opt_Listj.v>0:

                for m in range(len(Opt_List)):

                    Opt_Listm=Opt_List[m]

                    if Opt_Listm.v>0:
                        mn2=Opt_Listm.v
                        n2=Opt_Listm.n_traj
                        var2=(Opt_Listm.w-(n2*mn2**2))/(n2-1)
                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))
                        if mn1>mn2+w_val:
                            Opt_Listj.v=-1
                            break

                        elif mn2>mn1+w_val:
                            Opt_Listm.v=-1

        j=0
        while j<len(GA_Info):
            GA_Infoj=GA_Info[j]
            if GA_Infoj.v<0:
                soln_evals_tot+=GA_Infoj.n_traj
                soln_evals_num+=1
                GA_Info.remove(SequenceInfo(GA_Infoj.sequence,GA_Infoj.n_traj,GA_Infoj.v,GA_Infoj.queue_probs,GA_Infoj.w))
            else:
                step_new_logger.info('Retained sequence '+','+str(GA_Infoj)+','+' in GA_PopList'+'\n')
                j+=1

        j=0
        while j<len(Opt_List):
            Opt_Listj=Opt_List[j]
            if Opt_Listj.v<0:
                soln_evals_tot+=Opt_Listj.n_traj
                soln_evals_num+=1
                Opt_List.remove(SequenceInfo(Opt_Listj.sequence,Opt_Listj.n_traj,Opt_Listj.v,Opt_Listj.queue_probs,Opt_Listj.w))
            else:
                step_new_logger.info('Retained sequence '+','+str(Opt_Listj)+','+' in Opt_List'+'\n')
                j+=1

        if len(GA_Info)+len(Opt_List)<=Opt_Size:

            solns_left=len(GA_Info)+len(Opt_List)

            soln_evals_tot+=(solns_left*GA_LoopSize)
            soln_evals_num+=solns_left
            pruned=1

        if GA_counter>=GA_LoopSize:
            GA_CheckSize=GA_Check_Increment
            solns_left=len(GA_Info)+len(Opt_List)
            soln_evals_tot+=(solns_left*GA_LoopSize)
            soln_evals_num+=solns_left
        else:
            GA_CheckSize+=GA_Check_Increment
        step_new_logger.info('New GA_CheckSize: '+str(GA_CheckSize)+'\n'+'\n')

    Ac_added=[]

    if len(Arr_Pool)>0:

        if len(Opt_List)>0 and len(GA_Info)>0:
            if Opt_List[0].v < GA_Info[0].v:
                perm=Opt_List[0]
            else:
                perm=GA_Info[0]
        elif len(Opt_List)>0:
            perm=Opt_List[0]
        elif len(GA_Info)>0:
            perm=GA_Info[0]
        else:
            assert 1==2

        if perm.sequence[0] in Arr_Pool:

            counter=perm.n_traj
            qp=perm.queue_probs[0]

            if counter>=GA_Check_Increment or pruned==1:
                if qp>0: #0.05:
                    j=0
                    while perm.queue_probs[j] > 0: #0.05:
                        AC=perm.sequence[j]
                        Ac_added.append(AC)
                        step_new_logger.info('Counter is '+','+str(counter)+', ss_prob is '+','+str(perm.queue_probs[j])+', Adding AC '+','+str(AC)+'\n')
                        j+=1
                        if j==len(perm.sequence):
                            break

        else:
            counter=0
            qp=0

    else:
        counter=0
        qp=0


    return Ac_added, counter, qp, pruned, GA_CheckSize, GA_counter, soln_evals_tot, soln_evals_num
