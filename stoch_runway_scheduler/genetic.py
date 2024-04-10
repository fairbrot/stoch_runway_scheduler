from typing import List
import logging
import random
import math
import time
import numpy as np 

from .utils import weather, getcost
from .gamma import Gamma_GetServ, Gamma_GetServ_Future, Gamma_Conditional_GetServ

# JF: this is the main sim heuristic
def Genetic(Ac_Info,Arr_Pool,Arr_NotReady,Ac_queue,Left_queue,tm,NoA,k,prev_class,GA_PopList,GA_Info,GA_LoopSize,GA_CheckSize,GA_counter,basecost,wlb,wub, Opt_List, max_d, soln_evals_tot, soln_evals_num, gamma_cdf, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float, GA_Check_Increment: int, Opt_Size: int, w_rho: float, stepthrough: int, wiener_sig: float, weather_sig: float):
    # if len(Arr_Pool)+len(Ac_queue)>0:
    # 	print('Genetic')
    stepthrough_logger = logging.getLogger("stepthrough")
    step_summ_logger = logging.getLogger("step_summ")
    step_new_logger = logging.getLogger("step_new")

    output=0 #output=1 means we're printing results as we go along; output=2 means we're outputting results to "Detailed" csv file
    ee=0
    pruned=0 #indicator of whether or not the number of sequences has gone below the minimum number

    if stepthrough==1:
        ee=1
    start_time=time.time()

    stepthrough_logger.info('Now entering Genetic procedure')

    ArrTime=[0]*NoA
    ServTime=[0]*NoA
    Trav_Time=[0]*NoA

    #wiener_cdf_tau=wiener_cdf[int(10*tau)]

    #Generate arrival and service time percentiles for AC not yet in queue

    for AC in Arr_Pool:

        ArrTime[AC]=max(0,Ac_Info[AC].eta-tau)

        Trav_Time[AC]=np.random.wald(tau,(tau/wiener_sig)**2)

        ServTime[AC]=np.random.gamma(k,1)


    for AC in Arr_NotReady:

        sched=int(round(Ac_Info[AC].eta-(tm+tau),1))
        
        if sched<=0:
            ArrTime[AC]=tm
        else:
            ArrTime[AC]=np.random.wald(sched,(sched/wiener_sig)**2)+tm

        #st4.write(str(tm)+','+str(Ac_Info[AC].eta-(tm+tau))+','+str(ArrTime[AC])+',')

        Trav_Time[AC]=np.random.wald(tau,(tau/wiener_sig)**2)
        
        #st4.write(str(Trav_Time[AC])+'\n')

        #assert 1==2 


        ServTime[AC]=np.random.gamma(k,1)



    #Before proceeding, randomly generate wlb_gen and wub_gen
    chk=0
    while chk==0:
        if tm>=wlb:
            wlb_gen=wlb
        else:
            #Do wlb_gen
            sched=int(round(wlb-tm,1))
            # wald1.write(str(wlb_gen)+'\n')
            if sched<=0:
                wlb_gen=tm
            else:
                wlb_gen=np.random.wald(sched,(sched/weather_sig)**2)+tm
        if tm>=wub:
            wub_gen=wub
        else:
            #Do wub_gen
            sched=int(round(wub-tm,1))
            if sched<=0:
                wub_gen=tm
            else:
                wub_gen=np.random.wald(sched,(sched/weather_sig)**2)+tm
        if wlb_gen<=wub_gen:
            chk=1
        else:
            chk=1
            #print('tm: '+str(tm)+' wlb: '+str(wlb)+' wub: '+str(wub)+' wlb_gen: '+str(wlb_gen)+' wub_gen: '+str(wub_gen))

    stepthrough_logger.info("basecost is %s\n", basecost)
    stepthrough_logger.info('Generated results for ACs already in the queue are as follows:')
    stepthrough_logger.info('AC, Class, Time Sep , Release time, Travel time, Enters serv, Actual serv, Finish time, Pax weight, Cost')

    if len(Ac_queue)>0:

        #Need to generate service times for AC already in the queue; first consider the customer in position 0



        AC=Ac_queue[0]
        Ac_Infoi=Ac_Info[AC]
        rel_time=Ac_Infoi.release_time
        sv_time=Ac_Infoi.enters_service
        cur_class=Ac_Infoi.ac_class
        weather_state=Ac_Infoi.weather_state

        if tm>=Ac_Infoi.eta:
            trav_time=Ac_Infoi.travel_time #travel time has already finished
        else:
            sched=int(round(Ac_Infoi.eta-tm,1))
            if sched<=0:
                trav_time=0
            else:
                trav_time=np.random.wald(sched,(sched/wiener_sig)**2)

        queue_complete,straight_into_service=Gamma_Conditional_GetServ(k, Time_Sep, trav_time,rel_time,sv_time,prev_class,cur_class,tm,weather_state, gamma_cdf, w_rho)
        basecost+=getcost(Ac_Infoi.orig_sched_time,Ac_Infoi.pool_time,trav_time,queue_complete,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2)
        perm_prev_class=cur_class



        #Now consider the rest of the customers in the queue
        for j in range(1,len(Ac_queue)):

            AC=Ac_queue[j]
            Ac_Infoi=Ac_Info[AC]
            rel_time=Ac_Infoi.release_time
            cur_class=Ac_Infoi.ac_class
            weather_state=weather(rel_time,wlb_gen,wub_gen) #weather(queue_complete,wlb_gen,wub_gen)

            if trav_time<=0:
                trav_time=0
            else:
                trav_time=np.random.wald(tau,(tau/wiener_sig)**2)

            # if tm>=Ac_Infoi.eta: #this block of code is probably needed but wasn't included in the 5000 experiments for the paper
            # 	trav_time=Ac_Infoi.travel_time #travel time has already finished
            # else:
            # 	z=int(random.randrange(1,999))
            # 	sched=int(10*round(Ac_Infoi.eta-tm,1))
            # 	trav_time=wiener_cdf[sched][z]

            queue_complete,straight_into_service=Gamma_GetServ(k, Time_Sep, rel_time,trav_time,perm_prev_class,cur_class, tm, weather_state, gamma_cdf, w_rho)
            perm_prev_class=cur_class
            basecost+=getcost(Ac_Infoi.orig_sched_time,Ac_Infoi.pool_time,trav_time,queue_complete,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2)


    else:

        queue_complete=tm
        perm_prev_class=prev_class

    stored_prev_class=perm_prev_class
    stored_queue_complete=queue_complete

    #Try all the sequences in the population

    for j in range(len(GA_Info)):

        stepthrough_logger.info('Now trying sequence %s', GA_Info[j][0])
        stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost=basecost
        latest_tm=tm
        perm_prev_class=stored_prev_class
        perm_queue_complete=queue_complete
        #perm_weather_state=weather_state

        GA_Infoj=GA_Info[j]

        perm=GA_Infoj[0]
        GA_Infoj[1]+=1
        #gam=0.01
        #gam=2/(GA_Info[j][1]+1)
        gam=1/GA_Infoj[1]#
        #gam=0.1

        #print('GA_Info: '+str(GA_Info))

        no_ACs=min(Max_LookAhead,len(perm))
        for index in range(no_ACs):

            #index=perm[i]
            AC=perm[index]
            Ac_Infoi=Ac_Info[AC]
            perm_class=Ac_Infoi.ac_class
            reltime=max(latest_tm,ArrTime[AC])
            begin_serv=max(reltime,perm_queue_complete)
            weather_state=weather(reltime,wlb_gen,wub_gen) #weather(begin_serv,wlb_gen,wub_gen)

            if output==1:
                ee=1


            stepthrough_logger.info(str(AC)+','+str(perm_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(Trav_Time[AC])+','+str(perm_queue_complete)+',')

            AC_FinishTime, straight_into_service=Gamma_GetServ_Future(k, Time_Sep, reltime,ServTime[AC],Trav_Time[AC],perm_prev_class,perm_class,perm_queue_complete,weather_state, w_rho)

            GA_Infoj[3][index]=(1-gam)*GA_Infoj[3][index]+gam*straight_into_service

            permcost+=getcost(Ac_Infoi.orig_sched_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2) #Ac_Infoi.passenger_weight*(AC_FinishTime-(Ac_Infoi.ps_time+thres))**2
            latest_tm=reltime


            stepthrough_logger.info(str(AC_FinishTime)+','+str(Ac_Infoi.passenger_weight)+','+str(getcost(Ac_Infoi.ps_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2))+'\n')

            perm_queue_complete=AC_FinishTime
            perm_prev_class=perm_class

        stepthrough_logger.info('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j][2])+',')

        # if GA_Infoj[0]==[0,1,3,6,8,7,10,12,2,4,9,11,5,13,14,15,20,16,17,21,22,24,28,23,25,26,27,29,18,19]:
        # 	print('Evaluated cost for sequence '+str(GA_Infoj[0])+' as '+str(permcost))

        # if GA_Infoj[2]>0:
        # 	pct_diff=abs(permcost/GA_Infoj[2]-1)
        # 	if pct_diff>iter_max_d:
        # 		iter_max_d=pct_diff
        # else:
        # 	iter_max_d=1

        GA_Infoj[2]=(1-gam)*GA_Infoj[2]+gam*permcost
        GA_Infoj[4]+=permcost**2
        stepthrough_logger.info('Total cost: '+','+str(GA_Info[j][2])+','+'Queue probs: '+','+str(GA_Info[j][3])+'\n'+'\n')

    GA_Info.sort(key=lambda x: x[0])
    for j in range(len(GA_Info)):
        step_summ_logger.info(str(GA_Info[j][2])+',')
    step_summ_logger.info('\n')

    ######### NOW UPDATE OPT LIST SEQS ###################

    for j in range(len(Opt_List)):

        stepthrough_logger.info('\n'+'Now trying sequence '+','+str(Opt_List[j][0])+'\n')
        stepthrough_logger.info('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost=basecost
        latest_tm=tm
        perm_prev_class=stored_prev_class
        perm_queue_complete=queue_complete
        #perm_weather_state=weather_state

        Opt_Listj=Opt_List[j]

        #print('Opt_Listj[0]: '+str(Opt_Listj[0]))
        perm=Opt_Listj[0]
        Opt_Listj[1]+=1
        #gam=0.01
        #gam=2/(GA_Info[j][1]+1)
        gam=1/Opt_Listj[1]
        #gam=0.1

        #print('GA_Info: '+str(GA_Info))

        #print('Opt_List: '+str(Opt_List))
        no_ACs=min(Max_LookAhead,len(perm))
        for index in range(no_ACs):

            #index=perm[i]
            AC=perm[index]
            Ac_Infoi=Ac_Info[AC]
            perm_class=Ac_Infoi.ac_class
            reltime=max(latest_tm,ArrTime[AC])
            begin_serv=max(reltime,perm_queue_complete)
            weather_state=weather(reltime,wlb_gen,wub_gen) #weather(begin_serv,wlb_gen,wub_gen)

            if output==1:
                ee=1


            stepthrough_logger.info(str(AC)+','+str(perm_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(Trav_Time[AC])+','+str(perm_queue_complete)+',')


            #perm_weather_state=weather(reltime,wlb,wub)
            AC_FinishTime, straight_into_service=Gamma_GetServ_Future(k, Time_Sep, reltime,ServTime[AC],Trav_Time[AC],perm_prev_class,perm_class,perm_queue_complete,weather_state, w_rho)

            Opt_Listj[3][index]=(1-gam)*Opt_Listj[3][index]+gam*straight_into_service

            permcost+=getcost(Ac_Infoi.orig_sched_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2) #Ac_Infoi.passenger_weight*(AC_FinishTime-(Ac_Infoi.ps_time+thres))**2
            latest_tm=reltime

            stepthrough_logger.info(str(AC_FinishTime)+','+str(Ac_Infoi.passenger_weight)+','+str(getcost(Ac_Infoi.ps_time,ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi.passenger_weight,thres1,thres2, lam1, lam2))+'\n')

            perm_queue_complete=AC_FinishTime
            perm_prev_class=perm_class

        stepthrough_logger.info('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j][2])+',')

        # if Opt_Listj[2]>0:
        # 	pct_diff=abs(permcost/Opt_Listj[2]-1)
        # 	if pct_diff>iter_max_d:
        # 		iter_max_d=pct_diff
        # else:
        # 	iter_max_d=1

        Opt_Listj[2]=(1-gam)*Opt_Listj[2]+gam*permcost
        Opt_Listj[4]+=permcost**2 #Sum of squares

        stepthrough_logger.info('Total cost: '+','+str(Opt_Listj[2])+','+'Queue probs: '+','+str(Opt_Listj[3])+'\n'+'\n')

    Opt_List.sort(key=lambda x: x[0])
    for j in range(len(Opt_List)):
        step_summ_logger.info(str(Opt_List[j][2])+',')
    step_summ_logger.info('\n')

    # minperm=0
    # mincost=0
    # for j in range(len(GA_Info)):
    # 	if j==0 or GA_Info[j][2]<mincost:
    # 		mincost=GA_Info[j][2]
    # 		minperm=j

    GA_Info.sort(key=lambda x: x[2])
    Opt_List.sort(key=lambda x: x[2])

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

        for j in range(len(GA_Info)):

            GA_Infoj=GA_Info[j]

            if GA_Infoj[2]>0:

                mn1=GA_Infoj[2]
                n1=GA_Infoj[1]
                var1=(GA_Infoj[4]-(n1*mn1**2))/(n1-1)

                for m in range(len(GA_Info)):

                    GA_Infom=GA_Info[m]

                    if GA_Infom[2]>0:

                        mn2=GA_Infom[2]
                        n2=GA_Infom[1]
                        var2=(GA_Infom[4]-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:


                        # 	step_new_logger.info(str(GA_Infoj[0])+','+' loses to '+','+str(GA_Infom[0])+'\n')

                            GA_Infoj[2]=-1
                            break

                        elif mn2>mn1+w_val:


                        # 	step_new_logger.info(str(GA_Infom[0])+','+' loses to '+','+str(GA_Infoj[0])+'\n')

                            GA_Infom[2]=-1

            if GA_Infoj[2]>0:

                for m in range(len(Opt_List)):

                    Opt_Listm=Opt_List[m]

                    if Opt_Listm[2]>0:

                        mn2=Opt_Listm[2]
                        n2=Opt_Listm[1]
                        var2=(Opt_Listm[4]-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:

                        # 	step_new_logger.info(str(GA_Infoj[0])+','+' loses to '+','+str(Opt_Listm[0])+'\n')

                            GA_Infoj[2]=-1
                            break

                        elif mn2>mn1+w_val:


                        # 	step_new_logger.info(str(Opt_Listm[0])+','+' loses to '+','+str(GA_Infoj[0])+'\n')

                            Opt_Listm[2]=-1

        for j in range(len(Opt_List)):

            Opt_Listj=Opt_List[j]

            if Opt_Listj[2]>0:

                mn1=Opt_Listj[2]
                n1=Opt_Listj[1]
                var1=(Opt_Listj[4]-(n1*mn1**2))/(n1-1)

                for m in range(len(GA_Info)):

                    GA_Infom=GA_Info[m]

                    if GA_Infom[2]>0:

                        mn2=GA_Infom[2]
                        n2=GA_Infom[1]
                        var2=(GA_Infom[4]-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:

                        # 	step_new_logger.info(str(Opt_Listj[0])+','+' loses to '+','+str(GA_Infom[0])+'\n')

                            Opt_Listj[2]=-1
                            break

                        elif mn2>mn1+w_val:

                        # 	step_new_logger.info(str(GA_Infom[0])+','+' loses to '+','+str(Opt_Listj[0])+'\n')

                            GA_Infom[2]=-1

            if Opt_Listj[2]>0:

                for m in range(len(Opt_List)):

                    Opt_Listm=Opt_List[m]

                    if Opt_Listm[2]>0:

                        mn2=Opt_Listm[2]
                        n2=Opt_Listm[1]
                        var2=(Opt_Listm[4]-(n2*mn2**2))/(n2-1)

                        w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

                        if mn1>mn2+w_val:

                        # 	step_new_logger.info(str(Opt_Listj[0])+','+' loses to '+','+str(Opt_Listm[0])+'\n')

                            Opt_Listj[2]=-1
                            break

                        elif mn2>mn1+w_val:

                        # 	step_new_logger.info(str(Opt_Listm[0])+','+' loses to '+','+str(Opt_Listj[0])+'\n')
                            Opt_Listm[2]=-1

        j=0
        while j<len(GA_Info):
            GA_Infoj=GA_Info[j]
            if GA_Infoj[2]<0:
                soln_evals_tot+=GA_Infoj[1]
                soln_evals_num+=1
            # 	step_new_logger.info('soln_evals_tot: '+str(soln_evals_tot)+'\n')
            # 	step_new_logger.info('soln_evals_num: '+str(soln_evals_num)+'\n')
            # 	step_new_logger.info('soln_evals_avg: '+str(soln_evals_tot/soln_evals_num)+'\n')
                GA_Info.remove([GA_Infoj[0],GA_Infoj[1],GA_Infoj[2],GA_Infoj[3],GA_Infoj[4]])
            # 	step_new_logger.info('Removed entry '+str(j)+' from GA_PopList'+'\n')
            else:
                step_new_logger.info('Retained sequence '+','+str(GA_Infoj)+','+' in GA_PopList'+'\n')
                j+=1

        j=0
        while j<len(Opt_List):
            Opt_Listj=Opt_List[j]
            if Opt_Listj[2]<0:
                soln_evals_tot+=Opt_Listj[1]
                soln_evals_num+=1
            # 	step_new_logger.info('soln_evals_tot: '+str(soln_evals_tot)+'\n')
            # 	step_new_logger.info('soln_evals_num: '+str(soln_evals_num)+'\n')
            # 	step_new_logger.info('soln_evals_avg: '+str(soln_evals_tot/soln_evals_num)+'\n')
                Opt_List.remove([Opt_Listj[0],Opt_Listj[1],Opt_Listj[2],Opt_Listj[3],Opt_Listj[4]])
            # 	step_new_logger.info('Removed entry '+str(j)+' from Opt_List'+'\n')
            else:
                step_new_logger.info('Retained sequence '+','+str(Opt_Listj)+','+' in Opt_List'+'\n')
                j+=1

        if len(GA_Info)+len(Opt_List)<=Opt_Size:

            solns_left=len(GA_Info)+len(Opt_List)

            soln_evals_tot+=(solns_left*GA_LoopSize)
            soln_evals_num+=solns_left
            pruned=1

        # #Old R&S code is underneath this part

    # 	step_new_logger.info('GA_counter: '+','+str(GA_counter)+'\n'+'\n')
    # 	step_new_logger.info('GA_PopList:'+'\n')

        # #Find the lowest upper CI value
        # min_uppb=0
        # for j in range(len(GA_Info)):
        # 	GA_Infoj=GA_Info[j]
        # 	test_var=(GA_Infoj[4]-(GA_Infoj[1]*(GA_Infoj[2])**2))/(GA_Infoj[1]-1)
        # 	test_lowb=GA_Infoj[2]-1.96*math.sqrt(test_var/GA_Infoj[1])
        # 	test_uppb=GA_Infoj[2]+1.96*math.sqrt(test_var/GA_Infoj[1])
    # 		step_new_logger.info('Seq: '+str(GA_Infoj[0])+','+'Evals: '+str(GA_Infoj[1])+' Mn: '+str(GA_Infoj[2])+' qp[0]: '+str(GA_Infoj[3][0])+' LowB: '+str(test_lowb)+' UppB: '+str(test_uppb)+'\n')
        # 	if j==0 or test_uppb<min_uppb:
        # 		min_uppb=test_uppb

    # 	step_new_logger.info('\n'+'Opt_List:'+'\n')

        # for j in range(len(Opt_List)):
        # 	Opt_Listj=Opt_List[j]
        # 	test_var=(Opt_Listj[4]-(Opt_Listj[1]*(Opt_Listj[2])**2))/(Opt_Listj[1]-1)
        # 	test_lowb=Opt_Listj[2]-1.96*math.sqrt(test_var/Opt_Listj[1])
        # 	test_uppb=Opt_Listj[2]+1.96*math.sqrt(test_var/Opt_Listj[1])
    # 		step_new_logger.info('Seq: '+str(Opt_Listj[0])+','+'Evals: '+str(Opt_Listj[1])+' Mn: '+str(Opt_Listj[2])+' qp[0]: '+str(Opt_Listj[3][0])+' LowB: '+str(test_lowb)+' UppB: '+str(test_uppb)+'\n')
        # 	if test_uppb<min_uppb:
        # 		min_uppb=test_uppb

    # 	step_new_logger.info('\n'+' Min_uppb is '+str(min_uppb)+'\n'+'\n')

        # j=len(GA_Info)-1
        # while j>=0:
        # 	GA_Infoj=GA_Info[j]
        # 	test_var=(GA_Infoj[4]-(GA_Infoj[1]*(GA_Infoj[2])**2))/(GA_Infoj[1]-1)
        # 	test_lowb=GA_Infoj[2]-1.96*math.sqrt(test_var/GA_Infoj[1])
        # 	if test_lowb>min_uppb:
        # 		GA_Info.remove([GA_Infoj[0],GA_Infoj[1],GA_Infoj[2],GA_Infoj[3],GA_Infoj[4]])
        # 		soln_evals_tot+=GA_Infoj[1]
        # 		soln_evals_num+=1
    # 			step_new_logger.info('Removed entry '+str(j)+' from GA_PopList'+'\n')
        # 		if len(GA_Info)+len(Opt_List)==Opt_Size:
        # 			soln_evals_tot+=(Opt_Size*GA_LoopSize)
        # 			soln_evals_num+=Opt_Size
        # 			pruned=1
        # 			break
        # 	j+=-1

        # if len(GA_Info)+len(Opt_List)>Opt_Size:

        # 	j=len(Opt_List)-1
        # 	while j>=0:
        # 		Opt_Listj=Opt_List[j]
        # 		test_var=(Opt_Listj[4]-(Opt_Listj[1]*(Opt_Listj[2])**2))/(Opt_Listj[1]-1)
        # 		test_lowb=Opt_Listj[2]-1.96*math.sqrt(test_var/Opt_Listj[1])
        # 		if test_lowb>min_uppb:
        # 			Opt_List.remove([Opt_Listj[0],Opt_Listj[1],Opt_Listj[2],Opt_Listj[3],Opt_Listj[4]])
        # 			soln_evals_tot+=Opt_Listj[1]
        # 			soln_evals_num+=1
    # 				step_new_logger.info('Removed entry '+str(j)+' from Opt_List'+'\n')
        # 			if len(GA_Info)+len(Opt_List)==Opt_Size:
        # 				soln_evals_tot+=(Opt_Size*GA_LoopSize)
        # 				soln_evals_num+=Opt_Size
        # 				pruned=1
        # 				break
        # 		j+=-1

        if GA_counter>=GA_LoopSize:
            GA_CheckSize=GA_Check_Increment
            solns_left=len(GA_Info)+len(Opt_List)
            soln_evals_tot+=(solns_left*GA_LoopSize)
            soln_evals_num+=solns_left
        else:
            GA_CheckSize+=GA_Check_Increment
        step_new_logger.info('New GA_CheckSize: '+str(GA_CheckSize)+'\n'+'\n')

    Ac_added=[]

    # if len(Opt_List)>0:
    if len(Arr_Pool)>0:

        if len(Opt_List)>0 and len(GA_Info)>0:
            if Opt_List[0][2]<GA_Info[0][2]:
                perm=Opt_List[0]
            else:
                perm=GA_Info[0]
        elif len(Opt_List)>0:
            perm=Opt_List[0]
        elif len(GA_Info)>0:
            perm=GA_Info[0]
        else:
            assert 1==2

        if perm[0][0] in Arr_Pool:

            counter=perm[1]
            qp=perm[3][0]

            #perm=Opt_List[0]
            if counter>=GA_Check_Increment or pruned==1:
                if qp>0: #0.05:
                    j=0
                    while perm[3][j]>0: #0.05:
                        AC=perm[0][j]
                        Ac_added.append(AC)
                        step_new_logger.info('Counter is '+','+str(counter)+', ss_prob is '+','+str(perm[3][j])+', Adding AC '+','+str(AC)+'\n')
                        j+=1
                        if j==len(perm[0]):
                            break

        else:
            counter=0
            qp=0

    else:
        counter=0
        qp=0

    # end_time=time.time()
    # elap=(end_time-start_time)/conv_factor
    # #elap=0.001

    # if iter_max_d<max_d:
    # 	max_d=iter_max_d
    #print('max_d: '+str(max_d))

    #elap=fixed_elap/conv_factor

    # if len(Arr_Pool)+len(Ac_queue)>0:
    # 	print('Out of Genetic')

    return Ac_added,counter,qp,max_d,pruned,GA_CheckSize,GA_counter,soln_evals_tot,soln_evals_num
