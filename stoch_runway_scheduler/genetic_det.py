from typing import List
import time
from .utils import weather, FlightInfo, Cost
from .sequence import SequenceInfo
from .gamma import gamma_cond_exp

# JF: This is the main deterministic heuristic
# Name may not be best choice
def Genetic_determ(Ac_Info: List[FlightInfo], Arr_Pool: List[int], Arr_NotReady: List[int], 
                    Ac_queue: List[int], Left_queue: List[int], tm: float, NoA: int, k:int, 
                    prev_class: int, GA_Info: List[SequenceInfo], wlb, wub, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], cost_fn: Cost, 
                    tot_arr_cost: float, tot_dep_cost: float, w_rho: float, stepthrough:int, step_summ:int, step_new: int):

    output = 0 # output == 1 means we're printing results as we go along; output == 2 means we're outputting results to "Detailed" csv file
    ee = 0
    if stepthrough == 1:
        ee = 1
    start_time=time.time()

    if stepthrough==1:
        st.write('Now entering Genetic procedure'+'\n')

    ArrTime=[0]*NoA
    ServTime=[0]*NoA
    Trav_Time=[0]*NoA

    # Generate arrival and service time percentiles for AC not yet in queue

    ArrTime_Sorted=[]

    for AC in Arr_Pool:
        ArrTime[AC] = Ac_Info[AC].pool_time
        ArrTime_Sorted.append([ArrTime[AC],AC])

    for AC in Arr_NotReady:
        ArrTime[AC]=max(0, Ac_Info[AC].eta - tau)
        ArrTime_Sorted.append([ArrTime[AC],AC])

    ArrTime_Sorted.sort(key=lambda x: x[0])

    basecost=tot_arr_cost+tot_dep_cost
    if stepthrough==1:
        st.write('basecost is '+','+str(basecost)+'\n'+'\n')
        st.write('Generated results for ACs already in the queue are as follows:'+'\n')
        st.write('AC'+','+'Class'+','+'Time Sep'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

    if len(Ac_queue)>0:

        #Need to generate service times for AC already in the queue; first consider the customer in position 0



        AC=Ac_queue[0]
        Ac_Infoi=Ac_Info[AC]
        rel_time=Ac_Infoi.release_time
        sv_time=Ac_Infoi.enters_service
        cur_class=Ac_Infoi.ac_class

        t1=Ac_Infoi.eta

        #Get the conditional expectation of service time based on service time elapsed so far

        if Ac_Infoi.weather_state==1:
            beta=k/(w_rho*Time_Sep[prev_class][cur_class]/60)
        else:
            beta=k/(Time_Sep[prev_class][cur_class]/60)
        alpha=k
        cond_tm=gamma_cond_exp(tm-sv_time,alpha,beta)

        t2=sv_time+cond_tm #tm+(cond_tm-sv_time)
        #print('Time remaining based on mean value: '+str(sv_time+cond_tm-tm))

        queue_complete=max(t1,t2)

        basecost += cost_fn(Ac_Infoi.orig_sched_time,Ac_Infoi.pool_time,tau,queue_complete,Ac_Infoi.passenger_weight)
        if stepthrough==1:
            st.write(str(queue_complete)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time, Ac_Infoi.pool_time, tau, queue_complete, Ac_Infoi.passenger_weight))+'\n')

        perm_prev_class=cur_class\


        #Now consider the rest of the customers in the queue
        for j in range(1,len(Ac_queue)):

            AC=Ac_queue[j]
            Ac_Infoi=Ac_Info[AC]
            rel_time=Ac_Infoi.release_time
            cur_class=Ac_Infoi.ac_class

            if stepthrough==1:
                st.write(str(AC)+','+str(Ac_Infoi.ac_class)+','+str(Time_Sep[perm_prev_class][cur_class]/60)+','+str(Ac_Infoi.release_time)+','+str(Ac_Infoi.travel_time)+','+str(Ac_Infoi.enters_service)+',')

            t1=Ac_Infoi.eta
            weather_state=weather(rel_time,wlb,wub) #weather(queue_complete,wlb,wub)
            if weather_state==1: #Ac_Infoi.weather_state==1:
                t2=queue_complete+(w_rho*Time_Sep[perm_prev_class][cur_class]/60)
                if stepthrough==1:
                    st.write(str(w_rho*Time_Sep[perm_prev_class][cur_class]/60)+',')
            else:
                t2=queue_complete+(Time_Sep[perm_prev_class][cur_class]/60)
                if stepthrough==1:
                    st.write(str(Time_Sep[perm_prev_class][cur_class]/60)+',')
            queue_complete=max(t1,t2)

            perm_prev_class=cur_class
            basecost+=cost_fn(Ac_Infoi.orig_sched_time,Ac_Infoi.pool_time,tau,queue_complete,Ac_Infoi.passenger_weight)

            if stepthrough==1:
                st.write(str(queue_complete)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time,Ac_Infoi.pool_time,tau,queue_complete,Ac_Infoi.passenger_weight))+'\n')

    else:

        queue_complete=tm
        perm_prev_class=prev_class

    stored_prev_class=perm_prev_class
    stored_queue_complete=queue_complete

    #Try all the sequences in the population

    for j in range(len(GA_Info)):

        if stepthrough==1:
            st.write('\n'+'Now trying sequence '+','+str(GA_Info[j].sequence)+'\n')
            st.write('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

        permcost=basecost
        latest_tm=tm
        perm_prev_class=stored_prev_class
        perm_queue_complete=queue_complete
        #perm_weather_state=weather_state

        perm=GA_Info[j].sequence
        GA_Info[j].n_traj+=1
        #gam=0.01
        gam=1/GA_Info[j].n_traj

        #print('GA_Info: '+str(GA_Info))

        no_ACs=min(Max_LookAhead,len(perm)) #min(AC_List_Length,len(Arr_Pool)+len(Arr_NotReady))
        for index in range(no_ACs):

            #index=perm[i]
            AC=perm[index]
            Ac_Infoi=Ac_Info[AC]
            perm_class=Ac_Infoi.ac_class
            reltime=max(latest_tm,ArrTime[AC])

            if stepthrough==1:
                st.write(str(AC)+','+str(Ac_Infoi.ac_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(tau)+','+str(perm_queue_complete)+',')

            t1=reltime+tau
            if reltime>=wlb and reltime<=wub:
                exp_serv=w_rho*Time_Sep[perm_prev_class][perm_class]/60
            else:
                exp_serv=Time_Sep[perm_prev_class][perm_class]/60
            t2=perm_queue_complete+exp_serv
            if stepthrough==1:
                st.write(str(exp_serv)+',')

            AC_FinishTime=max(t1,t2)
            if t1>=t2:
                straight_into_service=1
            else:
                straight_into_service=0

            GA_Info[j].queue_probs[index]=(1-gam)*GA_Info[j].queue_probs[index]+gam*straight_into_service

            permcost += cost_fn(Ac_Infoi.orig_sched_time, ArrTime[AC], tau, AC_FinishTime, Ac_Infoi.passenger_weight)
            latest_tm = reltime

            if stepthrough==1:
                st.write(str(AC_FinishTime)+','+str(Ac_Infoi.passenger_weight)+','+str(cost_fn(Ac_Infoi.ps_time, ArrTime[AC], tau, AC_FinishTime, Ac_Infoi.passenger_weight))+'\n')

            perm_queue_complete = AC_FinishTime
            perm_prev_class = perm_class

        if stepthrough==1:
            st.write('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j].v)+',')

        GA_Info[j].v=(1-gam)*GA_Info[j].v+gam*permcost
        if stepthrough==1:
            st.write('Total cost: '+','+str(GA_Info[j].v)+','+'Queue probs: '+','+str(GA_Info[j].queue_probs)+'\n'+'\n')

    if step_summ==1:
        GA_Info.sort(key=lambda x: x.sequence)
        for j in range(len(GA_Info)):
            st2.write(str(GA_Info[j].v)+',')
        st2.write('\n')

    GA_Info.sort(key=lambda x: x.v)

    Ac_added=[]
    counter=0
    qp=0

    if len(Arr_Pool)>0:
        assert len(GA_Info) > 0
        perm = GA_Info[0]

        if perm[0][0] in Arr_Pool:

            if 1==1: #perm[1]>=100 and 1==1: #perm[3][0]>0: #0.05:
                j=0
                counter=perm[1]
                qp=perm[3][0]
                while 1==1: #perm[3][j]>0: #0.05:
                    AC=perm[0][j]
                    Ac_added.append(AC)
                    j+=1
                    if j==len(perm[0]):
                        break

    # end_time=time.time()
    # elap=(end_time-start_time)/conv_factor
    # #elap=0.01

    #elap=fixed_elap_vnsd/conv_factor

    return Ac_added,counter,qp,stored_queue_complete
