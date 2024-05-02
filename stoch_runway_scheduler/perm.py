from typing import List, TextIO
import random
import math
import itertools
from .utils import Cost
from .weather import WeatherProcess
from .annealing_cost import Annealing_Cost

def Perm_Heur_New(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process: WeatherProcess, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], cost_fn: Cost):

    #start_time=time.time()

    tm = 0

    totserv = 0
    totcost = 0
    prev_class = 4
    queue_complete = 0
    weather_state = 0

    Anneal_Seq=[0]*NoA
    for i in range(NoA):
        Anneal_Seq[i] = ArrTime_Sorted[i][1]

    NewCost = Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted, weather_process,0, NoA, w_rho, k, Time_Sep, cost_fn)
    OptCost = NewCost
    Opt_Seq = Anneal_Seq[:]

    #T=1000
    iter_no=0
    runtime_chk=0
    AC_remaining=NoA
    c=0
    perm_size=6 #no. of ACs to shuffle on one iteration

    while runtime_chk==0:

        #Select a starting point in the sequence 
        start_pos=int(random.random()*(AC_remaining-perm_size+1))

        perm=[]
        for j in range(start_pos,start_pos+perm_size):
            perm.append(Anneal_Seq[j])

        random.shuffle(perm)

        for j in range(perm_size):
            Anneal_Seq[start_pos+j]=perm[j]

        NewCost = Annealing_Cost(Anneal_Seq, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process, 0, NoA, w_rho, k, Time_Sep, cost_fn)

        #print('Opt_Seq: '+str(Opt_Seq)+' Cost: '+str(OptCost))
        #print('New Seq: '+str(Anneal_Seq)+' Cost: '+str(NewCost))

        if NewCost<OptCost:
            Opt_Seq=Anneal_Seq[:]
            OptCost=NewCost
            c=0
            #print('Accept!')
        else:
            Anneal_Seq=Opt_Seq[:]
            c+=1
            #print('Reject!')

        if c>=10000:
            runtime_chk=1

    return OptCost,c

def Perm_Heur(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process: WeatherProcess, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], cost_fn: Cost, f1: TextIO):

    tm=0

    AC_Used=[]

    totserv=0
    totcost=0
    prev_class=4
    queue_complete=0
    weather_state=0

    #n_max=6

    while totserv<NoA:

        #print('tm: '+str(tm))

        AC_List=[]

        i=0
        while len(AC_List)<pool_max:
            AC=i
            if ArrTime[AC][0]<=tm and AC not in AC_Used and AC not in AC_List:
                AC_List.append(AC)
            i+=1
            if i==NoA:
                break

        if len(AC_List)<list_min:
            i=0
            while len(AC_List)<list_min:
                AC=ArrTime_Sorted[i][1]
                if AC not in AC_Used and AC not in AC_List:
                    AC_List.append(AC)
                i+=1
                if i==NoA:
                    break

        #print('Perm_Heur AC_List: '+str(AC_List))

        assert len(AC_List)<9

        Perm_List=list(itertools.permutations(AC_List))

        mincost=0
        for ii in range(len(Perm_List)):

            perm=Perm_List[ii]

            perm_cost=0
            latest_tm=tm
            perm_prev_class=prev_class
            perm_queue_complete=queue_complete
            perm_weather_state=weather_state

            for j in range(len(AC_List)):

                AC=perm[j]
                release_time=max(latest_tm,ArrTime[AC][0])
                trav_time=Ac_Info[AC].travel_time
                perm_class=Ac_Info[AC].ac_class
                begin_serv=max(release_time,perm_queue_complete)
                perm_weather_state=weather_process(release_time)  #JF Question: why not begin_serv?


                if perm_weather_state==1:
                    ws=1/w_rho
                else:
                    ws=1
                rate=ws*k/(Time_Sep[perm_prev_class][perm_class]/60)
                serv=ServTime[AC]/rate

                t1=release_time+trav_time
                t2=perm_queue_complete+serv
                perm_queue_complete=max(t1,t2)

                perm_cost += cost_fn(Ac_Info[AC].orig_sched_time,ArrTime[AC][0], Ac_Info[AC].travel_time, perm_queue_complete, Ac_Info[AC].passenger_weight)

                latest_tm=release_time
                perm_prev_class=perm_class

            if perm==Perm_List[0] or perm_cost<mincost:
                minperm=perm
                mincost=perm_cost

            #print('Average cost for permutation '+str(perm)+' is '+str(perm_cost))

        #print('Optimal perm is '+str(minperm)+' with cost of '+str(mincost))

        AC=minperm[0]
        release_time=max(tm,ArrTime[AC][0])
        trav_time=Ac_Info[AC].travel_time
        cur_class=Ac_Info[AC].ac_class
        begin_serv=max(release_time,queue_complete)
        weather_state = weather_process(release_time) # JF Question: why not begin_serv?

        #print('AC: '+str(AC)+' pool arrival: '+str(ArrTime[AC][0])+' release_time: '+str(release_time))

        if weather_state==1:
            ws=1/w_rho
        else:
            ws=1
        rate=ws*k/(Time_Sep[prev_class][cur_class]/60)
        serv=ServTime[AC]/rate
        # serv=0
        # for m in range(k):
        # 	serv+=(-1/rate)*math.log(ServTime[AC][m])

        actual_serv=serv #for outputting
        begin_serv=queue_complete #for outputting

        t1=release_time+trav_time
        t2=queue_complete+serv
        queue_complete=max(t1,t2)

        f1.write(str(AC)+','+str(prev_class)+','+str(cur_class)+','+str(Time_Sep[prev_class][cur_class]/60)+','+str(Ac_Info[AC].ps_time)+','+str(Ac_Info[AC].eta)+','+str(Ac_Info[AC].pool_time)+','+str(release_time)+','+str(trav_time)+',')
        # for j in range(k):
        # 	f1.write(str(Ac_Info[AC].service_rns[j])+',')
        f1.write(str(actual_serv)+','+str(begin_serv)+','+str(queue_complete)+','+str(queue_complete-begin_serv)+','+str(Ac_Info[AC].passenger_weight)+','+str(cost_fn(Ac_Info[AC].ps_time,ArrTime[AC][0],Ac_Info[AC].travel_time,queue_complete,Ac_Info[AC].passenger_weight))+'\n')

        if queue_complete > Ac_Info[AC].ps_time + cost_fn.thres1:
            totcost += cost_fn(Ac_Info[AC].orig_sched_time, ArrTime[AC][0], Ac_Info[AC].travel_time, queue_complete, Ac_Info[AC].passenger_weight)
        tm=release_time
        prev_class=cur_class

        AC_Used.append(minperm[0])
        totserv+=1

    return totcost,AC_Used