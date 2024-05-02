from typing import List
import math
from .utils import Cost
from .weather import WeatherProcess

def Annealing_Cost(seq, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process: WeatherProcess, output, NoA: int, w_rho: float, k: int, Time_Sep: List[List[int]], cost_fn: Cost):

    perm=seq
    perm_cost=0
    latest_tm=0
    perm_prev_class=4
    perm_queue_complete=0
    perm_weather_state=0
    j=0

    while j < NoA:

        AC = perm[j]
        Ac_Infoi = Ac_Info[AC]
        release_time = max(latest_tm, ArrTime[AC][0])
        trav_time=Ac_Infoi.travel_time
        perm_class=Ac_Infoi.ac_class
        begin_serv=max(release_time, perm_queue_complete)
        perm_weather_state = weather_process(release_time) # weather(begin_serv) # JF Question - why not the latter?

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