from typing import List
import math
from .utils import weather, getcost

def Annealing_Cost(seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,output, NoA: int, NormalApprox, w_rho: float, k: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float):

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
        release_time=max(latest_tm,ArrTime[AC][0])
        trav_time=Ac_Infoi[6]
        perm_class=Ac_Infoi[1]
        begin_serv=max(release_time,perm_queue_complete)
        perm_weather_state=weather(release_time,wlb_tm,wub_tm) #weather(begin_serv,wlb_tm,wub_tm)

        if NormalApprox==0:
            if perm_weather_state==1:
                ws=1/w_rho
            else:
                ws=1
            rate=ws*k/(Time_Sep[perm_prev_class][perm_class]/60)
            serv=ServTime[AC]/rate
            # serv=0
            # for m in range(k):
            # 	serv+=(-1/rate)*math.log(ServTime[AC][m])
        else:
            Mn=Time_Sep[perm_prev_class][perm_class]/60
            if perm_weather_state==1:
                Mn*=w_rho
            SD=math.sqrt(Mn**2/k)
            serv=ServTime[AC]*SD+Mn
            # u=int(z*10000)
            # serv=normcdf[u]*SD+Mn

        t1=release_time+trav_time
        t2=perm_queue_complete+serv
        perm_queue_complete=max(t1,t2)

        perm_cost+=getcost(Ac_Infoi[18],ArrTime[AC][0],trav_time,perm_queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)

        if output==1:
            print('AC: '+str(AC)+' class: '+str(perm_class)+' release_time: '+str(release_time)+' trav_time: '+str(trav_time)+' begin_serv: '+str(begin_serv)+' t1: '+str(t1)+' t2: '+str(t2)+' finish time: '+str(perm_queue_complete)+' weather state: '+str(perm_weather_state)+' pax weight: '+str(Ac_Infoi[10])+' cost: '+str(getcost(Ac_Infoi[2],ArrTime[AC][0],trav_time,perm_queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)))

        latest_tm=release_time
        perm_prev_class=perm_class

        j+=1

    return perm_cost