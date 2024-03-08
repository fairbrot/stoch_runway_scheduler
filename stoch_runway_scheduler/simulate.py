from typing import List, TextIO
import logging
import math
import random
from .utils import weather, getcost
from .gamma import sample_cond_gamma
from .annealing_cost import Annealing_Cost

def generate_trajectory(Dep_time: float, Ps_time: float, tau: int, wiener_sig: float, freq: int = 100):
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
        if j > Dep_time: # only update ETA if we've gone beyond the AC's departure time - JF: I think logic is wrong here - j needs scaling (like below)
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

def Calculate_FCFS(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,pool_max,list_min,wlb_tm,wub_tm, NoA: int, NormalApprox, w_rho: float, k: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float):

    #start_time=time.time()

    tm=0

    totserv=0
    totcost=0
    prev_class=4
    queue_complete=0
    weather_state=0

    FCFS_Seq=[0]*NoA
    for i in range(NoA):
        FCFS_Seq[i]=ArrTime_Sorted[i][1]

    FCFS_cost=Annealing_Cost(FCFS_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0, NoA, NormalApprox, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)

    return FCFS_cost


def Posthoc_Check(seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,output, NoA: int, NormalApprox, w_rho: float, k: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float):

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
        release_time=Ac_Infoi[4] #max(latest_tm,ArrTime[AC][0])
        #print('j: '+str(j)+' AC: '+str(AC)+' release_time: '+str(release_time))
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

def Get_Actual_Serv(AC,prev_class,cur_class,weather_state,k, Time_Sep: List[List[int]], norm_approx_min, Ac_Info, w_rho: float):

    if k<norm_approx_min:
        #servtime=0
        serv_percs=Ac_Info[AC][7]
        if weather_state==1:
            ws=1/w_rho
        else:
            ws=1

        rate=ws*k/(Time_Sep[prev_class][cur_class]/60)

        servtime=serv_percs/rate #Transformation causes serv_percs to go from [mean k, var k] to [mean e_{ij}, var e_{ij}^2/k]

        # for j in range(k):
        # 	servtime+=(-1/rate)*math.log(serv_percs[j])

    else:
        if weather_state==1:
            Mn=w_rho*Time_Sep[prev_class][cur_class]/60
        else:
            Mn=Time_Sep[prev_class][cur_class]/60
        SD=math.sqrt(Mn**2/k)
        servtime=Ac_Info[AC][7]*SD+Mn
        # u=int(z*10000)
        # servtime=normcdf[u]*SD+Mn

    #print('For AC '+str(AC)+', prev class '+str(prev_class)+', current class '+str(cur_class)+', weather state '+str(weather_state)+', we calculated the actual service time as '+str(servtime))

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

# def GetServTime(trav_time,serv_time,rel_time,prev_class,cur_class,tm,ph_B,ee,weather_state): #old, disused version

# 	#print('tm: '+str(tm))

# 	t1=rel_time+trav_time

# 	#Next, get s1+Z2
# 	# if type(prev_class)==list:
# 	# 	print('prev_class: '+str(prev_class))
# 	# if type(cur_class)==list:
# 	# 	print('cur_class: '+str(cur_class))

# 	rate=k/(Time_Sep[prev_class][cur_class]/60)
# 	if weather_state==1:
# 		rate*=1/w_rho #=0.5

# 	t2=tm
# 	for m in range(ph_B):
# 		# if m>len(serv_time)-1:
# 		# 	print('ph_B: '+str(ph_B))
# 		t2+=(-1/rate)*math.log(serv_time[m])

# 	if ee==1 and stepthrough==1:
# 		stepthrough_logger.info(str(t2-tm)+',')

# 	# if len(Ac_queue)==2 and ee==1:
# 	# 	print('Trav time completed at '+str(t1)+', serv time completed at '+str(t2))

# 	if t1<t2:
# 		straight_into_service=0
# 		t_out=t2
# 	else:
# 		straight_into_service=1
# 		t_out=t1

# 	return t_out,straight_into_service



def Update_ETAs(Ac_Info,Arr_NotReady,Dep_NotReady,Ac_queue,tm,Brown_Motion, Arr_Pool, tau):
    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    step_new_logger = logging.getLogger('step_new')
    #print('Entered Update_ETAs')

    i=0
    while i<=len(Arr_NotReady)-1:
        AC=Arr_NotReady[i]
        Ac_Infoi=Ac_Info[AC]
        #print('AC: '+str(Arr_NotReadyi)+' tm: '+str(tm)+' Ac_Infoi[9]: '+str(Ac_Infoi[9]))
        if tm>=Ac_Infoi[9]:
            Arr_Pool.append(AC)
            #print('* Added aircraft '+str(AC)+' to the arrival pool at time '+str(tm)+' (new readiness time is '+str(Ac_Infoi[9]+tau)+')')
            stepthrough_logger.info('* Added aircraft '+str(AC)+' to the arrival pool at time '+str(tm)+' (new readiness time is '+str(Ac_Infoi[9]+tau)+')'+'\n'+'\n')
            step_summ_logger.info('* Added aircraft '+str(AC)+' to the arrival pool at time '+str(tm)+' (new readiness time is '+str(Ac_Infoi[9]+tau)+')'+'\n'+'\n')
            step_new_logger.info('* Added aircraft '+str(AC)+' to the arrival pool at time '+str(tm)+' (new readiness time is '+str(Ac_Infoi[9]+tau)+')'+'\n'+'\n')
            Arr_NotReady.remove(AC)
            Ac_Infoi[0]+=1
            Ac_Infoi[3]=Ac_Infoi[9]+tau
            i+=-1 #because we want to negate the "i+=1" below
        else:
            Ac_Infoi[3]=Brown_Motion[AC][int(tm*100)]
        i+=1

    i=0
    while i<=len(Ac_queue)-1:
        AC=Ac_queue[i]
        Ac_Infoi=Ac_Info[AC]
        if Ac_Infoi[11]==0:
            rel_time=Ac_Infoi[4]
            trav_so_far=tm-rel_time #amount of time spent travelling to the runway so far
            rounded_trav_so_far=round(trav_so_far,2)
            if rounded_trav_so_far>trav_so_far:
                rounded_trav_so_far+=-0.01
            trav_time=Ac_Infoi[6]
            if rounded_trav_so_far>=trav_time:
                #print('* Aircraft '+str(AC)+' has finished its travel time at '+str(tm)+'; still waiting for service time.')
                stepthrough_logger.info('* Aircraft '+str(AC)+' has finished its travel time at '+str(tm)+'; still waiting for service time.'+'\n'+'\n')
                step_summ_logger.info('* Aircraft '+str(AC)+' has finished its travel time at '+str(tm)+'; still waiting for service time.'+'\n'+'\n')
                Ac_Infoi[11]=1
            else:
                Ac_Infoi[3]=Brown_Motion[AC][int((Ac_Infoi[9]+rounded_trav_so_far)*100)]
        #print('tm: '+str(tm)+' AC: '+str(AC)+' Ac_Infoi[9]: '+str(Ac_Infoi[9])+' rel_time: '+str(rel_time)+' rounded_trav_so_far: '+str(rounded_trav_so_far)+' Ac_Infoi[3]: '+str(Ac_Infoi[3]))
        i+=1

    # i=0
    # while i<=len(Dep_NotReady)-1:
    # 	Dep_NotReadyi=Dep_NotReady[i]
    # 	random.seed(int((Dep_NotReadyi+1)*(repn+1)*(tm+delta)*100000)) #note this is included in order to enable consistency in comparisons between FCFS, Perm and PI
    # 	Ac_Infoi=Ac_Info[Dep_NotReadyi]
    # 	new_readiness_time=random.gauss(Ac_Infoi[3],delta*wiener_sig)
    # 	if new_readiness_time-tau<=tm+delta:
    # 		Dep_Pool.append(Dep_NotReadyi)
    # 		print('* Added aircraft '+str(Dep_NotReadyi)+' to the departure pool at time '+str(tm+delta)+' (new readiness time is '+str(new_readiness_time)+')')
    # 		Dep_NotReady.remove(Dep_NotReadyi)
    # 		Ac_Infoi[0]+=1
    # 		i+=-1 #because we want to negate the "i+=1" below
    # 	Ac_Infoi[3]=new_readiness_time
    # 	i+=1

    #print('Exited Update_ETAs')

def Update_Stats(tm,AC,Ac_Info,Ac_queue,real_queue_complete,wlb_tm,wub_tm,latest_class,Ov_GA_counter,next_completion_time, k: int, Time_Sep: List[List[int]], norm_approx_min, w_rho: float, SubPolicy, counter, qp):

    Ac_Infoi=Ac_Info[AC]

    Ac_Infoi[0]+=1 #update status
    release_time=tm #release time
    begin_serv=max(release_time,real_queue_complete) #time that service begins
    cur_class=Ac_Infoi[1]

    get_weather_state=weather(release_time,wlb_tm,wub_tm) #weather(begin_serv,wlb_tm,wub_tm)

    actual_serv=Get_Actual_Serv(AC,latest_class,cur_class,get_weather_state,k,Time_Sep, norm_approx_min, Ac_Info, w_rho)
    trav_time=Ac_Infoi[6]

    t1=release_time+trav_time
    t2=begin_serv+actual_serv

    finish_time=max(t1,t2)
    real_queue_complete=finish_time

    Ac_Infoi[4]=release_time
    Ac_Infoi[5]=begin_serv
    Ac_Infoi[12]=get_weather_state
    Ac_Infoi[8]=actual_serv
    Ac_Infoi[16]=finish_time

    if SubPolicy in ('GA','GAD','VNS','VNSD'):
        Ac_Infoi[13]=Ov_GA_counter
        Ov_GA_counter=0
    else:
        Ac_Infoi[13]=counter
    Ac_Infoi[14]=qp

    latest_class=cur_class

    if len(Ac_queue)==1:
        next_completion_time=finish_time

    return real_queue_complete,next_completion_time,latest_class,Ov_GA_counter

def Serv_Completions(Ac_Info,Ac_queue,prev_class,totserv,Ac_finished,tm,next_completion_time, thres1: int, thres2: int, lam1: float, lam2: float, f: TextIO, SubPolicy, rep, Time_Sep: List[List[int]], Left_queue):

    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    # print('Entered Serv_Completions')
    # print('ESC tm: '+str(tm)+' next_completion_time: '+str(next_completion_time))

    arr_cost=0
    dep_cost=0

    j=0

    while len(Ac_queue)>0:

        AC=Ac_queue[0]
        Ac_Infoi=Ac_Info[AC]

        finish_time=Ac_Infoi[16]
        current_class=Ac_Infoi[1]

        #print('finish_time: '+str(finish_time))

        if tm>=finish_time: #release_time+trav_time and phase==k:
            #print('* Service phase '+str(phase)+' completed for aircraft '+str(Ac_queue[0])+' at time '+str(tm+delta))
            Ac_finished[AC]=finish_time
            #print('* Service completion finished for aircraft '+str(AC))
            stepthrough_logger.info('* Service completion finished for aircraft '+str(AC)+'\n'+'\n')
            step_summ_logger.info('* Service completion finished for aircraft '+str(AC)+'\n'+'\n')
            if Ac_Infoi[0]==2:
                Ac_Infoi[0]=3
                arr_cost+=getcost(Ac_Infoi[18],Ac_Infoi[9],Ac_Infoi[6],finish_time,Ac_Infoi[10],thres1,thres2, lam1, lam2)
                #print('* Cost incurred is '+str(arr_cost))
                totserv+=1
            else:
                Ac_Infoi[0]=6
                arr_cost+=getcost(Ac_Infoi[18],Ac_Infoi[9],Ac_Infoi[6],finish_time,Ac_Infoi[10],thres1,thres2, lam1, lam2)

            f.write(str(SubPolicy)+','+str(rep)+','+str(AC)+','+str(Ac_Infoi[19])+','+str(prev_class)+','+str(current_class)+','+str(Time_Sep[prev_class][current_class]/60)+','+str(Ac_Infoi[18])+','+str(Ac_Infoi[2])+','+str(Ac_Infoi[9])+','+str(Ac_Infoi[4])+','+str(Ac_Infoi[6])+','+str(Ac_Infoi[12])+','+str(Ac_Infoi[5])+','+str(Ac_Infoi[8])+','+str(Ac_Infoi[16])+','+str(max(0,finish_time-(Ac_Infoi[2]+thres1)))+','+str(finish_time-(Ac_Infoi[9]+Ac_Infoi[6]))+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],Ac_Infoi[9],Ac_Infoi[6],finish_time,Ac_Infoi[10],thres1,thres2, lam1, lam2))+',')
            f.write(str(Ac_Infoi[13])+','+str(Ac_Infoi[14])+',')

            if SubPolicy in ('GA','GAD','VNS','VNSD'):
                f.write(str(Ac_Infoi[15])+',')
            # if SubPolicy in ('VNS'):
            # 	f.write(str(Loop_Evals/Loop_Nums)+',')
            # elif SubPolicy in ('VNSD'):
            # 	f.write(',')
            # if SubPolicy in ('GA','GAD','VNS','VNSD'):
            # 	f.write(str(elap_tot/elap_num)+','+str(Repop_elap_tot/Repop_elap_num)+','+str(Pop_elap_tot/Pop_elap_num)+',')
                #f.write(str(elap)+','+str(Repop_elap)+','+str(Pop_elap)+',')
            f.write('\n')

            prev_class=current_class

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
