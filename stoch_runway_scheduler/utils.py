from typing import List
from enum import Enum
from dataclasses import dataclass
import csv
import math
import random

#0: not ready yet (arrival), 1: in arrival pool, 2: added to arrival queue, 3: not ready yet (departure), 4: in departure pool, 5: added to departure queue, 6: finished.
class FlightStatus(Enum):
    NOT_READY = 0
    IN_POOL = 1
    IN_QUEUE = 2
    DEP_NOT_READY = 3 # JF Question: should these final two be combined?
    FINISHED = 6

@dataclass
class FlightInfo:
    status: FlightStatus
    ac_class: int 
    ps_time: float # pre-scheduled time (plus pre-tactical delay)
    eta: float # latest ETA
    release_time: float # time at which aircraft is released from pool
    enters_service: float # the time at which aircraft enters service
    travel_time: float # travel time (generated in advance), # JF Question note from entering pool to runway?
    service_rns: List[float] # list of random numbers used to calculate service times
    service_time: float # actual service time s1+Z2 (worked out after class information is known)
    pool_time: float # actual time that they join the pool (generated in advance)
    passenger_weight: float # g_i
    travel_time_indicator: bool #  indicator to show whether or not the AC's travel time has already been completed
    weather_state: int # weather state at the time of release
    counter: int # counter (for Perm only)
    qp: float # Perm only
    pred_cost: float # predicted total cost at time of release
    service_completion_time: float # actual service completion time
    sched_dep_time: float # scheduled departure time
    orig_sched_time: float # original pre-scheduled arrival time before adding pre-tactical delay
    flight_id: int

def read_flight_data(data_fn: str, min_time: int, max_time: int, wiener_sig: float):
    Ac_class = [] # this will store the weight class for each aircraft
    Orig_Ps = [] # original pre-scheduled times of aircraft, before applying the pre-tactical delay
    Dep_Ps = []  # h (i.e. time at which tactical uncertainty begins)
    flight_id = []
    Alpha_Ps = [] # Alpha parameters for Gamma dist for sampling pretactical delays
    Beta_Ps = []  # Beta parameters
    late_means = []
    with open(data_fn, 'r') as csvfile: # data file includes the on-time performance data such as means, variance of lateness based on 1 year of historical info from FlightRadar [DON'T CHANGE THIS FILE]
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata = list(datareader)
        n_lines = len(inputdata)
        for i in range(1, n_lines): # start from 1 because there's a title row
            ps_time = float(inputdata[i][1]) # pre-scheduled time
            # When reading in the flight data from the data file we only want to include flights with a pre-scheduled time between 6AM (360 mins) and 2PM (840 mins), including 6AM but not including 2PM
            if ps_time >= min_time and ps_time < max_time:
                ## Read Data ##
                ac_class = int(inputdata[i][2]) # weight class
                flight_name = str(inputdata[i][3]) # flight number
                dep_time = int(inputdata[i][4]) # departure time from the origin airport
                # scheduled flight dur, i.e. difference between scheduled departure and arrival time
                # JF Question: This doesn't seem to be used at the moment?
                sched_dur = int(inputdata[i][5]) 
                lateness_mn = float(inputdata[i][6]) # mean lateness based on historical data
                lateness_var = float(inputdata[i][7]) # variance of lateness based on historical data

                h = dep_time - 15 # JF Question: this 15 probably shouldn't be hard-coded - I think this is q_i in the paper? Probably yes

                ## Transform times ##
                # HERE WE RE-SCALE TIME SO THAT TIME '6AM' IS COUNTED AS TIME (ZERO+60).
                # WE START SIMULATING FROM TIME ZERO, I.E. AN HOUR BEFORE 6AM.
                h += -min_time+60
                ps_time += -min_time+60
                
                # The equations for xi_bar, si2, h_i, alpha and beta below are for calculating the parameters of the gamma distribution 
                # used for the pre-tactical delay. Details of this method are in Section 4 of the paper.
                xibar = ps_time + lateness_mn
                si2 = lateness_var 

                alpha = ((xibar-h)**2)/(si2-(wiener_sig**2)*(xibar-h)) # Can be negative via denominator
                beta = (xibar-h)/(si2-(wiener_sig**2)*(xibar-h))

                Ac_class.append(ac_class)
                Orig_Ps.append(ps_time) # original pre-scheduled time (as opposed to scheduled time following pre-tactical delay)
                flight_id.append(flight_name)
                Dep_Ps.append(h)
                Alpha_Ps.append(alpha)
                Beta_Ps.append(beta)
                late_means.append(lateness_mn)
            elif ps_time >= max_time: # Indicates that we have got to the end of the set of flights scheduled to arrive by 2PM
                break
    return flight_id, Ac_class, Orig_Ps, Dep_Ps, Alpha_Ps, Beta_Ps, late_means


def weather(tm, wlb, wub): # wlb is starting time for bad weather period, wub is ending time

    if wlb<wub:

        if tm<wlb:
            get_weather_state=0
        elif tm<wub:
            get_weather_state=1
        else:
            get_weather_state=2

    else:

        get_weather_state=0

    return get_weather_state

def getcost(ps_time,pool_time,trav_time,landing_time,pax_weight,thres1,thres2, lam1: float, lam2: float):

    cost=0
    #lam1=0.5 #weight for punctuality
    #lam2=0.5 #weight for queueing HMMM

    if landing_time>ps_time+thres1:
        cost+=lam1*pax_weight*(landing_time-(ps_time+thres1))**2

    if landing_time>pool_time+trav_time+thres2:
        cost+=lam2*pax_weight*(landing_time-(pool_time+trav_time+thres2))**2

    return cost



def Normal_GetServ(rel_time,trav_time,prev_class,cur_class,tm,weather_state, Time_Sep: List[List[int]], w_rho: float, k: int):

    #This is for ACs that are already in the queue but not yet in service

    t1=rel_time+trav_time

    sep=Time_Sep[prev_class][cur_class]/60

    if weather_state==1:
        Mn=w_rho*sep
    else:
        Mn=sep

    SD=math.sqrt(Mn**2/k)

    t2=tm+random.gauss(0,1)*SD+Mn

    if t1<t2:
        straight_into_service=0
        t_out=t2
    else:
        straight_into_service=1
        t_out=t1

    return t_out,straight_into_service


def Normal_GetServ_Future(rel_time,serv_time,trav_time,prev_class,cur_class,tm,weather_state, Time_Sep: List[List[int]], w_rho: float, k: int):

    #This is for ACs that have not yet been added to the queue

    t1=rel_time+trav_time

    sep=Time_Sep[prev_class][cur_class]/60

    if weather_state==1:
        Mn=w_rho*sep
    else:
        Mn=sep

    SD=math.sqrt(Mn**2/k)

    #z=int(random.random()*10000)

    t2=tm+serv_time*SD+Mn #Look up the stored N(0,1) value, serv_time, and convert it to the correct scale

    if t1<t2:
        straight_into_service=0
        t_out=t2
    else:
        straight_into_service=1
        t_out=t1

    return t_out,straight_into_service



def Normal_Conditional_GetServ(trav_time,rel_time,sv_time,prev_class,cur_class,tm,weather_state, Time_Sep: List[List[int]], w_rho: float, k: int, norm_cdf):

    #This is for the AC currently in service

    t1=rel_time+trav_time

    #Now do the service time
    sep=Time_Sep[prev_class][cur_class]/60

    if weather_state==1:
        Mn=w_rho*sep
    else:
        Mn=sep

    SD=math.sqrt(Mn**2/k)

    t2=truncnorm(sv_time+Mn,SD,tm,norm_cdf)

    #max_t=max(t1,t2)

    if t1<t2:
        straight_into_service=0
        t_out=t2
    else:
        straight_into_service=1
        t_out=t1

    return t_out,straight_into_service


def ncdf(x): #Standard Normal dist CDF
    return 0.5+0.5*math.erf(x/(2**0.5))

def npdf(x): #Standard Normal dist PDF
    return math.exp(-(x**2/2))/(math.sqrt(2*math.pi))

def truncnorm(Mn,SD,min_x,norm_cdf): #Generate a value x from the truncated normal distribution
    Min_U=math.ceil(ncdf((min_x-Mn)/SD)*1000)
    Max_U=1000 #int(ncdf((max_x-Mn)/SD)*10000)
    U1=int(Min_U+int(random.random()*(Max_U-Min_U+1)))
    x=norm_cdf[U1]*SD+Mn
    return x

def truncexp(Mn,SD,min_x): #Mean of the truncated normal distribution
    phi1=npdf((min_x-Mn)/SD)
    #phi2=npdf((max_x-Mn)/SD)
    phi2=1
    Phi1=ncdf((min_x-Mn)/SD)
    #Phi2=ncdf((max_x-Mn)/SD)
    Phi2=1
    if Phi1==1:
        texp=min_x
    else:
        texp=Mn+SD*((phi1-phi2)/(Phi2-Phi1))
    return texp
