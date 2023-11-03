#This is set up to run numerical experiments with a preset value of wiener_sig (also used as the weather variance parameter) and randomly-
#-generated values of k, thres1, pax weights, start & end times for bad weather (& need to also include random lam1 & lam2)

from __future__ import print_function, division
import sys
import math
import logging
import random
import csv
import time
import os

from stoch_runway_scheduler import weather, Genetic, Genetic_determ, Populate, Repopulate_VNS, sample_cond_gamma, getcost, Annealing_Cost, Perm_Heur, Perm_Heur_New, Calculate_FCFS, sample_gamma, gamma_create_cdf, Posthoc_Check, Update_Stats, Update_ETAs, Serv_Completions

# JF: these three options seem to be for logging purposes
#     they are broken for now as a separate output stream is created for each one, and all these currently
#     need to be passed a function where they are used.
stepthrough=0
step_summ=0
step_new=0

# Logs can similarly be set up for step_summ and step_new
stepthrough_logger = logging.getLogger('stepthrough')
stepthrough_logger.setLevel(level=logging.ERROR)
c_handler = logging.FileHandler("SRSP_stepthrough.log", mode='w')
# c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
# c_handler.setFormatter(c_format)
stepthrough_logger.addHandler(c_handler)

step_summ_logger = logging.getLogger('step_summ')
step_new_logger = logging.getLogger('step_new')

DATA_DIR = '/home/jamie/Insync/fairbrot@lancaster.ac.uk/OneDrive Biz - Shared/SRSP data files'

# JF: policy is an scheduling policy algorithm - may not affect anything now
# Alternate means sway between VNS and VNSD
Policy='Alternate' #FCFS, Perm, SA, GA or Alternate

Use_VNS=1
Use_VNSD=1
# Use_FCFS=1

Policies=[]
if Use_VNS==1:
    Policies.append('VNS')
if Use_VNSD==1:
    Policies.append('VNSD')
# if Use_FCFS==1:
# 	Policies.append('FCFS')

NoA=700 # number of aircraft - temporary value which will get changed later
S=40 # number of time slots (Rob not sure whether this is needed)
#thres1=15 #will get changed later
#thres2=0 #will get changed later

# fixed_elap=0.00015 #fixed CPU time for one "elap" for VNS algorithm
# fixed_elap_vnsd=0.00005 #fixed CPU time for one "elap" for VNSD algorithm
# fixed_repop_elap=0.0000825 #fixed CPU time for one "repop_elap"
# fixed_pop_elap=0.0003 #fixed CPU time for one "pop_elap"

# NoA=15
# S=2
# thres=0

max_FCFS=NoA #int(NoA/2)
conv_factor=1 #1 #0.1 #3 # no. of seconds in a minute for conversion purposes
norm_approx_min=100 #100 - Erlang can be approximated well for large k by Normal
Max_LookAhead=15 #15 #15 #30 #5 #30 #NoA #This is the length of a sequence, equivalent to parameter l in paper  - in paper this is 15
w_rho=10/9 #separation multiplier for bad weather (multiplies mean, not rate, so should be >1); we have used a reduction of 10% due to bad weather based on Odoni et al (2011) cited in Shone et al (2021)

#max_FCFS=0

f=open("SRSP_out_%s_%s.csv" % (Policy, conv_factor), "w") # detailed results - for each rep what happened for each aircraft
g=open("AC_predictions_%s_%s.csv" % (Policy, conv_factor), "w")
#h=open("cost_log_%s_%s.csv" % (Policy, conv_factor), "w")
gg=open("SRSP_rep_results_%s_%s.csv" % (Policy, conv_factor), "w")  # summaries for each rep

f1=open("Perm_Heur_out_%s_%s.csv" % (Policy, conv_factor), "w") # output from a particular heuristic - may not be needed (can remove if you like)
#f2=open("Detailed_Out_%s_%s.csv" % (Policy, conv_factor), "w")

#rr=open("Runtime_Tracker_%s_%s.csv" % (Policy, conv_factor), "w")

#f5=open("Wiener_Test_%s_%s.csv" % (Policy, conv_factor), "w")

# if stepthrough==1:
#     # JF: CSV extension is used, but prints don't always follow this format
#     st=open("SRSP_stepthrough.csv", "w")
# if step_summ==1:
#     st2=open("SRSP_step_summarised.csv", "w")
# if step_new==1:
#     st3=open("SRSP_step_new.csv", "w")

# not needed anymore
st4=open("elap_out.csv", "w")

f.write('Policy'+',''Rep'+','+'AC'+','+'Flight Num'+','+'Prev Class'+','+'Cur Class'+','+'Time Sep'+','+'Orig PS time'+','+'PS time'+','+'Pool Arrival'+','+'Release Time'+','+'Travel Time'+','+'Weather Coeff'+','+'Enters Serv'+','+'Actual Serv'+','+'Ends Serv'+','+'Lateness'+','+'Queue Delay'+','+'Pax Weight'+','+'Cost'+','+'counter'+','+'qp'+','+'Predicted total'+'\n')








#random.seed(100)

wiener_sig=0.1 #0.1 #1 #0.1 #standard deviation for Brownian motion
weather_sig=wiener_sig #this assumption is being made in the paper for simplicity

#Import the Wiener cdf - JF: could simplify
print('*** Importing the Wiener array...')
#wiener_cdf=[[0]*(1000) for _ in range(12000)]
wiener_cdf=[[0]*(1000) for _ in range(12000)] #This array is used to store quantiles of the Inverse Gaussian distribution (see equation (7) in 1st revision of paper). Rows represent different values for the difference between current time and ETA. Columns are for different quantiles (from 0 to 1 in increments of 0.001).
weather_cdf=[[0]*(1000) for _ in range(12000)] #This array will be identical to wiener_cdf for simplicity.
if wiener_sig==0.1:
    with open(os.path.join(DATA_DIR, 'wiener_array_sig0point1.csv'), 'r') as csvfile: #this data file contains pre-generated quantiles of the Inverse Gaussian distribution with sigma=0.1 (Rob has Python scripts to generate these)
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata=list(datareader)
        for i in range(12000):
            for j in range(1000):
                wiener_cdf[i][j]=float(inputdata[i][j])
                weather_cdf[i][j]=float(inputdata[i][j])
elif wiener_sig==0.3:
    with open(os.path.join(DATA_DIR, 'wiener_array_sig0point3.csv'), 'r') as csvfile: #Inverse Gaussian with sigma=0.3
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata=list(datareader)
        for i in range(12000):
            for j in range(1000):
                wiener_cdf[i][j]=float(inputdata[i][j])
                weather_cdf[i][j]=float(inputdata[i][j])
elif wiener_sig==0.5:
    with open(os.path.join(DATA_DIR, 'wiener_array_sig0point5.csv'), 'r') as csvfile: #Inverse Gaussian with sigma=0.5
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata=list(datareader)
        for i in range(12000):
            for j in range(1000):
                wiener_cdf[i][j]=float(inputdata[i][j])
                weather_cdf[i][j]=float(inputdata[i][j])
elif wiener_sig==0.7:
    with open(os.path.join(DATA_DIR, 'wiener_array_sig0point7.csv'), 'r') as csvfile: #Inverse Gaussian with sigma=0.7
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata=list(datareader)
        for i in range(12000):
            for j in range(1000):
                wiener_cdf[i][j]=float(inputdata[i][j])
                weather_cdf[i][j]=float(inputdata[i][j])
elif wiener_sig==0.9:
    with open(os.path.join(DATA_DIR, 'wiener_array_sig0point9.csv'), 'r') as csvfile: #Inverse Gaussian with sigma=0.9
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata=list(datareader)
        for i in range(12000):
            for j in range(1000):
                wiener_cdf[i][j]=float(inputdata[i][j])
                weather_cdf[i][j]=float(inputdata[i][j])
else:
    assert 1==2 # JF: raise an exception instead

# #Import the weather cdf
# print('*** Importing the weather array...')
# weather_cdf=[[0]*(1000) for _ in range(12000)]
# with open('wiener_array_sig0point1.csv', 'r') as csvfile:
# 	datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
# 	inputdata=list(datareader)
# 	for i in range(12000):
# 		for j in range(1000):
# 			weather_cdf[i][j]=float(inputdata[i][j])

#Import the normal cdf
normcdf=[0]*(10001)
with open(os.path.join(DATA_DIR, 'norm_cdf.csv'), 'r') as csvfile: #data file contains quantiles from the standard normal distribution
    datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
    inputdata=list(datareader)
    for i in range(10001):
        normcdf[i]=float(inputdata[i][0])

#Set the parameters

NoC=4 #no. of aircraft classes

Ac_finished=[0]*NoA # 'finished' means aircraft has completed landing, i.e. service time has finished
pred_serv=[0]*NoA #seems to be not used anymore

tau=30 #Determines when an aircraft is deemed to have entered the pool. E.g. if tau=30, aircraft enters pool when its ETA is 30 minutes away from the current time.
t=15 #length of a time slot in minutes

print('wiener_sig: '+str(wiener_sig))

pool_max=6 #Used as a parameter for "perm heuristic" which searches for the best landing sequence under perfect information, i.e. assumes all random information already known
list_min=6 #Also used only for the "perm heuristic""

GA_PopSize=20 #Initial number of sequences in the population, written as S in paper (see Section 3.1)
Ov_GA_counter=0 #only used for stepthrough purposes; to do with counting how many times the 'Genetic' function has been called
Tabu_Size=50 #only used by a Tabu search heuristic which we're not using anymore
VNS_limit=25 #important parameter, determines how many non-improving heuristic moves need to be made before a mutation is carried out; this is written as m_{mut} in paper (see the flow chart, Figure 3)
VNS_counter=0 #this counts how many non-improving heuristic moves we've made since the last reset
tot_mut=0 #counts how many total mutations we've made; really just for output purposes

AC_List_Length=6 #doesn't seem to be used anymore
perm_length=4 #doesn't seem to be used anymore

# print('Testing the wiener array')
# sched=240
# for j in range(1000):
# 	#First, use the wiener array
# 	z=int(random.randrange(1,999))
# 	trav=wiener_cdf[10*sched][z]
# 	f5.write(str(trav)+',')
# 	#Next, manually generate the trajectory
# 	i=0
# 	dt=0.01
# 	ETA=sched
# 	while 1==1:
# 		i+=dt
# 		ETA=random.gauss(ETA,0.1*wiener_sig)
# 		if i>=ETA:
# 			trav=i
# 			break
# 	f5.write(str(trav)+'\n')

# Dis_Sep=[[2,3,4],[3,4,5],[4,5,6]] #Distance separation requirements in miles - not sued anymore
#Time_Sep=[[96,138,240],[60,72,162],[60,72,102],[60,72,102]] #Time separations in seconds taken from Solak et al (2018) appendix; the 4th array is for the situation where there is no leading aircraft
# sep also includes separation when there is not "leading" aircraft - maybe last row should be 0s
Time_Sep=[[97,121,121,145],[72,72,72,97,97],[72,72,72,72],[72,72,72,72],[72,72,72,72]] #Time separations in seconds taken from Bennell et al (2017) with H, U, M, S as the 4 classes; the 5th array is for the situation where there is no leading aircraft
# JF: Time_Sep is List[List[int]]

Ac_Info=[0]*NoA #this will be a multi-dimensional list storing lots of information about each aircraft; see later
Ac_class=[0]*NoA #this will store the weight class for each aircraft
Arr_Ps=[0]*NoA #this stores the adjusted scheduled times for aircraft after applying the random pre-tactical delay
Dep_Ps=[0]*NoA #not needed anymore because we do not consider departures
Orig_Ps=[0]*NoA #original pre-scheduled times of aircraft, before applying the pre-tactical delay
flight_id=[0]*NoA #stores the aircraft flight numbers, for identification purposes
pax_weight=[0]*NoA #stores the randomly-generated cost weightings for aircraft, based on (hypothetical) numbers of passengers carried; written as g_i in the paper (see objective function (13))

no_reps=10000 #total number of random scenarios that we will simulate; in each scenario we evaluate the performances of different algorithms such as SimHeur, DetHeur, FCFS

# if Policy=='Alternate':
# 	SubPolicy='Perm'
# else:
# 	SubPolicy=Policy

rep=0 #counter of which scenario we're currently on
policy_index=0 #indicates which policy we're currently evaluating, e.g. SimHeur, DetHeur etc (if this is zero then we take the first policy from the list of policies to be evaluated)

#for rep in range(no_reps):
while rep<no_reps:

    repn=rep #int(rep/100+1)
    random.seed(repn*100) #set the random seed for generating random parameter values; seed in set according to the replication (scenario) number

    #Import the flight data
    print('*** Importing the flight data...')
    #When reading in the flight data from the data file we only want to include flights with a pre-scheduled time between 6AM (360 mins) and 2PM (840 mins), including 6AM but not including 2PM
    min_ps_time=360 #inclusive
    max_ps_time=840 #non-inclusive

    AC=0 #counts how many flights have been read in so far
    #earliest_ps_time=0
    with open(os.path.join(DATA_DIR, 'flight_pretac_data.csv'), 'r') as csvfile: #data file includes the on-time performance data such as means, variance of lateness based on 1 year of historical info from FlightRadar [DON'T CHANGE THIS FILE]
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        inputdata=list(datareader)
        for i in range(1,688): #start from 1 because there's a title row
            ps_time=float(inputdata[i][1]) #pre-scheduled time
            # ft_time=float(inputdata[i][5])
            # if i==0 or ps_time<earliest_ps_time:
            # 	earlest_ps_time=ps_time
            if ps_time>=min_ps_time and ps_time<max_ps_time:

                arr_time=int(inputdata[i][1]) #initially we set this equal to the pre-scheduled time but it will get adjusted later by applying a random pre-tactical delay
                dep_time=int(inputdata[i][4]) #departure time from the origin airport
                flight_name=str(inputdata[i][3]) #flight number
                sched_time=int(inputdata[i][5]) #scheduled flight time, i.e. difference between scheduled departure and arrival time
                lateness_mn=float(inputdata[i][6]) #mean lateness based on historical data
                lateness_var=float(inputdata[i][7]) #variance of lateness based on historical data

                Ac_class[AC]=int(inputdata[i][2]) #weight class
                Orig_Ps[AC]=arr_time #original pre-scheduled time (as opposed to scheduled time following pre-tactical delay)

                #The equations for xi_bar, si2, h_i, alpha and beta below are for calculating the parameters of the gamma distribution used for the pre-tactical delay. Details of this method are in Section 4 of the paper.

                xibar=arr_time+lateness_mn
                si2=lateness_var
                h_i=dep_time-15

                alpha=((xibar-h_i)**2)/(si2-(wiener_sig**2)*(xibar-h_i))
                beta=(xibar-h_i)/(si2-(wiener_sig**2)*(xibar-h_i))

                if alpha>0 and beta>0:
                    pretac_delay=sample_gamma(alpha,1/beta)-(arr_time-h_i) #Here we sample from a gamma distribution to get the pre-tactical delay for the flight under consideration.
                else:
                    pretac_delay=lateness_mn #In this case the pre-tactical delay is set equal to the average lateness rather than being sampled randomly.

                Arr_Ps[AC]=arr_time+pretac_delay #Stores the initial ETA of the current aircraft, after adjustment based on pre-tactical delay
                Dep_Ps[AC]=h_i #This is actually 15 minutes before the flight's scheduled departure time, and indicates the point at which we assume the ETA starts varying according to Brownian motion (see paper Section 2)
                flight_id[AC]=flight_name #Store the flight number

                AC+=1
            elif ps_time>=max_ps_time: #Indicates that we have got to the end of the set of flights scheduled to arrive by 2PM
                break

    print('No. of ACs: '+str(AC))
    NoA=AC
    #NoA=8

    for i in range(NoA): #HERE WE RE-SCALE TIME SO THAT TIME '6AM' IS COUNTED AS TIME (ZERO+60). WE START SIMULATING FROM TIME ZERO, I.E. AN HOUR BEFORE 6AM.
        Arr_Ps[i]+=-min_ps_time+60
        Dep_Ps[i]+=-min_ps_time+60
        Orig_Ps[i]+=-min_ps_time+60

    for i in range(NoA):
        #Ac_class[i]=int(random.random()*3)
        #Ac_class[i]=0
        #z=tau/(S*t)
        #Arr_Ps[i]=(z+(1-z)*random.random())*(S*t)
        #Passenger weight
        if Ac_class[i]==0:
            pax_weight[i]=0.2*random.random()+0.8 #Flights in the 'heavy' class have a passenger weight between 0.8 and 1
        elif Ac_class[i]==1 or Ac_class[i]==2:
            pax_weight[i]=0.2*random.random()+0.6 #Flights in the 'upper medium' or 'lower medium' class have a passenger weight between 0.6 and 0.8 
        else:
            pax_weight[i]=0.2*random.random()+0.4 #Flights in the 'small' class have a passenger weight between 0.4 and 0.6

    SubPolicy=Policies[policy_index] #SubPolicy indicates the policy we are currently considering (e.g. SimHeur, DetHeur)

    if SubPolicy=='GA' or SubPolicy=='VNS':
        GA_LoopSize=500 # This is
    elif SubPolicy=='GAD' or SubPolicy=='VNSD':
        GA_LoopSize=1

    #random.seed(100)
    #p=rep/no_reps

    # #Randomly generate the Erlang service time parameter
    # z=int(random.random()*5)+1
    # z*=0.05 #z is either 0.05, 0.1, 0.15, 0.2 or 0.25; this corresponds to the coefficient of variation for service times
    # k=int(1/(z**2))

    z=random.random()*5 #*5
    if z<1:
        k=16
    elif z<2:
        k=25
    elif z<3:
        k=44
    elif z<4:
        k=100
    else:
        k=400

    #k=50

    print('k: '+str(k))

    #Randomly generate lam1 and lam2
    z=random.random()*5 #*5
    if z<1:
        lam1=0.1
    elif z<2:
        lam1=0.3
    elif z<3:
        lam1=0.5
    elif z<4:
        lam1=0.7
    else:
        lam1=0.9

    lam1=0.1
    lam2=1-lam1

    print('lam1: '+str(lam1)+' lam2: '+str(lam2))

    #Randomly generate thres1
    z=random.random()
    if z<0.5:
        thres1=0
    else:
        thres1=15
    thres2=0

    if k>=norm_approx_min:
        NormalApprox=1
        gamma_cdf=[]
    else:
        NormalApprox=0
        print('*** Creating the gamma CDF...')
        gamma_cdf=gamma_create_cdf(k)

    # if NormalApprox==0:
    # 	print('*** Generating the conditional serv time arrays...')
    # 	Serv_Pr=[[0]*(NoC) for _ in range(NoC+1)]
    # 	last_epoch=[[0]*(NoC) for _ in range(NoC+1)]
    # 	IFR_Serv_Pr=[[0]*(NoC) for _ in range(NoC+1)]
    # 	IFR_last_epoch=[[0]*(NoC) for _ in range(NoC+1)]
    # 	for i in range(NoC+1):
    # 		for j in range(NoC):
    # 			Serv_Pr[i][j],last_epoch[i][j]=conditional_phase_serv(S*t,0,k,i,j,1,Time_Sep)
    # 			IFR_Serv_Pr[i][j],IFR_last_epoch[i][j]=conditional_phase_serv(S*t,0,k,i,j,1/w_rho,Time_Sep)
    # 	print('last_epoch: '+str(last_epoch))
    # 	print('IFR_last_epoch: '+str(IFR_last_epoch))
    # else:
    # 	Serv_Pr=[]
    # 	IFR_Serv_Pr=[]
    # 	last_epoch=[]
    # 	IFR_last_epoch=[]

    print('*** Generating initial aircraft info...')
    Ac_Info=[0]*NoA

    for i in range(NoA):
        Status=0 #0: not ready yet (arrival), 1: in arrival pool, 2: added to arrival queue, 3: not ready yet (departure), 4: in departure pool, 5: added to departure queue, 6: finished.

        z=int(random.random()*1001)

        if NormalApprox==0:
            z=int(random.random()*1000)
            ServPercs=gamma_cdf[z] #ServPercs is now the service time sampled from a standard Gamma(k,1) dist
        else:
            z=int(random.random()*10000)
            ServPercs=normcdf[z] #ServPercs is now the service time sampled from a standard N(0,1) dist

        # if NormalApprox==0:
        # 	ServPercs=[0]*k
        # 	for j in range(k):
        # 		ServPercs[j]=random.random()
        # else:
        # 	ServPercs=random.random()

        #print('i: '+str(i)+' TravTime: '+str(TravTime))
        #Arr_Ps=(1-random.random()*random.random())*(S*t)
        Ac_Info[i]=[Status,Ac_class[i],Arr_Ps[i],Arr_Ps[i],0,0,0,ServPercs,0,0,pax_weight[i],0,1,0,0,0,0,Dep_Ps[i],Orig_Ps[i],flight_id[i]] #index 2 is pre-scheduled time (plus pre-tactical delay), index 3 is latest ETA, index 4 is the time at which aircraft is released from pool, index 5 is the time at which aircraft enters service, index 6 is the travel time (generated in advance), index 7 is the list of random numbers used for service phase completions, index 8 is the actual service time s1+Z2 (worked out after class information is known, index 9 is the actual time that they join the pool (generated in advance)), index 10 is the passenger weight, index 11 is an indicator to show whether or not the AC's travel time has already been completed, index 12 is the weather state at the time of release, index 13 is the counter (for Perm only), index 14 is qp (Perm only), index 15 is predicted total cost at time of release, index 16 is the actual service completion time, index 17 is the scheduled departure time, index 18 is the original pre-scheduled arrival time before adding pre-tactical delay, index 19 is the flight number

    Ac_Info.sort(key=lambda x: x[2])

    print('*** Generating the ETA trajectory array...')
    Brown_Motion=[[0]*int(S*t*2*100) for _ in range(NoA)]
    for i in range(NoA):
        j=0
        Dep_time=Ac_Info[i][17]
        Ps_time=Ac_Info[i][2]
        if Dep_time<0:
            ETA=random.gauss(Ps_time,math.sqrt(0-Dep_time)*wiener_sig) #Update the latest ETAs for ACs that already had their dep time before time zero
        else:
            ETA=Ac_Info[i][2] #ETA = pre-scheduled time
        Brown_Motion[i][0]=ETA
        if 0>=ETA-tau:
            Ac_Info[i][9]=0 #index 9 is actual pool arrival
            chk=1
            ETA=tau
            threshold_time=0
        else:
            chk=0
        while j<S*t*1.5*100:
            j+=1 #step forward in hundredths of a minute
            if j>Ac_Info[i][17]: #only update ETA if we've gone beyond the AC's departure time
                ETA=random.gauss(ETA,0.1*wiener_sig)
            # if i==6:
            # 	print('j: '+str(j)+' ETA: '+str(ETA))
            Brown_Motion[i][j]=ETA
            if j/100>=ETA-tau and chk==0:
                threshold_time=round(j/100,2) #j/100 is the 'current time'
                Ac_Info[i][9]=threshold_time
                chk=1
            elif j/100>=ETA and chk==1:
                runway_time=round(j/100,2)
                Ac_Info[i][6]=runway_time-threshold_time #index 6 is travel time
                chk=2
                break

    stepthrough_logger.info('AC'+','+'Class'+','+'PS time'+','+'Pool arrival'+','+'Travel time'+','+'Runway time'+'\n')
    for AC in range(NoA):
        Ac_Infoi=Ac_Info[AC]
        stepthrough_logger.info(str(AC)+','+str(Ac_Infoi[1])+','+str(Ac_Infoi[2])+','+str(Ac_Infoi[9])+','+str(Ac_Infoi[6])+','+str(Ac_Infoi[9]+Ac_Infoi[6])+'\n')
    stepthrough_logger.info('\n')

    print('*** Generating the weather transition array...')
    weather_lb=[0]*(int(S*t*8*100)+1)
    weather_ub=[0]*(int(S*t*8*100)+1)

    #4 possible cases: no bad weather, 30 minutes of bad weather, 60 minutes of bad weather or 120 minutes of bad weather (bad weather is always forecast for the middle of the day)
    z=random.random()*4
    #z=0
    if z<1:
        wlb=0 #120 #Forecast for bad weather start at the beginning of day (can change this to non-zero)
        wub=0 #180 #Forecast for bad weather end at the beginning of day (can change this to non-zero)
    elif z<2:
        wlb=285
        wub=315 #Note: 300 is 10AM according to the rescaling of time used earlier
    elif z<3:
        wlb=270
        wub=330
    else:
        wlb=240
        wub=360

    wlb_tm=0 #Actual (randomly generated) time at which bad weather starts; leave this as zero
    wub_tm=0 #Actual (randomly generated) time at which bad weather ends; leave this as zero
    weather_lb[0]=wlb
    weather_ub[0]=wub
    old_lb=wlb
    old_ub=wub

    j=0
    while j<S*t*2*100:
        j+=1
        new_lb=random.gauss(old_lb,0.1*weather_sig) #random.gauss(old_lb,0.01)
        new_ub=random.gauss(old_ub,0.1*weather_sig)
        if new_lb>new_ub:
            new_ub=new_lb
        if j/100>=new_lb and wlb_tm==0:
            wlb_tm=j/100
        if j/100>=new_ub and wub_tm==0:
            wub_tm=j/100
        weather_lb[j]=new_lb
        weather_ub[j]=new_ub
        old_lb=new_lb
        old_ub=new_ub

    stepthrough_logger.info('wlb_tm:'+','+str(wlb_tm)+','+'wub_tm'+','+str(wub_tm)+'\n'+'\n')

    # if SubPolicy=='SA':
    # 	Anneal_Seq=[0]*NoA
    # 	for i in range(NoA):
    # 		Anneal_Seq[i]=i

    if SubPolicy in ('GA','GAD','VNS','VNSD'):
        #Generate the initial population of sequences
        print('Generating initial population of sequences...')

        no_ACs=min(Max_LookAhead,NoA)
        FSFS_seq=[0]*no_ACs
        for i in range(no_ACs):
            FSFS_seq[i]=i

        GA_PopList,GA_Info=Populate(Ac_Info, FSFS_seq,[],FSFS_seq,GA_PopSize,Max_LookAhead, stepthrough, step_summ, step_new)
        # print('GA_PopList: '+str(GA_PopList))
        # print('GA_Info: '+str(GA_Info))

        Opt_Seq=FSFS_seq[:]
        OptCost=1000000
        queue_probs=[0]*NoA
        Opt_List=[]
        Opt_Seqs=[]
        # Opt_List.append([FSFS_seq[:],0,0,queue_probs,0])
        # Opt_Seqs=[FSFS_seq[:]]
        Opt_Size=10 #10 #3 #10

        while len(Opt_List)<Opt_Size:
            new_seq=FSFS_seq[:]
            random.shuffle(new_seq)
            if new_seq not in GA_PopList and new_seq not in Opt_Seqs:
                Opt_List.append([new_seq[:],0,0,queue_probs,0])
                Opt_Seqs.append(new_seq[:])

        GA_Check_Increment=GA_LoopSize/10

        GA_counter=0
        GA_CheckSize=GA_Check_Increment
        print('GA_CheckSize: '+str(GA_CheckSize))
        Ov_GA_counter=0

        stepthrough_logger.info('Initial population of sequences:'+'\n')
        for j in range(len(GA_PopList)):
            stepthrough_logger.info(str(j)+','+str(GA_PopList[j])+'\n')
        stepthrough_logger.info('\n')

        print('Ov_GA_counter: '+str(Ov_GA_counter))

    for AC in range(NoA):
        print('AC: '+str(AC)+' Class: '+str(Ac_Info[AC][1])+' Orig Ps time: '+str(Ac_Info[AC][18])+' Ps time: '+str(Ac_Info[AC][2])+' Arrives in pool: '+str(Ac_Info[AC][9])+' Travel time: '+str(Ac_Info[AC][6]))

    Ac_queue=[]
    Left_queue=[]
    Arr_Pool=[]
    Dep_Pool=[]
    Arr_NotReady=[]
    Dep_NotReady=[]
    Old_Perms=[]
    Old_Perm_Info=[]

    tm=0 #-30 #earliest_ps_time-30 #-30 #time counter
    old_tm=tm
    totserv=0 #counter of number of aircraft served
    latest_class=4 #class of latest aircraft to be added to the queue, initially set to 3
    prev_class=4 #class of previous aircraft to be served, initially set to 3
    Ac_added=[]
    counter=0
    qp=0
    weather_state=0 #0 means it's good weather, 1 means bad weather, 2 means good weather again
    real_queue_complete=0
    next_completion_time=0
    Pop_elap=0
    max_d=1
    pruned=0
    soln_evals_tot=0
    soln_evals_num=0
    tot_mut=0

    # elap_tot=0
    # elap_num=0
    # Repop_elap_tot=0
    # Repop_elap_num=0
    # Pop_elap_tot=0
    # Pop_elap_num=0

    tot_arr_cost=0
    tot_dep_cost=0

    Loop_Nums=0
    Loop_Evals=0

    Long_Perm_List=[]
    Perm_Info=[]

    #Generate the initial pool and notready arrays
    print('*** Generating the initial pool...')
    for i in range(NoA):
        Ac_Infoi=Ac_Info[i]
        if Ac_Infoi[3]-tau<=0:
            Arr_Pool.append(i)
            print('Aircraft '+str(i)+' initially included in pool (ETA is '+str(Ac_Infoi[3])+')')
            Ac_Infoi[0]+=1
        else:
            Arr_NotReady.append(i)

    #print('Ac_Info: '+str(Ac_Info))
    print('Arr_Pool: '+str(Arr_Pool))
    print('Ac_queue: '+str(Ac_queue))
    print('Arr_NotReady: '+str(Arr_NotReady))
    # for i in range(NoA):
    # 	print('i: '+str(i)+' Ac_Infoi[9]: '+str(Ac_Info[i][9]))

    print('*** Into main loop for rep '+str(rep)+' and policy '+str(SubPolicy)+'...')
    begin_time=time.time()
    mv_time=0
    if policy_index==0:
        gg.write(str(rep)+','+str(conv_factor)+','+str(wiener_sig)+','+str(k)+','+str(lam1)+','+str(lam2)+','+str(wlb)+','+str(wub)+','+str(wlb_tm)+','+str(wub_tm)+','+str(thres1)+','+str(thres2)+',')

    initial_time=time.time()

    while totserv<NoA:

        # if tm>=10:
        # 	stepthrough=0

        #print('tm: '+str(tm))

        #print('tm: '+str(tm)+' Best sequence is '+str(GA_Info[0][0])+' Estimated cost is '+str(GA_Info[0][2]))

        # latest_timer=time.time()

        stepthrough_logger.info('tm is '+','+str(tm)+'\n')
        stepthrough_logger.info('Arr_NotReady is '+','+str(Arr_NotReady)+'\n')
        stepthrough_logger.info('Arr_Pool is '+','+str(Arr_Pool)+'\n')
        stepthrough_logger.info('Ac_queue is '+','+str(Ac_queue)+'\n')
        stepthrough_logger.info('Left_queue is '+','+str(Left_queue)+'\n')
        stepthrough_logger.info('Ac_added is '+','+str(Ac_added)+'\n'+'\n')

        # if SubPolicy=='GAD' and tm>28.5:
        # 	GA_Info.sort(key=lambda x: x[2])
        # 	print('Time: '+str(tm)+' Opt sequence: '+str(GA_Info[0][0]))

        if tm>=0 and len(Ac_added)>0 and Ac_added[0] in Arr_Pool:
            if SubPolicy in ('GA','GAD','VNS','VNSD'):
                Opt_List.sort(key=lambda x: x[2])
                GA_Info.sort(key=lambda x: x[2])
                if len(Opt_List)>0 and len(GA_Info)>0:
                    if Opt_List[0][2]<GA_Info[0][2]:
                        base_seq=Opt_List[0][0][:]
                        pred_cost=Opt_List[0][2]
                    else:
                        base_seq=GA_Info[0][0][:]
                        pred_cost=GA_Info[0][2]
                elif len(Opt_List)>0:
                    base_seq=Opt_List[0][0][:]
                    pred_cost=Opt_List[0][2]
                elif len(GA_Info)>0:
                    base_seq=GA_Info[0][0][:]
                    pred_cost=GA_Info[0][2]
                else:
                    assert 1==2
            for AC in Ac_added:
                if AC in Arr_Pool:
                    Arr_Pool.remove(AC)
                    if SubPolicy=='SA':
                        Anneal_Seq.remove(AC)
                    Ac_queue.append(AC)
                    # if SubPolicy=='GAD':
                    # 	f.write('AC '+str(AC)+' added to queue, Ac_queue is '+str(Ac_queue)+', stored_queue_complete is '+str(stored_queue_complete)+'\n')
                    if SubPolicy=='Perm':
                        continue #print('Added AC '+str(AC)+' to the queue, counter is '+str(counter)+', LPLL is '+str(len(Long_Perm_List))+', qp is '+str(qp))
                    elif SubPolicy in ('GA','GAD','VNS','VNSD'):
                        #print('Added AC '+str(AC)+' to the queue, counter is '+str(Ov_GA_counter)+', qp is '+str(qp))
                        Ac_Info[AC][15]=pred_cost
                        stepthrough_logger.info('Added AC '+str(AC)+' to the queue, counter is '+str(Ov_GA_counter)+', qp is '+str(qp)+'\n')
                        step_summ_logger.info('Added AC '+str(AC)+' to the queue, counter is '+str(Ov_GA_counter)+', qp is '+str(qp)+'\n')
                        step_new_logger.info('Added AC '+str(AC)+' to the queue, counter is '+str(Ov_GA_counter)+', qp is '+str(qp)+'\n')
                        base_seq.remove(AC)

                    else:
                        print('Added AC '+str(AC)+' to the queue')

                    real_queue_complete,next_completion_time,latest_class,Ov_GA_counter=Update_Stats(tm,AC,Ac_Info,Ac_queue,real_queue_complete,wlb_tm,wub_tm,latest_class,Ov_GA_counter,next_completion_time, k, Time_Sep, norm_approx_min, w_rho, SubPolicy, counter, qp)

                else:
                    break

            if SubPolicy in ('GA','GAD','VNS','VNSD'):
                # if tm>=0 and int(tm*100)!=int(old_tm*100):
                # 	rr.write('Populate'+',')
                GA_PopList,GA_Info=Populate(Ac_Info, base_seq,Arr_Pool,Arr_NotReady,GA_PopSize,Max_LookAhead,stepthrough, step_summ, step_new)
                queue_probs=[0]*(len(Arr_Pool)+len(Arr_NotReady))
                # Pop_elap_tot+=Pop_elap 
                # Pop_elap_num+=1

                Opt_List=[]
                Opt_Seqs=[]
                # Opt_List.append([base_seq[:],0,0,queue_probs,0])
                # Opt_Seqs.append(base_seq[:])
                c=0
                while len(Opt_List)<Opt_Size and c<25:
                    new_seq=base_seq[:]
                    random.shuffle(new_seq)
                    if new_seq not in GA_PopList and new_seq not in Opt_Seqs:
                        Opt_List.append([new_seq[:],0,0,queue_probs,0])
                        Opt_Seqs.append(new_seq[:])
                        c=0
                    else:
                        c+=1

                # Opt_List=[]
                # Opt_List.append([base_seq[:],0,0,queue_probs,0])
                GA_counter=0
                GA_CheckSize=GA_Check_Increment
            else:
                Pop_elap=0

            Long_Perm_List=[]
            Perm_Info=[]

        if SubPolicy in ('GA','GAD','VNS','VNSD'):
            #print('tm: '+str(tm)+' GA_counter: '+str(GA_counter))
            Repop_elap=0
            if len(Arr_Pool)+len(Arr_NotReady)>0: #4:
                if GA_counter>=GA_LoopSize or pruned==1: #or max_d<0.01:
                    #print('tm: '+str(tm)+' Best sequence is '+str(GA_Info[0][0])+' Estimated cost is '+str(GA_Info[0][2]))
                    #print('tm: '+str(tm)+' GA_counter: '+str(GA_counter)+', Repop')
                    # if SubPolicy=='GA' or SubPolicy=='GAD':
                    # 	GA_PopList,GA_Info,Repop_elap,Tabu_List,Opt_Seq,OptCost,Opt_List=Repopulate(GA_PopList,GA_Info,Arr_Pool,Arr_NotReady,GA_PopSize,Tabu_List,Tabu_Size,Opt_Seq,OptCost,Opt_List,Opt_Size,stepthrough, step_summ, step_new)
                    # else:
                    # if tm>=0 and int(tm*100)!=int(old_tm*100):
                    # 	rr.write('Repopulate_VNS'+',')
                    Loop_Nums+=1
                    #print('Loop_Nums: '+str(Loop_Nums))
                    Loop_Evals+=GA_counter
                    GA_PopList,GA_Info,Opt_Seq,OptCost,Opt_List,VNS_counter,tot_mut=Repopulate_VNS(GA_PopList,GA_Info,Arr_Pool,Arr_NotReady,GA_PopSize,Opt_Seq,OptCost,Opt_List,Opt_Size,Max_LookAhead,VNS_counter,VNS_limit,tot_mut, stepthrough, step_summ, step_new)
                    GA_counter=0
                    GA_CheckSize=GA_Check_Increment
                    mv_time=1
                    if SubPolicy=='GA':
                        max_d=1
                    # Repop_elap_tot+=Repop_elap
                    # Repop_elap_num+=1
            else:
                mv_time=1
        else:
            Repop_elap=0

        #First, get the AC List

        if tm>=0 and tm<=wub_tm and int(tm*100)!=int(old_tm*100):
            #Permute and update the weather
            wlb=weather_lb[int(tm*100)] #random.gauss(wlb,0.05)
            wub=weather_ub[int(tm*100)]
            #rr.write('\n'+str(SubPolicy)+','+str(tm)+','+str(time.time()-begin_time)+',')

        if len(Arr_Pool)+len(Arr_NotReady)>0:
            # if SubPolicy=='FCFS':
            # 	Ac_added,elap=FCFS_rule(Ac_Info,Arr_Pool,Dep_Pool)
            if SubPolicy=='VNS':
                #print('soln_evals_tot: '+str(soln_evals_tot))
                # if tm>=0 and int(tm*100)!=int(old_tm*100):
                # 	rr.write('Genetic'+',')
                Ac_added,counter,qp,max_d,pruned,GA_CheckSize,GA_counter,soln_evals_tot,soln_evals_num=Genetic(Ac_Info,Arr_Pool,Arr_NotReady,Ac_queue,Left_queue,max(tm,0),NoA,k,prev_class,GA_PopList,GA_Info,wiener_cdf,GA_LoopSize,GA_CheckSize,GA_counter,tot_arr_cost+tot_dep_cost,wlb,wub,weather_cdf,Opt_List,max_d,soln_evals_tot,soln_evals_num,gamma_cdf,normcdf, norm_approx_min, tau, Max_LookAhead, Time_Sep, thres1, thres2, lam1, lam2, GA_Check_Increment, Opt_Size, w_rho, stepthrough)
                Ov_GA_counter+=1
                #GA_counter+=1
                stepthrough_logger.info('GA_counter is '+','+str(GA_counter)+'\n')
                #Tabu_List.sort(key=lambda x: x[2])
                #print('tm: '+str(tm)+' Opt_Seq: '+str(Tabu_List[0][0])+' Cost: '+str(Tabu_List[0][2]))
                #print('tm: '+str(tm)+' Opt_Seq: '+str(Opt_Seq)+' Cost: '+str(OptCost))
            elif SubPolicy=='VNSD':
                Ac_added,counter,qp,stored_queue_complete=Genetic_determ(Ac_Info,Arr_Pool,Arr_NotReady,Ac_queue,Left_queue,max(tm,0),NoA,k,prev_class,GA_PopList,GA_Info,wiener_cdf,wlb,wub,Opt_List, norm_approx_min, tau, Max_LookAhead, Time_Sep, thres1, thres2, lam1, lam2, tot_arr_cost, tot_dep_cost, w_rho, stepthrough, step_summ, step_new)
                Ov_GA_counter+=1
                GA_counter+=1
                stepthrough_logger.info('GA_counter is '+','+str(GA_counter)+'\n')
                #Tabu_List.sort(key=lambda x: x[2])
                #print('tm: '+str(tm)+' Opt_Seq: '+str(Tabu_List[0][0])+' Cost: '+str(Tabu_List[0][2]))
                #print('tm: '+str(tm)+' Opt_Seq: '+str(Opt_Seq)+' Cost: '+str(OptCost))
            # elap_tot+=elap
            # elap_num+=1
            #st4.write(str(tm)+','+str(elap)+'\n')

        else:
            Ac_added,elap,counter,qp=[],0.1,0,0

        # if len(Arr_Pool)+len(Arr_NotReady)==18:
        # 	print('tm: '+str(tm)+' Arr_Pool: '+str(Arr_Pool)+' Opt seq: '+str(GA_Info[0][0])+' Ac_added: '+str(Ac_added)+' counter: '+str(counter)+' qp: '+str(qp))

        latest_time=(time.time()-initial_time)/conv_factor

        # if tm>=0:

        #start_update_time=time.time()

        if int(tm*100)!=int(old_tm*100):
            Update_ETAs(Ac_Info,Arr_NotReady,Dep_NotReady,Ac_queue,tm,Brown_Motion, Arr_Pool, tau)

        if len(Ac_queue)>0 and tm>=next_completion_time: #len(Ac_queue)>0:
            arr_cost,dep_cost,totserv,prev_class,Ac_finished,next_completion_time=Serv_Completions(Ac_Info,Ac_queue,prev_class,totserv,Ac_finished,latest_time,next_completion_time, thres1, thres2, lam1, lam2, f, SubPolicy, rep, Time_Sep, Left_queue)
            tot_arr_cost+=arr_cost
            tot_dep_cost+=dep_cost

        #update_elap=(time.time()-start_update_time)/conv_factor

        # else:

        # 	update_elap=0

        # print('update_elap: '+str(update_elap)+' other elap: '+str(Repop_elap+Pop_elap+elap))
        # if update_elap>Repop_elap+Pop_elap+elap:
        # 	print('=========================>>>>>>>>>>>>>>>>>>>>>>>>>>> update_elap won')

        # if int(tm*100)!=int(old_tm*100):
        # 	print(str(SubPolicy)+' tm: '+str(tm)+' elap: '+str(elap))

        old_tm=tm
        tm=latest_time

        #tm+=max(elap+Repop_elap+Pop_elap,0.00001) #max(Repop_elap+Pop_elap+elap,0.00001)
        #tm+=0.01
        #tm+=Repop_elap+Pop_elap+elap
        # if SubPolicy=='FCFS' or SubPolicy=='CLS':
        # 	tm+=0.001
        # elif len(Arr_Pool)+len(Arr_NotReady)>0:
        # 	#tm+=((time.time()-latest_timer)/conv_factor)
        # 	if mv_time==1:
        # 		mv_time=0
        # 		#print('tm: '+str(tm)+' Ov_GA_counter: '+str(Ov_GA_counter))
        # 		tm+=0.01
        # else:
        # 	tm=min(next_completion_time,tm+0.01)

    print('Final cost is '+str(tot_arr_cost+tot_dep_cost))
    gg.write(str(SubPolicy)+','+str(tot_arr_cost+tot_dep_cost)+',')
    if Policy!='Alternate':
        gg.write('Weather transitions'+','+str(wlb_tm)+','+str(wub_tm)+',')
    gg.write(str(time.time()-begin_time)+',')
    if SubPolicy=='VNS' or SubPolicy=='VNSD':
        gg.write(str(Loop_Nums)+','+str(tot_mut)+',')
    if SubPolicy=='VNS':
        gg.write(str(soln_evals_tot/soln_evals_num)+',')

    ArrTime=[0]*NoA
    ArrTime_Sorted=[0]*NoA
    ServTime=[0]*NoA
    for i in range(NoA):
        ArrTime[i]=[Ac_Info[i][9],i]
        ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        ServTime[i]=Ac_Info[i][7]

    ArrTime_Sorted.sort(key=lambda x: x[0])

    posthoc_cost=Posthoc_Check(Left_queue,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0, NoA, NormalApprox, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)
    gg.write('Posthoc Check'+','+str(posthoc_cost)+',')

    for i in range(NoA):
        g.write(str(i)+','+str(pred_serv[i])+','+str(Ac_finished[i])+'\n')

    #h.write(str(rep)+','+str(tot_arr_cost+tot_dep_cost)+'\n')

    print('Done!')

    policy_index+=1
    if policy_index==len(Policies):

        #Do Perm Heuristic

        ArrTime=[0]*NoA
        ArrTime_Sorted=[0]*NoA
        ServTime=[0]*NoA
        for i in range(NoA):
            ArrTime[i]=[Ac_Info[i][9],i]
            ArrTime_Sorted[i]=[Ac_Info[i][9],i]
            ServTime[i]=Ac_Info[i][7]

        ArrTime_Sorted.sort(key=lambda x: x[0])

        FCFS_cost=Calculate_FCFS(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,pool_max,list_min,wlb_tm,wub_tm, NoA, NormalApprox, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)
        gg.write('FCFS'+','+str(FCFS_cost)+',')

        perm_heur_cost,AC_Used=Perm_Heur(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,pool_max,list_min,wlb_tm,wub_tm, NoA, NormalApprox, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2, f1)
        gg.write('Perm Heuristic'+','+str(perm_heur_cost)+',')

        perm_heur_cost,AC_Used=Perm_Heur_New(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,pool_max,list_min,wlb_tm,wub_tm, NoA, NormalApprox, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)
        gg.write('New Perm Heuristic'+','+str(perm_heur_cost)+',')

        # #Do GA Heuristic
        # print('Starting GA Heuristic...')

        # ArrTime=[0]*NoA
        # ArrTime_Sorted=[0]*NoA
        # ServTime=[0]*NoA
        # for i in range(NoA):
        # 	ArrTime[i]=[Ac_Info[i][9],i]
        # 	ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        # 	ServTime[i]=Ac_Info[i][7]

        # ArrTime_Sorted.sort(key=lambda x: x[0])

        # GA_heur_cost,iter_no=GA_Heur(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,GA_PopSize,wlb_tm,wub_tm)

        # gg.write('GA Heuristic'+','+str(GA_heur_cost)+','+str(iter_no)+',')

        # #Do GA Heuristic (Crossover)
        # print('Starting GA Heuristic Crossover...')

        # ArrTime=[0]*NoA
        # ArrTime_Sorted=[0]*NoA
        # ServTime=[0]*NoA
        # for i in range(NoA):
        # 	ArrTime[i]=[Ac_Info[i][9],i]
        # 	ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        # 	ServTime[i]=Ac_Info[i][7]

        # ArrTime_Sorted.sort(key=lambda x: x[0])

        # GA_heur_cost,iter_no=GA_Heur_Crossover(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,GA_PopSize,wlb_tm,wub_tm)

        # gg.write('GA Heuristic Crossover'+','+str(GA_heur_cost)+','+str(iter_no)+',')

        # #Do ILS Heuristic
        # print('Starting ILS Heuristic...') #iterated local search

        # ArrTime=[0]*NoA
        # ArrTime_Sorted=[0]*NoA
        # ServTime=[0]*NoA
        # for i in range(NoA):
        # 	ArrTime[i]=[Ac_Info[i][9],i]
        # 	ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        # 	ServTime[i]=Ac_Info[i][7]

        # ArrTime_Sorted.sort(key=lambda x: x[0])

        # ILS_heur_cost,iter_no=ILS(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm)

        # gg.write('ILS Heuristic'+','+str(ILS_heur_cost)+','+str(iter_no)+',')

        # print('Starting VNS Heuristic...') #variable neighbourhood iterated local search

        # ArrTime=[0]*NoA
        # ArrTime_Sorted=[0]*NoA
        # ServTime=[0]*NoA
        # for i in range(NoA):
        # 	ArrTime[i]=[Ac_Info[i][9],i]
        # 	ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        # 	ServTime[i]=Ac_Info[i][7]

        # ArrTime_Sorted.sort(key=lambda x: x[0])

        # VNS_heur_cost,iter_no=VNS(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm)

        # gg.write('VNS Heuristic'+','+str(VNS_heur_cost)+','+str(iter_no)+',')

        # #Do SA Heuristic
        # print('Starting SA Heuristic...')

        # ArrTime=[0]*NoA
        # ArrTime_Sorted=[0]*NoA
        # ServTime=[0]*NoA
        # for i in range(NoA):
        # 	ArrTime[i]=[Ac_Info[i][9],i]
        # 	ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        # 	ServTime[i]=Ac_Info[i][7]

        # ArrTime_Sorted.sort(key=lambda x: x[0])

        # SA_heur_cost,iter_no=SA_Heur(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm)

        # gg.write('SA Heuristic'+','+str(SA_heur_cost)+','+str(iter_no)+',')

        # if 1==1: #Max_LookAhead==NoA:

        # 	#Do Tabu Heuristic
        # 	print('Starting Tabu Heuristic...')

        # 	ArrTime=[0]*NoA
        # 	ArrTime_Sorted=[0]*NoA
        # 	ServTime=[0]*NoA
        # 	for i in range(NoA):
        # 		ArrTime[i]=[Ac_Info[i][9],i]
        # 		ArrTime_Sorted[i]=[Ac_Info[i][9],i]
        # 		ServTime[i]=Ac_Info[i][7]

        # 	ArrTime_Sorted.sort(key=lambda x: x[0])

        # 	Tabu_heur_cost,iter_no=Tabu(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm)

        # 	gg.write('Tabu Heuristic'+','+str(Tabu_heur_cost)+','+str(iter_no)+',')

        gg.write('\n')
        gg.flush()

        policy_index=0

        # if Max_LookAhead==NoA:
        # 	Max_LookAhead=30
        # else:
        # 	Max_LookAhead=NoA
        # 	rep+=1

        # if rep<49:
        # 	rep+=1
        # else:
        # 	conv_factor+=0.5
        # 	rep=0

        # if wiener_sig<0.9:
        # 	wiener_sig+=0.2
        # else:
        # 	wiener_sig=0.1

        rep+=1

# print('Outputting Ac_Info array')
# print(str(Ac_Info))

f.close()
g.close()
#h.close()
gg.close()

f1.close()
#f2.close()

#rr.close()

#f5.close()

# if stepthrough==1:
#     st.close()
# if step_summ==1:
#     st2.close()
# if step_new==1:
#     st3.close()

st4.close()
