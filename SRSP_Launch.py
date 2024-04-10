# This is set up to run numerical experiments with a preset value of wiener_sig (also used as the weather variance parameter) and randomly-
# -generated values of k, thres1, pax weights, start & end times for bad weather (& need to also include random lam1 & lam2)

from __future__ import print_function, division
import math
import logging
import random
import time
import os

import numpy as np 

from stoch_runway_scheduler import generate_weather, generate_trajectory, read_flight_data, sample_pretac_delay, weather, Genetic, Genetic_determ, Populate, Repopulate_VNS, sample_cond_gamma, getcost, Annealing_Cost, Perm_Heur, Perm_Heur_New, Calculate_FCFS, sample_gamma, gamma_create_cdf, Posthoc_Check, Update_Stats, Update_ETAs, Serv_Completions, FlightInfo, FlightStatus, SequenceInfo

#################
# CONIFIGUATION #
#################

no_reps = 10000 # total number of random scenarios that we will simulate; in each scenario we evaluate the performances of different algorithms such as SimHeur, DetHeur, FCFS

#-----------------#
# Problem Options #
#-----------------#

DATA_DIR = '/home/jamie/Insync/fairbrot@lancaster.ac.uk/OneDrive Biz - Shared/SRSP data files'

# JF: should remove this when we can as NoA is set after reading data file
NoA = 700 # number of aircraft - temporary value which will get changed later
S = 40 # number of time slots (Rob not sure whether this is needed)

w_rho = 10/9 # separation multiplier for bad weather (multiplies mean, not rate, so should be >1); we have used a reduction of 10% due to bad weather based on Odoni et al (2011) cited in Shone et al (2021)
wiener_sig = 0.1 # standard deviation for Brownian motion
weather_sig = wiener_sig # this assumption is being made in the paper for simplicity
print('wiener_sig: '+str(wiener_sig))

tau = 30 # Determines when an aircraft is deemed to have entered the pool. E.g. if tau=30, aircraft enters pool when its ETA is 30 minutes away from the current time.
t = 15 # length of a time slot in minutes

NoC = 4 # no. of aircraft classes
# Time separations in seconds taken from Bennell et al (2017) with H, U, M, S as the 4 classes; 
# the 5th array is for the situation where there is no leading aircraft
Time_Sep = [[97,121,121,145],[72,72,72,97,97],[72,72,72,72],[72,72,72,72],[72,72,72,72]] 
# JF: Time_Sep is List[List[int]]

# k controls variances of service times - larger means less variance
# Chosen values correspond to specified values of coefficients of variation
pot_k = [16, 25, 44, 100, 400] # Potential set of k values to sample from

# lam1 and lam2 are the weights of scheduling delay and airborne holding delays - these are called theta^S and theta^W in the paper
pot_lam1 = [0.1, 0.3, 0.5, 0.7, 0.9] # Potential set of values for lambda 1 to sample from

# These correspond to gamma^S (thres1) and gamma^W (thres2) in paper
pot_thres1 = [0,15] # potential set of values of thres1 to sample from
thres2 = 0

# Min and max prescheduled time to consider in solution approach
# Measured in minutes from from midnight? (6AM-2PM is when simulation runs between)
min_ps_time = 360 # inclusive
max_ps_time = 840 # non-inclusive

#-------------------#
# Algorithm Options #
#-------------------#

# JF: policy is an scheduling policy algorithm - may not affect anything now
# Alternate means sway between VNS and VNSD
Policy = 'Alternate' #FCFS, Perm, SA, GA or Alternate

Use_VNS = 1
Use_VNSD = 1
# Use_FCFS=1

Policies = []
if Use_VNS == 1:
    Policies.append('VNS')
if Use_VNSD == 1:
    Policies.append('VNSD')
# if Use_FCFS==1:
#   Policies.append('FCFS')

# JF Question: what is this? Would be good to avoid setting it to NoA before data has been read
max_FCFS = NoA # int(NoA/2)

Max_LookAhead = 15 # NoA # This is the length of a sequence, equivalent to parameter l in paper  - in paper this is 15

pool_max = 6 # Used as a parameter for "perm heuristic" which searches for the best landing sequence under perfect information, i.e. assumes all random information already known
list_min = 6 # Also used only for the "perm heuristic""

GA_PopSize = 20 # Initial number of sequences in the population, written as S in paper (see Section 3.1)
VNS_limit = 25 # important parameter, determines how many non-improving heuristic moves need to be made before a mutation is carried out; this is written as m_{mut} in paper (see the flow chart, Figure 3)

# JF Note: Important - should these two things be linked?
conv_factor = 1 # no. of seconds in a minute for conversion purposes
freq = 100 # number of updates for each minute

###########
# LOGGING #
###########

# JF: these three options seem to be for logging purposes
#     they are broken for now as a separate output stream is created for each one, and all these currently
#     need to be passed a function where they are used.
stepthrough = 0
step_summ = 0
step_new = 0

# Logs can similarly be set up for step_summ and step_new
stepthrough_logger = logging.getLogger('stepthrough')
stepthrough_logger.setLevel(level=logging.ERROR)
c_handler = logging.FileHandler("SRSP_stepthrough.log", mode='w')
# c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
# c_handler.setFormatter(c_format)
stepthrough_logger.addHandler(c_handler)

step_summ_logger = logging.getLogger('step_summ')
step_new_logger = logging.getLogger('step_new')

f = open("SRSP_out_%s_%s.csv" % (Policy, conv_factor), "w") # detailed results - for each rep what happened for each aircraft
g = open("AC_predictions_%s_%s.csv" % (Policy, conv_factor), "w")
#h = open("cost_log_%s_%s.csv" % (Policy, conv_factor), "w")
gg = open("SRSP_rep_results_%s_%s.csv" % (Policy, conv_factor), "w") # summaries for each rep

f1 = open("Perm_Heur_out_%s_%s.csv" % (Policy, conv_factor), "w") # output from a particular heuristic - may not be needed (can remove if you like)

f.write('Policy'+',''Rep'+','+'AC'+','+'Flight Num'+','+'Prev Class'+','+'Cur Class'+','+'Time Sep'+','+'Orig PS time'+','+'PS time'+','+'Pool Arrival'+','+'Release Time'+','+'Travel Time'+','+'Weather Coeff'+','+'Enters Serv'+','+'Actual Serv'+','+'Ends Serv'+','+'Lateness'+','+'Queue Delay'+','+'Pax Weight'+','+'Cost'+','+'counter'+','+'qp'+','+'Predicted total'+'\n')

##################
# INITIALISATION #
##################

Ov_GA_counter = 0 # only used for stepthrough purposes; to do with counting how many times the 'Genetic' function has been called
VNS_counter = 0 # this counts how many non-improving heuristic moves we've made since the last reset
tot_mut = 0 # counts how many total mutations we've made; really just for output purposes

# Presumably indicates for each aircraft whether it has landed yet
# 'finished' means aircraft has completed landing, i.e. service time has finished
Ac_finished = [0]*NoA
Ac_Info = [0]*NoA
pax_weight = [0]*NoA # stores the randomly-generated cost weightings for aircraft, based on (hypothetical) numbers of passengers carried; written as g_i in the paper (see objective function (13))

# if Policy=='Alternate':
#   SubPolicy='Perm'
# else:
#   SubPolicy=Policy

rep = 0 #counter of which scenario we're currently on
policy_index = 0 # indicates which policy we're currently evaluating, e.g. SimHeur, DetHeur etc (if this is zero then we take the first policy from the list of policies to be evaluated)

#for rep in range(no_reps):
# JF Question Why is for loop not used?: Rob not sure - may be fine to change back to for loop
while rep < no_reps:

    repn = rep # int(rep/100+1)
    random.seed(repn*100) #set the random seed for generating random parameter values; seed in set according to the replication (scenario) number
    print('*** Importing the flight data...')
    #---------------------------------------------------------#
    # Read flight data file and initialise pretactical delays #
    #---------------------------------------------------------#
    flight_id, Ac_class, Orig_Ps, Dep_Ps, Alpha_Ps, Beta_Ps, late_means = read_flight_data(DATA_DIR + '/flight_pretac_data.csv',
                                                                                            min_ps_time, max_ps_time, wiener_sig)
    pretac_delays = [sample_pretac_delay(a, b, ps_t, hi, l_mn) for (a, b, ps_t, hi, l_mn) in zip(Alpha_Ps, Beta_Ps, Orig_Ps, Dep_Ps, late_means)]
    # this stores the adjusted scheduled times for aircraft after applying the random pre-tactical delay
    Arr_Ps = [orig_ps + pretac_d for (orig_ps, pretac_d) in zip(Orig_Ps, pretac_delays)]

    NoA = len(flight_id)
    print('No. of ACs: '+str(NoA))

    # --------------------------- #
    # Random parameter generation #
    #-----------------------------#
    for i in range(NoA):
        # Passenger weight: this is called g_i in the paper
        # In objective? Shouldn't these be the same for every run? Rob says perhaps, but decided to use random ones for each replication.
        # Rob views this as similar to changing the weather
        if Ac_class[i] == 0:
            pax_weight[i] = 0.2*random.random() + 0.8 # Flights in the 'heavy' class have a passenger weight between 0.8 and 1
        elif Ac_class[i] == 1 or Ac_class[i] == 2:
            pax_weight[i] = 0.2*random.random() + 0.6 # Flights in the 'upper medium' or 'lower medium' class have a passenger weight between 0.6 and 0.8 
        else:
            pax_weight[i] = 0.2*random.random() + 0.4 # Flights in the 'small' class have a passenger weight between 0.4 and 0.6

    SubPolicy = Policies[policy_index] # SubPolicy indicates the policy we are currently considering (e.g. SimHeur, DetHeur)

    if SubPolicy == 'GA' or SubPolicy == 'VNS':
        GA_LoopSize = 500 # This is
    elif SubPolicy == 'GAD' or SubPolicy == 'VNSD':
        GA_LoopSize = 1

    # Randomly generate the Erlang service time parameter

    # Results can be stratified by k (roughlt 1 fifth of runs for each value of k)
    k = random.choice(pot_k)
    print('k: '+str(k))

    # These are random for similar reasons pax_weight (g_i in paper) - results may be stratified by this as well
    # lam1 and lam2 are the weights of scheduling delay and airborne holding delays - these are called theta^S and theta^W in the paper
    lam1 = random.choice(pot_lam1)     # Random lam1
    lam2 = 1-lam1
    print('lam1: '+str(lam1)+' lam2: '+str(lam2))

    # Randomly generate thres1 (thres2 is set above)
    # Random so results could potentially be stratified
    thres1 = random.choice(pot_thres1) # 15 means allow 15 minutes schedule delay

    print('*** Creating the Gamma CDF...')
    gamma_cdf = gamma_create_cdf(k)

    print('*** Generating initial aircraft info...')
    # this will be a multi-dimensional list storing lots of information about each aircraft; see later
    Ac_Info = [0]*NoA

    for i in range(NoA):
        # Generating service times for arrivals - these are scheduled to have the right mean later
        # JF Question: What is ServPercs? Below RS says this is RNs used for service time
        ServPercs=np.random.gamma(k,1)

        Ac_Info[i] = FlightInfo(FlightStatus.NOT_READY, Ac_class[i], Arr_Ps[i], Arr_Ps[i],
                      0, 0, 0, ServPercs,
                      0, 0, pax_weight[i], False,
                      1, 0, 0, 0,
                      0, Dep_Ps[i], Orig_Ps[i], flight_id[i])

    Ac_Info.sort(key=lambda x: x.ps_time)

    print('*** Generating the ETA trajectory array...')
    # Trajectories are generated for whole 8 hour period for each flight

    Brown_Motion = []
    for i in range(NoA):
        Ps_time, Dep_time,  = Ac_Info[i].ps_time, Ac_Info[i].sched_dep_time
        pool_arr_time, travel_time, brown_motion = generate_trajectory(Dep_time, Ps_time, tau, wiener_sig, freq=freq)
        Ac_Info[i].travel_time = travel_time
        Ac_Info[i].pool_time = pool_arr_time
        Brown_Motion.append(brown_motion)

    header = ['AC', 'Class', 'PS time', 'Pool arrival', 'Travel time', 'Runway time']
    stepthrough_logger.info(', '.join(header))
    for AC in range(NoA):
        Ac_Infoi = Ac_Info[AC]
        stepthrough_logger.info('%s, %s, %s, %s, %s, %s, %s', AC, Ac_Infoi.ac_class, Ac_Infoi.ps_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time)
    stepthrough_logger.info('\n')

    print('*** Generating the weather transition array...')

    # 4 possible cases: no bad weather, 30 minutes of bad weather, 60 minutes of bad weather or 120 minutes of bad weather 
    # (bad weather is always forecast for the middle of the day)
    # E.g. 300 is 10AM according to the rescaling of time used earlier
    # *Predicted* start and end times of bad weather
    # Random for result stratification wlb is T_0(0) and wub is T_1(0) in paper
    wlb, wub = random.choice([(0,0), (285, 315), (270, 330), (240, 360)])

    # Long time period over which to generate weather predictions
    # We make longer in case bad weather finish time falls outside of time horizon
    T = (max_ps_time - min_ps_time) * 2 # The factor 2 here is probably a bit over cautious
    wlb_tm, wub_tm, weather_lb, weather_ub = generate_weather(wlb, wub, T, weather_sig, freq=freq)

    stepthrough_logger.info('wlb_tm: %d wub_tm: %d', wlb_tm, wub_tm)


    # VNS = SimHeur, VNSD = DetHeur
    if SubPolicy in ('VNS','VNSD'): 
        # Generate the initial population of sequences
        print('Generating initial population of sequences...')

        # Should probably just be Max_LookAhead as values haven't changed at this point
        no_ACs = min(Max_LookAhead, NoA) # JF: perhaps change this name to be clearer (it is sequence length)
        # JF Question: should this be FCFS? No - this is first scheduled first served
        FSFS_seq = [i for i in range(no_ACs)]

        # Arg 3 is Arr_Pool which is initially empty at this point
        GA_PopList, GA_Info = Populate(Ac_Info, FSFS_seq, [], FSFS_seq, GA_PopSize, Max_LookAhead)

        Opt_Seq = FSFS_seq[:]
        OptCost = 1000000 # initial cost
        queue_probs = [0]*NoA # JF Question: this is shared between solutions? Should this really be the case? Possibly is independent of sequence
        Opt_List = [] # For storing SequenceInfo of best solutions
        Opt_Seqs = []
        Opt_Size = 10 # Length of shortlist JF - Perhaps move to parameters

        # JF Note: we may now not need separate lists Opt_Seqs and GA_PopList - Rob thinks a list of best sequences probably isn't needed now
        # Only really need to keep best sequence and mutate from that.
        # When Rob wrote this code originally he had idea of using best n sequences for generating new sequences
        # Create an initial population of sequences for Opt_Seq
        # which shouldn't intersect GA_PopList
        while len(Opt_List) < Opt_Size:
            new_seq = random.sample(FSFS_seq, k=len(FSFS_seq)) # shuffle fcfs sequence
            if new_seq not in GA_PopList and new_seq not in Opt_Seqs: # This is just to initialise opt_list - we want to compare with distinct GA_PopList
                Opt_List.append(SequenceInfo(new_seq[:],0,0,queue_probs,0)) # analogous to GA_Info
                Opt_Seqs.append(new_seq[:])

        GA_Check_Increment = GA_LoopSize / 10 # Called r in paper - how often to do ranking and selection
        # Is done for every iteration for VNSD

        GA_counter = 0 # called n in paper (number of trajectories sampled so far)
        GA_CheckSize = GA_Check_Increment
        print(f'GA_CheckSize: {GA_CheckSize}')
        Ov_GA_counter = 0 # Rob thinks for logging purposes

        stepthrough_logger.info('Initial population of sequences:')
        for (j, poplist) in enumerate(GA_PopList):
            stepthrough_logger.info("%d, %s\n", j, poplist)

        print(f'Ov_GA_counter: {Ov_GA_counter}')

    for AC in range(NoA):
        print(f'AC: {AC} Class: {Ac_Info[AC].ac_class} Orig Ps time: {Ac_Info[AC].orig_sched_time} Ps time: {Ac_Info[AC].ps_time} Arrives in pool: {Ac_Info[AC].pool_time} Travel time: {Ac_Info[AC].travel_time}')

    Ac_queue = [] # flights currently in queue for landing
    Left_queue = [] # Aircraft which have landed
    Arr_Pool = [] # Aircraft in Pool
    Arr_NotReady = [] # Aircraft that haven't joined pool yet
    Dep_NotReady = [] # Is this needed? It is used once below - Possibly not - check how it is used

    # JF Question: I need a lot of these parameters explaining

    tm = 0  # tm is current time (in minutes)
    old_tm = tm # previous time
    totserv = 0 # counter of number of aircraft served (landed)
    latest_class = 4 # class of latest aircraft to be added to the queue, initially set to 4 (dummy value)
    prev_class = 4 # class of previous aircraft to be served, initially set to 4 (dummy value)
    Ac_added = [] # aircraft to be added to the queue
    counter = 0
    qp = 0 # Possibly related to an old permutation based heuristic - may be redundant now
    weather_state = 0 # 0 means it's good weather, 1 means bad weather, 2 means good weather again
    real_queue_complete = 0 # stores time last place was serviced
    next_completion_time = 0
    Pop_elap = 0 # Used to track time spent in certain functions - used to update simulation clock
    max_d = 1 # Possibly redundant
    pruned = 0 # indicator of whether or not number of sequences has gone below minimum
    soln_evals_tot = 0
    soln_evals_num = 0
    tot_mut = 0 # number of times we have done mutation step

    tot_arr_cost = 0
    tot_dep_cost = 0

    Loop_Nums = 0
    Loop_Evals = 0



    # Generate the initial pool and notready arrays
    print('*** Generating the initial pool...')
    for i in range(NoA):
        Ac_Infoi = Ac_Info[i]
        if Ac_Infoi.eta - tau <= 0:
            Arr_Pool.append(i)
            print('Aircraft '+str(i)+' initially included in pool (ETA is '+str(Ac_Infoi.eta)+')')
            Ac_Infoi.status = FlightStatus.IN_POOL
        else:
            Arr_NotReady.append(i)

    print('Arr_Pool:', Arr_Pool)
    print('Ac_queue:', Ac_queue)
    print('Arr_NotReady:', Arr_NotReady)

    print(f'*** Into main loop for rep {rep} and policy {SubPolicy}...')
    begin_time = time.time()
    mv_time = 0 # some sort of indicator
    if policy_index == 0:
        gg.write(str(rep)+','+str(conv_factor)+','+str(wiener_sig)+','+str(k)+','+str(lam1)+','+str(lam2)+','+str(wlb)+','+str(wub)+','+str(wlb_tm)+','+str(wub_tm)+','+str(thres1)+','+str(thres2)+',')

    initial_time = time.time()

    while totserv < NoA:

        stepthrough_logger.info('tm is %d', tm)
        stepthrough_logger.info('Arr_NotReady is %s', Arr_NotReady)
        stepthrough_logger.info('Arr_Pool is %s', Arr_Pool)
        stepthrough_logger.info('Ac_queue is %s', Ac_queue)
        stepthrough_logger.info('Left_queue is %s', Left_queue)
        stepthrough_logger.info('Ac_added is %s', Ac_added)

        # If one aircraft needs releasing and the first aircraft to be released is in the pool
        # tm >= may be redundant
        if tm >= 0 and len(Ac_added) > 0 and Ac_added[0] in Arr_Pool:
            if SubPolicy in ('VNS','VNSD'):
                Opt_List.sort(key=lambda x: x.v) # sort by mean V_s^n
                GA_Info.sort(key=lambda x: x.v)
                # Set up base sequence and pred_cost (latest predicted cost)
                if len(Opt_List) > 0 and len(GA_Info) > 0: # If both lists non-empty then compare best in one with best inother
                    if Opt_List[0][2] < GA_Info[0][2]:
                        base_seq = Opt_List[0].sequence[:]
                        pred_cost = Opt_List[0].v
                    else:
                        base_seq = GA_Info[0].sequence[:]
                        pred_cost = GA_Info[0].v
                elif len(Opt_List) > 0:
                    base_seq = Opt_List[0].sequence[:]
                    pred_cost = Opt_List[0].v
                elif len(GA_Info) > 0:
                    base_seq = GA_Info[0].sequence[:]
                    pred_cost = GA_Info[0].v
                else:
                    assert 1 == 2
            # For each aircraft being released...
            for AC in Ac_added:
                if AC in Arr_Pool:
                    Arr_Pool.remove(AC)
                    Ac_queue.append(AC)
                    if SubPolicy in ('VNS','VNSD'):
                        #print('Added AC '+str(AC)+' to the queue, counter is '+str(Ov_GA_counter)+', qp is '+str(qp))
                        Ac_Info[AC].pred_cost = pred_cost
                        stepthrough_logger.info('Added AC %d to the queue, counter is %d, qp is %.2f\n', AC, Ov_GA_counter, qp)
                        step_summ_logger.info('Added AC %d to the queue, counter is %d, qp is %s', AC, Ov_GA_counter, qp)
                        step_new_logger.info('Added AC %d to the queue, counter is %d, qp is %s', AC, Ov_GA_counter, qp)
                        base_seq.remove(AC)

                    else:
                        print('Added AC '+str(AC)+' to the queue')

                    # Gets important statistics about aircraft being released and serviced
                    # some of these outputs are used for simulation itself
                    real_queue_complete, next_completion_time, latest_class, Ov_GA_counter = Update_Stats(tm,AC,Ac_Info,Ac_queue,real_queue_complete,wlb_tm,wub_tm,latest_class,Ov_GA_counter,next_completion_time, k, Time_Sep, w_rho, SubPolicy, counter, qp)

                else:
                    # Important: if aircraft not in the pool then don't consider any others (order in Ac_Added comes from a sequence and is important)
                    break

            # We have released flights so we need to reinitialise population (Step 4A in paper)
            if SubPolicy in ('VNS','VNSD'): # This is always the case FCFS policy tested elsewhere so this is redundant - we are always using one of these policies
                GA_PopList, GA_Info = Populate(Ac_Info, base_seq, Arr_Pool, Arr_NotReady, GA_PopSize, Max_LookAhead)
                queue_probs = [0]*(len(Arr_Pool)+len(Arr_NotReady))

                Opt_List = []
                Opt_Seqs = []
                c = 0 # JF Question: Is this actually used?
                # This looks very similar to initialisation of Opt_List above
                while len(Opt_List)<Opt_Size and c<25:
                    new_seq = base_seq[:]
                    random.shuffle(new_seq)
                    if new_seq not in GA_PopList and new_seq not in Opt_Seqs:
                        Opt_List.append(SequenceInfo(new_seq[:],0,0,queue_probs,0))
                        Opt_Seqs.append(new_seq[:])
                        c = 0
                    else:
                        c += 1
                GA_counter = 0 # reset counter
                GA_CheckSize = GA_Check_Increment # next point to check
            else:
                Pop_elap = 0 # redundant case now since always using VNS or VNSD


        # JF Question: what does this part do?
        if SubPolicy in ('VNS','VNSD'):
            Repop_elap = 0
            if len(Arr_Pool) + len(Arr_NotReady) > 0: #4: # JF Question - Should this condition be automatically satisfied? See condition on outermost loop
                # Condition for step 4C in paper
                if GA_counter >= GA_LoopSize or pruned == 1: #or max_d<0.01:
                    Loop_Nums+=1
                    Loop_Evals += GA_counter
                    # Create more seqeuences from best current sequence
                    GA_PopList, GA_Info, Opt_Seq,OptCost, Opt_List, VNS_counter, tot_mut = Repopulate_VNS(GA_PopList, GA_Info, Arr_Pool, Arr_NotReady, GA_PopSize, Opt_Seq, 
                                                                                                          OptCost, Opt_List,Opt_Size,Max_LookAhead,VNS_counter,VNS_limit,tot_mut, 
                                                                                                          stepthrough, step_summ, step_new)
                    GA_counter = 0
                    GA_CheckSize = GA_Check_Increment
                    mv_time = 1
            else:
                mv_time = 1 # May be redundant
        else:
            Repop_elap = 0

        # First, get the AC List

        # JF Question: what is the condition?
        if tm >= 0 and tm <= wub_tm and int(tm*freq) != int(old_tm*freq):
            # Permute and update the weather
            wlb = weather_lb[int(tm*freq)] # random.gauss(wlb, 0.05)
            wub = weather_ub[int(tm*freq)]

        if len(Arr_Pool) + len(Arr_NotReady) > 0:
            if SubPolicy == 'VNS':
                # JF Question: what is happening here?
                Ac_added, counter, qp, max_d, pruned, GA_CheckSize, GA_counter, soln_evals_tot, soln_evals_num = Genetic(Ac_Info, Arr_Pool, Arr_NotReady, Ac_queue, Left_queue, max(tm,0), NoA, k, prev_class, GA_PopList, GA_Info, GA_LoopSize, GA_CheckSize, GA_counter, tot_arr_cost + tot_dep_cost, wlb, wub, Opt_List, max_d, soln_evals_tot, soln_evals_num, gamma_cdf, tau, Max_LookAhead, Time_Sep, thres1, thres2, lam1, lam2, GA_Check_Increment, Opt_Size, w_rho, stepthrough, wiener_sig, weather_sig)
                Ov_GA_counter+=1
                stepthrough_logger.info('GA_counter is %d', GA_counter)
            elif SubPolicy=='VNSD':
                Ac_added, counter, qp, stored_queue_complete = Genetic_determ(Ac_Info, Arr_Pool, Arr_NotReady, Ac_queue, Left_queue, max(tm,0), NoA, k, prev_class, GA_PopList, GA_Info, wlb, wub, Opt_List, tau, Max_LookAhead, Time_Sep, thres1, thres2, lam1, lam2, tot_arr_cost, tot_dep_cost, w_rho, stepthrough, step_summ, step_new)
                Ov_GA_counter += 1
                GA_counter += 1
                stepthrough_logger.info('GA_counter is %d', GA_counter)

        else:
            Ac_added, elap, counter, qp = [], 0.1, 0, 0 # JF Question: is this a syntax error?

        latest_time = (time.time() - initial_time)/conv_factor

        # JF Question: what is this condition? freq replaced 100 here
        if int(tm*freq) != int(old_tm*freq):
            Update_ETAs(Ac_Info, Arr_NotReady, Dep_NotReady, Ac_queue, tm, Brown_Motion, Arr_Pool, tau, freq)

        if len(Ac_queue) > 0 and tm >= next_completion_time: #len(Ac_queue)>0:
            arr_cost, dep_cost, totserv, prev_class, Ac_finished, next_completion_time = Serv_Completions(Ac_Info, Ac_queue, prev_class, totserv, Ac_finished, latest_time, next_completion_time, thres1, thres2, lam1, lam2, f, SubPolicy, rep, Time_Sep, Left_queue)
            tot_arr_cost += arr_cost
            tot_dep_cost += dep_cost

        old_tm = tm
        tm = latest_time

    print('Final cost is '+str(tot_arr_cost+tot_dep_cost))
    gg.write(str(SubPolicy)+','+str(tot_arr_cost+tot_dep_cost)+',')
    if Policy != 'Alternate':
        gg.write('Weather transitions'+','+str(wlb_tm)+','+str(wub_tm)+',')
    gg.write(str(time.time()-begin_time)+',')
    if SubPolicy == 'VNS' or SubPolicy == 'VNSD':
        gg.write(str(Loop_Nums)+','+str(tot_mut)+',')
    if SubPolicy == 'VNS':
        gg.write(str(soln_evals_tot/soln_evals_num)+',')

    ArrTime = [0]*NoA
    ArrTime_Sorted = [0]*NoA
    ServTime = [0]*NoA
    for i in range(NoA):
        ArrTime[i] = [Ac_Info[i].pool_time,i]
        ArrTime_Sorted[i] = [Ac_Info[i].pool_time,i]
        ServTime[i] = Ac_Info[i].service_rns

    ArrTime_Sorted.sort(key=lambda x: x[0])

    posthoc_cost = Posthoc_Check(Left_queue,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0, NoA, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)
    gg.write('Posthoc Check'+','+str(posthoc_cost)+',')

    for i in range(NoA):
        g.write(str(i)+','+','+str(Ac_finished[i])+'\n')

    print('Done!')

    policy_index += 1
    if policy_index == len(Policies):
        # Do Perm Heuristic
        ArrTime = [0]*NoA
        ArrTime_Sorted = [0]*NoA
        ServTime = [0]*NoA
        for i in range(NoA):
            ArrTime[i] = [Ac_Info[i].pool_time,i]
            ArrTime_Sorted[i] = [Ac_Info[i].pool_time,i]
            ServTime[i] = Ac_Info[i].service_rns

        ArrTime_Sorted.sort(key=lambda x: x[0])

        FCFS_cost = Calculate_FCFS(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, wlb_tm, wub_tm, NoA, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)
        gg.write('FCFS'+','+str(FCFS_cost)+',')

        perm_heur_cost, AC_Used = Perm_Heur(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, wlb_tm, wub_tm, NoA, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2, f1)
        gg.write('Perm Heuristic'+','+str(perm_heur_cost)+',')

        perm_heur_cost, AC_Used = Perm_Heur_New(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, wlb_tm, wub_tm, NoA, w_rho, k, Time_Sep, thres1, thres2, lam1, lam2)
        gg.write('New Perm Heuristic'+','+str(perm_heur_cost)+',')

        gg.write('\n')
        gg.flush()

        policy_index=0

        rep += 1


f.close()
g.close()
gg.close()
f1.close()
