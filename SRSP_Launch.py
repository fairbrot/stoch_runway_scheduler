# This is set up to run numerical experiments with a preset value of wiener_sig (also used as the weather variance parameter) and randomly-
# -generated values of k, thres1, pax weights, start & end times for bad weather (& need to also include random lam1 & lam2)

from __future__ import print_function, division
import logging
import random

from stoch_runway_scheduler import Simulation, BrownianWeatherProcess, BrownianTrajectory, read_flight_data, sample_pretac_delay, SimHeur, Genetic_determ, Perm_Heur, Perm_Heur_New, Calculate_FCFS, Posthoc_Check, Cost, ErlangSeparation

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
pot_thres1 = [0, 15] # potential set of values of thres1 to sample from
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


pool_max = 6 # Used as a parameter for "perm heuristic" which searches for the best landing sequence under perfect information, i.e. assumes all random information already known
list_min = 6 # Also used only for the "perm heuristic""

l = 15 # NoA # This is the length of a sequence, equivalent to parameter l in paper  - in paper this is 15
n_rel = 50
n_repop = 500
r = 50
S = 20 # Initial number of sequences in the population, written as S in paper (see Section 3.1)
S_min = 10 # Length of shortlist JF - Perhaps move to parameters
m_mut = 25 # important parameter, determines how many non-improving heuristic moves need to be made before a mutation is carried out; this is written as m_{mut} in paper (see the flow chart, Figure 3)

# JF Note: Important - should these two things be linked?
conv_factor = 1 # no. of seconds in a minute for conversion purposes
freq = 100 # number of updates for each minute
resolution = 0.1

###########
# LOGGING #
###########

##################
# INITIALISATION #
##################

# if Policy=='Alternate':
#   SubPolicy='Perm'
# else:
#   SubPolicy=Policy

def main():

    rep = 0 #counter of which scenario we're currently on
    policy_index = 0 # indicates which policy we're currently evaluating, e.g. SimHeur, DetHeur etc (if this is zero then we take the first policy from the list of policies to be evaluated)


    flight_data, Alpha_Ps, Beta_Ps, late_means = read_flight_data(DATA_DIR + '/flight_pretac_data.csv',
                                                                min_ps_time, max_ps_time, wiener_sig)

    #for rep in range(no_reps):
    # JF Question Why is for loop not used?: Rob not sure - may be fine to change back to for loop
    while rep < no_reps:

        repn = rep # int(rep/100+1)
        random.seed(repn*100) #set the random seed for generating random parameter values; seed in set according to the replication (scenario) number
        print('*** Importing the flight data...')
        #---------------------------------------------------------#
        # Read flight data file and initialise pretactical delays #
        #---------------------------------------------------------#
        # Orig_Ps = sched arrival time
        # Dep_Ps = time at which tactical delay begins (called h in paper)

        pretac_delays = [sample_pretac_delay(alpha, beta, fdata.arr_sched, fdata.dep_sched, x_bar) for (alpha, beta, fdata, x_bar) in zip(Alpha_Ps, Beta_Ps, flight_data, late_means)]
        # this stores the adjusted scheduled times for aircraft after applying the random pre-tactical delay
        ps_time = [fdata.arr_sched + pretac_d for (fdata, pretac_d) in zip(flight_data, pretac_delays)]

        NoA = len(flight_data)
        print('No. of ACs: '+str(NoA))

        # --------------------------- #
        # Random parameter generation #
        #-----------------------------#
        pax_weight = [0]*NoA # stores the randomly-generated cost weightings for aircraft, based on (hypothetical) numbers of passengers carried; written as g_i in the paper (see objective function (13))
        for i in range(NoA):
            # Passenger weight: this is called g_i in the paper
            # In objective? Shouldn't these be the same for every run? Rob says perhaps, but decided to use random ones for each replication.
            # Rob views this as similar to changing the weather
            match flight_data[i].ac_class:
                case 0:
                    pax_weight[i] = 0.2*random.random() + 0.8 # Flights in the 'heavy' class have a passenger weight between 0.8 and 1
                case 1 | 2:
                    pax_weight[i] = 0.2*random.random() + 0.6 # Flights in the 'upper medium' or 'lower medium' class have a passenger weight between 0.6 and 0.8 
                case _:
                    pax_weight[i] = 0.2*random.random() + 0.4 # Flights in the 'small' class have a passenger weight between 0.4 and 0.6

        SubPolicy = Policies[policy_index] # SubPolicy indicates the policy we are currently considering (e.g. SimHeur, DetHeur)

        # Randomly generate the Erlang service time parameter

        # Results can be stratified by k (roughlt 1 fifth of runs for each value of k)
        k = random.choice(pot_k)
        print('k: '+str(k))

        sep = ErlangSeparation(Time_Sep, k, w_rho)

        # These are random for similar reasons pax_weight (g_i in paper) - results may be stratified by this as well
        # lam1 and lam2 are the weights of scheduling delay and airborne holding delays - these are called theta^S and theta^W in the paper
        lam1 = random.choice(pot_lam1) # Random lam1
        lam2 = 1-lam1
        print('lam1: '+str(lam1)+' lam2: '+str(lam2))

        # Randomly generate thres1 (thres2 is set above)
        # Random so results could potentially be stratified
        thres1 = random.choice(pot_thres1) # 15 means allow 15 minutes schedule delay

        cost_fn = Cost(thres1, thres2, lam1, lam2)

        print('*** Generating the ETA trajectory array...')
        # Trajectories are generated for whole 8 hour period for each flight

        trajecs = []
        for i in range(NoA):
            trajec = BrownianTrajectory(flight_data[i].dep_sched, ps_time[i], tau, wiener_sig, freq=freq)
            trajecs.append(trajec)

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
        weather_process = BrownianWeatherProcess(wlb, wub, T, weather_sig, freq=freq)

        release_policy = SimHeur(trajecs, sep, weather_process, cost_fn, S, l, 
                                n_rel, r, n_repop, S_min, m_mut)
        simulation = Simulation(flight_data, ps_time, pax_weight, trajecs, sep, weather_process, tau, release_policy,
                                conv_factor, resolution)
        print(f'*** Into main loop for rep {rep} and policy {SubPolicy}...')
        simulation.run()


        # print('Final cost is '+str(tot_cost))


        # ArrTime = [0]*NoA
        # ArrTime_Sorted = [0]*NoA
        # ServTime = [0]*NoA
        # for i in range(NoA):
        #     ArrTime[i] = [Ac_Info[i].pool_time,i]
        #     ArrTime_Sorted[i] = [Ac_Info[i].pool_time,i]
        #     ServTime[i] = Ac_Info[i].service_rns

        # ArrTime_Sorted.sort(key=lambda x: x[0])
        # # Posthoc check validates that flight statistics are consistent with costs
        # # posthoc_cost = Posthoc_Check(Left_queue, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process, 0, NoA, w_rho, k, Time_Sep, cost_fn)

        # print('Done!')

        # JF Question: What is happening here?
        policy_index += 1
        if policy_index == len(Policies):
            # Do Perm Heuristic
            # ArrTime = [0]*NoA
            # ArrTime_Sorted = [0]*NoA
            # ServTime = [0]*NoA
            # for i in range(NoA):
            #     ArrTime[i] = [Ac_Info[i].pool_time,i]
            #     ArrTime_Sorted[i] = [Ac_Info[i].pool_time,i]
            #     ServTime[i] = Ac_Info[i].service_rns

            # ArrTime_Sorted.sort(key=lambda x: x[0])
            # FCFS_cost = Calculate_FCFS(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process, NoA, w_rho, k, Time_Sep, cost_fn)
            # perm_heur_cost, AC_Used = Perm_Heur(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process, NoA, w_rho, k, Time_Sep, cost_fn, f1)
            # perm_heur_cost, AC_Used = Perm_Heur_New(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process, NoA, w_rho, k, Time_Sep, cost_fn)
            # policy_index=0

            rep += 1

if __name__ == '__main__':
    main()