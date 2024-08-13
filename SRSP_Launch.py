# This is set up to run numerical experiments with a preset value of wiener_sig (also used as the weather variance parameter) and randomly-
# -generated values of k, thres1, pax weights, start & end times for bad weather (& need to also include random lam1 & lam2)

#from __future__ import print_function, division
import random
import hydra
from omegaconf import DictConfig, OmegaConf

from stoch_runway_scheduler import Simulation, BrownianWeatherProcess, BrownianTrajectory, read_flight_data, sample_pretac_delay, SimHeur, Genetic_determ, Perm_Heur, Perm_Heur_New, Calculate_FCFS, Posthoc_Check, Cost, ErlangSeparation

#################
# CONIFIGUATION #
#################

NO_REPS = 10000 # total number of random scenarios that we will simulate; in each scenario we evaluate the performances of different algorithms such as SimHeur, DetHeur, FCFS

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

pool_max = 6 # Used as a parameter for "perm heuristic" which searches for the best landing sequence under perfect information, i.e. assumes all random information already known
list_min = 6 # Also used only for the "perm heuristic""

# JF Note: Important - should these two things be linked?
conv_factor = 1 # no. of seconds in a minute for conversion purposes
resolution = 0.1


@hydra.main(version_base=None, config_path="conf", config_name="config")
def main(cfg: DictConfig):
    print(OmegaConf.to_yaml(cfg))

    rep = 0 #counter of which scenario we're currently on
    policy_index = 0 # indicates which policy we're currently evaluating, e.g. SimHeur, DetHeur etc (if this is zero then we take the first policy from the list of policies to be evaluated)


    flight_data, Alpha_Ps, Beta_Ps, late_means = read_flight_data(cfg.problem.flight_data,
                                                                  cfg.problem.min_ps_time,
                                                                  cfg.problem.max_ps_time, 
                                                                  cfg.trajectory.wiener_sig)

    sep_cfg = cfg.separation
    sep = ErlangSeparation(sep_cfg.time_sep, sep_cfg.k, sep_cfg.w_rho)

    #for rep in range(no_reps):
    # JF Question Why is for loop not used?: Rob not sure - may be fine to change back to for loop
    while rep < NO_REPS:

        repn = rep # int(rep/100+1)
        random.seed(repn*100) #set the random seed for generating random parameter values; seed in set according to the replication (scenario) number

        print('*** Importing the flight data...')
        #---------------------------------------------------------#
        # Read flight data file and initialise pretactical delays #
        #---------------------------------------------------------#
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

        cost_fn = Cost(cfg.cost.threshold_1, cfg.cost.threshold_2, 
                       cfg.cost.lambda_1, cfg.cost.lambda_2)

        print('*** Generating the ETA trajectory array...')
        # Trajectories are generated for whole 8 hour period for each flight

        trajecs = []
        for i in range(NoA):
            trajec = BrownianTrajectory(flight_data[i].dep_sched, ps_time[i], cfg.problem.tau, 
                                        cfg.trajectory.wiener_sig, 
                                        freq=cfg.trajectory.freq)
            trajecs.append(trajec)

        print('*** Generating the weather transition array...')

        # Long time period over which to generate weather predictions
        # We make longer in case bad weather finish time falls outside of time horizon
        T = (cfg.problem.max_ps_time - cfg.problem.min_ps_time) * 2 # The factor 2 here is probably a bit over cautious
        weather_process = BrownianWeatherProcess(cfg.weather.wlb, cfg.weather.wub, T,
                                                 cfg.weather.weather_sig, freq=cfg.weather.freq)

        release_policy = SimHeur(trajecs, sep, weather_process, cost_fn, cfg.sim_heur.s, cfg.sim_heur.l,
                                cfg.sim_heur.n_rel, cfg.sim_heur.r, cfg.sim_heur.n_repop,
                                cfg.sim_heur.s_min, cfg.sim_heur.m_mut)
        simulation = Simulation(flight_data, ps_time, pax_weight, trajecs, sep, weather_process, cfg.problem.tau,
                                release_policy, conv_factor, resolution)
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
        # # posthoc_cost = Posthoc_Check(Left_queue, Ac_Info, ArrTime, ServTime, ArrTime_Sorted, weather_process, 0, NoA, cfg.separation.w_rho, k, Time_Sep, cost_fn)

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
            # FCFS_cost = Calculate_FCFS(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process, NoA, cfg.separation.w_rho, k, Time_Sep, cost_fn)
            # perm_heur_cost, AC_Used = Perm_Heur(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process, NoA, cfg.separation.w_rho, k, Time_Sep, cost_fn, f1)
            # perm_heur_cost, AC_Used = Perm_Heur_New(Ac_Info, ArrTime, ServTime, ArrTime_Sorted, pool_max, list_min, weather_process, NoA, cfg.separation.w_rho, k, Time_Sep, cost_fn)
            # policy_index=0

            rep += 1

if __name__ == '__main__':
    main()