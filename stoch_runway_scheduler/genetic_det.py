from typing import List
import logging
from .utils import FlightInfo, Cost, FlightStatus
from .weather import WeatherProcess, WeatherStatus, StochasticWeatherProcess
from .sequence import SequenceInfo
from .trajectory import StochasticTrajectory
from .populate import Populate, Repopulate_VNS
from .state import State
from .separation import StochasticSeparation, landing_time
from .gamma import gamma_cond_exp

log = logging.getLogger(__name__)

class DetHeur:

    def __init__(self, sep: StochasticSeparation, weather: StochasticWeatherProcess, tau: float, cost_fn: Cost,
                S: int, l: int, S_min: int, m_mut: int):
        # For simulation
        self.sep = sep
        self.weather = weather
        self.tau = tau

        # Algorithm params
        self.cost_fn = cost_fn
        self.S = S                
        self.S_min = S_min
        self.l = l
        self.m_mut = m_mut

        # Internal variables
        self.GA_counter = 0
        self.VNS_counter = 0
        self.tot_mut = 0
        self.reset_pop = True
        self.GA_Info = []

        no_ACs = l # JF Note: this will fail if l is than initial number of aircraft
        self.base_seq = [i for i in range(no_ACs)] # Initial value assumes AcInfo is ordered by Ps_time
        "Base sequence from which a new population is generated"

    def run(self, state: State) -> list[int]:
        # Can't run DetHeur if no flights remain
        assert len(state.Arr_Pool) + len(state.Arr_NotReady) > 0

        # Step 4A - Type 1 Repopulation (after flights are released)
        if self.reset_pop:
            self.GA_Info = Populate(state.Ac_Info, self.base_seq, state.Arr_Pool, state.Arr_NotReady, self.S, self.l)
            self.GA_counter = 0 # reset counter
            self.reset_pop = False

        # Steps 4B and 4C - Filter and Repopulate (Type 2)
        self.VNS_counter, self.tot_mut = Repopulate_VNS(self.GA_Info, self.S, self.S_min,
                                                        self.VNS_counter, self.m_mut, self.tot_mut)

        # Step 2A - Simulation and Evaluation
        exp_weather = self.weather.expected_process(state.tm)
        costs = expected_sequence_costs(self.GA_Info, state.tm, state.Ac_Info, state.Ac_queue,
                                        self.sep, exp_weather, state.prev_completion,
                                        state.prev_class, self.tau, self.cost_fn)
        for info, cost in zip(self.GA_Info, costs):
            info.v = cost

        self.GA_Info.sort(key=lambda x: x.v)

        Ac_added = []
        assert len(self.GA_Info) > 0
        seq_info = self.GA_Info[0]
        log.debug("Checking flights for release...")
        log.debug("%s", seq_info)
        log.debug("Status of flights in sequence: %s", ', '.join(str(state.Ac_Info[i].status) for i in seq_info.sequence))
        for AC in seq_info.sequence:
            if state.Ac_Info[AC].status != FlightStatus.IN_POOL:
                break
            Ac_added.append(AC)
        
        # If flights will be removed we need to update base_seq in order to reset sequence population
        if Ac_added:
            self.reset_pop = True
            self.base_seq = self.GA_Info[0].sequence[:]
            for AC in Ac_added:
                self.base_seq.remove(AC)
            # Previously reported - perhaps update this later or remove
            pred_cost = self.GA_Info[0].v
            state.Ac_Info[AC].pred_cost = pred_cost

        return Ac_added

def expected_sequence_costs(GA_Info: list[SequenceInfo], tm: float, Ac_Info: list[FlightInfo], Ac_queue: list[int], 
                            sep: StochasticSeparation,  weather: WeatherProcess, prev_ld: float,
                            prev_class: int,tau: int, cost_fn: Cost):

     # The case where queue is empty needs fixing
     # This function needs to know time of last landing in order to set prev_landing here
    new_cost = 0
    perm_prev_class = prev_class

    for i, AC in enumerate(Ac_queue):
        # Need to generate service times for AC already in the queue; first consider the customer in position 0
        info = Ac_Info[AC]

        # Get the conditional expectation of service time based on service time elapsed so far
        if i == 0:
            sv_time = info.enters_service
            min_sep = sep.expected_conditional_seperation(tm - sv_time, prev_class, info.ac_class,
                                                        info.weather_state)
        else:
            min_sep = sep.expected_separation(prev_class, info.ac_class, info.weather_state)

        landing_complete, _ = landing_time(prev_ld, min_sep, info.eta)
        # JF Note: Do we need conditional travel time below?
        new_cost += cost_fn(info.orig_sched_time, info.pool_time, tau, landing_complete, info.passenger_weight)
        perm_prev_class = info.ac_class
        prev_ld = landing_complete

    stored_prev_ld = prev_ld
    stored_prev_class = perm_prev_class

    seq_ACs = set(ac for seq in GA_Info for ac in seq.sequence)
    pool_times = dict()
    for AC in seq_ACs:
        info = Ac_Info[AC]
        match (status := info.status):
            case FlightStatus.NOT_READY:
                pool_times[AC] = info.eta - tau
            case FlightStatus.IN_POOL:
                pool_times[AC] = info.pool_time
            case _:
                raise RuntimeError(f"Invalid flight status {status}")

    # Try all the sequences in the population
    costs = len(GA_Info)*[new_cost]
    for i, info in enumerate(GA_Info):
        prev_class = stored_prev_class
        prev_ld = stored_prev_ld

        for AC in info.sequence:
            info = Ac_Info[AC]
            pool_time = pool_times[AC]
            reltime = max(tm, pool_time)
            cur_class = info.ac_class
            
            exp_serv = sep.expected_separation(prev_class, cur_class, weather(reltime))
            landing_complete, _ = landing_time(prev_ld, exp_serv, reltime + tau)
            costs[i] += cost_fn(info.orig_sched_time, pool_time, tau, landing_complete, info.passenger_weight)
            
            prev_ld = landing_complete
            prev_class = cur_class

    return costs