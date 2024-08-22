import logging
from .utils import Cost, FlightStatus
from .state import State
from .weather import StochasticWeatherProcess
from .sequence import SequenceInfo
from .separation import StochasticSeparation
from .trajectory import StochasticTrajectory
from .simulate import simulate_sequences
from .populate import Populate, Repopulate_VNS


# A logger for this file
log = logging.getLogger(__name__)

class SimHeur:

    def __init__(self, trajecs: list[StochasticTrajectory], sep: StochasticSeparation, weather: StochasticWeatherProcess, cost_fn: Cost, S: int, l:int, n_rel: int, r: int, n_repop: int, S_min: int, m_mut: int):

        # For simulation
        self.trajecs = trajecs
        self.sep = sep
        self.weather = weather

        # Algorithm params
        self.cost_fn = cost_fn
        self.S = S
        self.n_repop = n_repop
        "Number of samples between filter and repopulate (type 2)"
        self.r = r
        "Number of samples between running rank and select"
        self.n_rel = n_rel
        "Minimum number of samples before allowing flights to be released to queue"
        self.S_min = S_min
        self.l = l
        self.m_mut = m_mut
        
        # Internal variables
        self.GA_counter = 0
        self.VNS_counter = 0
        self.tot_mut = 0
        self.reset_pop = True
        self.GA_Info = []

        NoA = len(self.trajecs)
        no_ACs = min(l, NoA)
        self.base_seq = [i for i in range(no_ACs)] # Initial value assumes AcInfo is ordered by Ps_time
        "Base sequence from which a new population is generated"

    def run(self, state: State) -> list[int]:
        # Can't run SimHeur if no flights remain
        assert len(state.Arr_Pool) + len(state.Arr_NotReady) > 0

        # Step 4A - Type 1 Repopulation
        if self.reset_pop:
            self.GA_Info = Populate(state.Ac_Info, self.base_seq, state.Arr_Pool, state.Arr_NotReady, self.S, self.l)
            self.GA_counter = 0 # reset counter
            self.reset_pop = False

        # Steps 4B and 4C - Filter and Repopulate (Type 2)
        if self.GA_counter >= self.n_repop or len(self.GA_Info) < self.S_min:
            self.VNS_counter, self.tot_mut = Repopulate_VNS(self.GA_Info, self.S, self.S_min,
                                                  self.VNS_counter, self.m_mut, self.tot_mut)
            self.GA_counter = 0

        # Step 2A - Simulation and Evaluation
        # JF Note: tm really should be replaced with time last service ended - this might have been an issue with the original code
        costs, xi_lists = simulate_sequences(self.GA_Info, state.tm, state.Ac_Info, 
                                            state.Ac_queue, self.trajecs, self.sep, 
                                            self.weather, state.prev_completion,
                                            state.prev_class, self.cost_fn)
        for info, cost, xi_list in zip(self.GA_Info, costs, xi_lists):
            info.add_observation(cost, xi_list)

        self.GA_Info.sort(key=lambda x: x.v)
        self.GA_counter += 1

        

        # Step 2B - Rank and Select
        if self.GA_counter % self.r == 0:
            log.debug("Best sequence (out of %s) is %s (v = %.2f, var = %.2f, n_traj = %d)", 
                    len(self.GA_Info), self.GA_Info[0].sequence, self.GA_Info[0].v,
                    self.GA_Info[0].var_cost(), self.GA_Info[0].n_traj)
            to_remove = SequenceInfo.rank_and_select(self.GA_Info)
            to_remove.sort(reverse=True)
            log.debug("Rank and select removes sequences with costs: %s",
                     [self.GA_Info[i].v for i in to_remove])
            for i in to_remove:
                info = self.GA_Info.pop(i)

        # Step 2C - Release Flights from Pool
        Ac_added = []
        assert len(self.GA_Info) > 0
        # JF Question: why is second condition necessary here?
        if self.GA_counter % self.n_rel == 0 or (len(self.GA_Info) <= self.S_min):
            seq_info = self.GA_Info[0]
            log.debug("Checking flights for release...")
            log.debug("%s", seq_info)
            log.debug("Status of flights in sequence: %s", ', '.join(str(state.Ac_Info[i].status) for i in seq_info.sequence))
            for (j, AC) in enumerate(seq_info.sequence):
                if seq_info.queue_probs[j] <= 0: # JF Question: in paper this check is only done on first element - modifying this to be the case could allow for simplifications
                    break
                if state.Ac_Info[AC].status != FlightStatus.IN_POOL:
                    break
                #log.debug("Mark flight %d for release ()", j)
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
            log.info("Flights %s marked for release. New base sequence %s", Ac_added, self.base_seq)

        return Ac_added