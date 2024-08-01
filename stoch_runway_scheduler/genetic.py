from .utils import FlightInfo, Cost, FlightStatus
from .weather import StochasticWeatherProcess
from .sequence import SequenceInfo
from .separation import StochasticSeparation
from .trajectory import StochasticTrajectory
from .simulate import simulate_sequences
from .populate import Populate, Repopulate_VNS

# Params internal to Genetic
# GA_Checksize is reset to GA_Check_Increment at beginning of each full run
# It is also reset to GA_Check_Increment when a flight is released to the queue
# and when set of sequences is repopulated. It is increased or reset after
# running rank_and_select in Genetic
# 
# GA_Check_Increment is initialised to GA_LoopSize / 10 at the beginning of every
# full run but then isn't changed.
# GA_Check_Increment used for releasing flights

#
# GA_LoopSize is initialized at beginning of each run depending on what
# policy is being used but then never changes
# 
# GA_Counter is incremented every time Genetic is run (inside Genetic).
# It is initialised to 0 at beginning of every full run.
# It is reset to zero after running Populate and Repopulate
# It is also incremented by 1 after Genetic_determ is run
#
# basecost is an input which isn't changed by Genetic
#
# tau is a constant input used for simulating sequences
# 
# Max_Lookahead is no longer used
# 
# cost_fn is a fixed input
# 
# S_min is fixed as is wiener_sig

# Outputs

# Ac_added is list of flights to be released to queue
# counter is number of trajectories sampled for a sequence
# from which flights are released. 0 if no flights are released.
# It doesn't seem to be used anywhere

class SimHeur:

    def __init__(self, Ac_Info: list[FlightInfo], Arr_NotReady: list[int], Arr_Pool: list[int], Ac_queue: list[int], trajecs: list[StochasticTrajectory], sep: StochasticSeparation, weather: StochasticWeatherProcess, cost_fn: Cost, GA_PopSize: int, Max_LookAhead:int, n_rel: int, r: int, n_repop: int, S_min: int, VNS_limit: int):
        # State
        self.Ac_Info = Ac_Info
        self.Arr_NotReady = Arr_NotReady
        self.Arr_Pool = Arr_Pool
        self.Ac_queue = Ac_queue
        self.trajecs = trajecs
        # For simulation
        self.sep = sep
        self.weather = weather
        self.cost_fn = cost_fn
        # Algorithm params
        self.GA_PopSize = GA_PopSize
        self.n_repop = n_repop
        "Number of samples between filter and repopulate (type 2)"
        self.r = r
        "Number of samples between running rank and select"
        self.n_rel = n_rel
        "Minimum number of samples before allowing flights to be released to queue"
        self.S_min = S_min
        self.Max_LookAhead = Max_LookAhead
        self.VNS_limit = VNS_limit
        # Internal variables
        self.GA_counter = 0
        self.VNS_counter = 0
        self.tot_mut = 0
        self.reset_pop = True
        self.GA_Info = []

        NoA = len(self.Ac_Info)
        no_ACs = min(Max_LookAhead, NoA)
        self.base_seq = [i for i in range(no_ACs)] # Initial value assumes AcInfo is ordered by Ps_time
        "Base sequence from which a new population is generated"

    def run(self, tm: float, prev_class: int, basecost: float) -> list[int]:
        # Can't run SimHeur if no flights remain
        assert len(self.Arr_Pool) + len(self.Arr_NotReady) > 0

        # Step 4A - Type 1 Repopulation
        if self.reset_pop:
            self.GA_Info = Populate(self.Ac_Info, self.base_seq, self.Arr_Pool, self.Arr_NotReady, self.GA_PopSize, self.Max_LookAhead)
            self.GA_counter = 0 # reset counter
            self.reset_pop = False

        # Steps 4B and 4C - Filter and Repopulate (Type 2)
        if self.GA_counter >= self.n_repop or len(self.GA_Info) < self.S_min:
            self.VNS_counter, self.tot_mut = Repopulate_VNS(self.GA_Info, self.GA_PopSize, self.S_min,
                                                  self.VNS_counter, self.VNS_limit, self.tot_mut)
            self.GA_counter = 0

        # Step 2A - Simulation and Evaluation
        # JF Note: tm really should be replaced with time last service ended - this might have been an issue with the original code
        costs, xi_lists = simulate_sequences(self.GA_Info, tm, self.Ac_Info, self.Ac_queue, self.trajecs, self.sep, self.weather, prev_class, self.cost_fn)
        for info, cost, xi_list in zip(self.GA_Info, costs, xi_lists):
            info.add_observation(cost + basecost, xi_list)

        self.GA_Info.sort(key=lambda x: x.v)
        self.GA_counter += 1

        # Step 2B - Rank and Select
        if self.GA_counter % self.r == 0:
            to_remove = SequenceInfo.rank_and_select(self.GA_Info)
            to_remove.sort(reverse=True)
            for i in to_remove:
                info = self.GA_Info.pop(i)

        # Step 2C - Release Flights from Pool
        Ac_added = []
        assert len(self.GA_Info) > 0
        # JF Question: why is second condition necessary here?
        if self.GA_counter >= self.n_rel or (len(self.GA_Info) <= self.S_min):
            perm = self.GA_Info[0]
            for (j, AC) in enumerate(perm.sequence):
                if perm.queue_probs[j] <= 0: # JF Question: in paper this check is only done on first element - modifying this to be the case could allow for simplifications
                    break
                if self.Ac_Info[AC].status != FlightStatus.IN_POOL:
                    break
                Ac_added.append(AC)

        # If flights will be removed we need to update base_seq in order to reset sequence population
        if Ac_added:
            self.reset_pop = True
            self.GA_Info.sort(key=lambda x: x.v)
            self.base_seq = self.GA_Info[0].sequence[:]
            for AC in Ac_added:
                self.base_seq.remove(AC)
            # Previously reported - perhaps update this later or remove
            pred_cost = self.GA_Info[0].v
            self.Ac_Info[AC].pred_cost = pred_cost

        return Ac_added