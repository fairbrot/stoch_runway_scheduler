from dataclasses import dataclass, field
import logging
from .utils import FlightInfo, FlightStatus, Cost
from .simulate import Update_Stats
from .trajectory import StochasticTrajectory
from .weather import WeatherStatus, WeatherProcess
from .separation import StochasticSeparation

@dataclass
class State:
    Ac_Info: list[FlightInfo]
    "State of all flights"
    tm: float = field(default=0)
    "Current Time (minutes)"
    Arr_NotReady: list[int] = field(default_factory=list)
    "Flights with status NOT_READY"
    Arr_Pool: list[int] = field(default_factory=list)
    "Flights with status IN_POOL"
    Ac_queue: list[int] = field(default_factory=list)
    "Flights with status IN_QUEUE"
    prev_class: int = field(default=4) # initially set to 4 (dummy value)
    "Class of previous aircraft to be served in queue"
    prev_completion: float = field(default=-1e6)
    "Time of latest service completion"
    weather: WeatherStatus = field(default=WeatherStatus.GOOD)
    "Current state of weather"
    # The following attributes can be removed later once scheduling is properly implemented
    real_queue_complete: float = field(default=0)
    "Time that last released flight is scheduled to complete service"
    latest_class: int = field(default = 4)
    "class of latest aircraft to be added to the queue"
    next_completion_time: int = field(default = 0)
    "time next flight is scheduled to land (finish service)"

    def initialise_pool(self, tau):
        print('*** Generating the initial pool...')
        for (i, Ac_Infoi) in enumerate(self.Ac_Info):
            if Ac_Infoi.eta - tau <= 0:
                self.Arr_Pool.append(i)
                print('Aircraft '+str(i)+' initially included in pool (ETA is '+str(Ac_Infoi.eta)+')')
                Ac_Infoi.status = FlightStatus.IN_POOL
            else:
                self.Arr_NotReady.append(i)

    # Note that second arguments won't be necessary once scheduling is implemented
    def release_flights(self, Ac_added: list[int], weather_process: WeatherProcess, sep: StochasticSeparation):
        for AC in Ac_added:
            assert self.Ac_Info[AC].status == FlightStatus.IN_POOL
            self.Arr_Pool.remove(AC)
            self.Ac_queue.append(AC)
            # Gets important statistics about aircraft being released and serviced
            # some of these outputs are used for simulation itself
            self.real_queue_complete, self.next_completion_time, self.latest_class = Update_Stats(self.tm, AC, self.Ac_Info, self.Ac_queue, self.real_queue_complete, 
                                                                                                weather_process, self.latest_class, self.next_completion_time, sep)
    def num_remaining_flights(self) -> int:
        return len(self.Ac_queue) + len(self.Arr_NotReady) + len(self.Arr_Pool)

def round_down(tm: float, freq: int) -> float:
    """Rounds tm down to beginning of interval defined by freq.

    For example, if freq is 100, then we round down for step size of 0.01,
    that is 0.32453 would become 0.32.
    """
    int_size = 1/freq
    i = tm // int_size
    return i*int_size



def Update_ETAs(state: State, trajecs: list[StochasticTrajectory], tau: float, freq: int):
    # JF Question: this only updates ETAs for aircraft not ready, which join the pool, or are in the queue.
    # It does not seemingly update the ETAs for flights aleady in the pool.
    
    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    step_new_logger = logging.getLogger('step_new')

    i = 0
    # Updates flights not aleady in Pool
    to_remove = []
    for i, AC in enumerate(state.Arr_NotReady):
        Ac_Infoi = state.Ac_Info[AC]
        if state.tm >= Ac_Infoi.pool_time:
            state.Arr_Pool.append(AC)
            to_remove.append(i)
            Ac_Infoi.status = FlightStatus.IN_POOL
            Ac_Infoi.eta = Ac_Infoi.pool_time + tau # JF Question: should this be tm + tau?
        else:
            Ac_Infoi.eta = trajecs[AC].expected_eta(state.tm, Ac_Infoi)
    for i in reversed(to_remove):
        state.Arr_NotReady.pop(i)

    # Updates flights which are in the queue
    for i, AC in enumerate(state.Ac_queue):
        Ac_Infoi = state.Ac_Info[AC]
        if not Ac_Infoi.travel_time_indicator:
            rel_time = Ac_Infoi.release_time
            trav_so_far = state.tm - rel_time # amount of time spent travelling to the runway so far
            # JF Question: why round? Should this be round down? Is this related to freq?
            # I think we need rounding down to the beginning of interval defined by freq
            rounded_trav_so_far = round_down(trav_so_far, freq)
            if rounded_trav_so_far >= Ac_Infoi.travel_time:
                Ac_Infoi.travel_time_indicator = True
            else:
                Ac_Infoi.eta = trajecs[AC].expected_eta(state.tm, Ac_Infoi)

def Serv_Completions(state: State, cost_fn: Cost):

    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    # print('Entered Serv_Completions')
    # print('ESC tm: '+str(tm)+' next_completion_time: '+str(next_completion_time))

    arr_cost=0
    dep_cost=0

    j=0

    while len(state.Ac_queue) > 0:

        AC = state.Ac_queue[0]
        Ac_Infoi = state.Ac_Info[AC]
        assert Ac_Infoi.status == FlightStatus.IN_QUEUE

        finish_time = Ac_Infoi.service_completion_time
        current_class = Ac_Infoi.ac_class

        if state.tm >= finish_time: #release_time+trav_time and phase==k:
            #print('* Service phase '+str(phase)+' completed for aircraft '+str(Ac_queue[0])+' at time '+str(tm+delta))
            #print('* Service completion finished for aircraft '+str(AC))
            stepthrough_logger.info('* Service completion finished for aircraft '+str(AC)+'\n'+'\n')
            step_summ_logger.info('* Service completion finished for aircraft '+str(AC)+'\n'+'\n')
            print(f"Flight {AC} finished")
            if Ac_Infoi.status == FlightStatus.IN_QUEUE:
                Ac_Infoi.status = FlightStatus.FINISHED
                arr_cost += cost_fn(Ac_Infoi.orig_sched_time, Ac_Infoi.pool_time, Ac_Infoi.travel_time, finish_time, Ac_Infoi.passenger_weight)

            state.prev_class = current_class

            state.Ac_queue.remove(AC)

            if len(state.Ac_queue) > 0:
                New_AC = state.Ac_queue[0]
                state.next_completion_time = state.Ac_Info[New_AC].service_completion_time

        else:
            break

    return arr_cost, dep_cost
