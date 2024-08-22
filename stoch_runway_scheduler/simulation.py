from dataclasses import dataclass, field
import heapq
import logging
from enum import Enum
import time
from .utils import FlightInfo, FlightStatus, FlightData
from .trajectory import StochasticTrajectory
from .separation import StochasticSeparation, landing_time
from .clock import Clock
from .weather import StochasticWeatherProcess
from .state import State
from .fcfs import ReleasePolicy

# A logger for this file
log = logging.getLogger(__name__)

class EventType(Enum):
    JOINS_POOL = 0
    RELEASE_FLIGHT = 1
    TRAVEL_TIME_COMPLETE = 2
    SERVICE_BEGINS = 3
    SERVICE_ENDS = 4

@dataclass(order=True)
class Event:
    time: float
    type: EventType = field(compare=False)
    i: int = field(compare=False)
    "Index of flight"
    kwargs: dict[str, float] = field(compare=False, default_factory=dict)

class Simulation:
    def __init__(self, flight_data: list[FlightData], ps_time: list[float],
                 pax_weight: list[float], trajecs: list[StochasticTrajectory], 
                 sep: StochasticSeparation, weather_process: StochasticWeatherProcess,
                 tau: float):
        
        assert len(flight_data) == len(ps_time) == len(trajecs) == len(pax_weight)

        self.weather_process = weather_process
        self.flight_data = flight_data
        self.ps_time = ps_time
        self.pax_weight = pax_weight
        self.sep = sep
        self.tau = tau
        self.event_list = []

        self.trajecs = trajecs
        self.norm_service_time = [sep.sample_normalized_separation()
                                    for info in self.flight_data]



    def run_scheduled_events(self, state: State):
        "Process all events in schedule up to current time"
        flg = False
        while self.event_list:
            if self.event_list[0].time > state.tm:
                break
            event = heapq.heappop(self.event_list)
            self.process_event(event, state)
            flg = True
        if flg: # Only log if something has happened
            log.info("Time: %.2f, Arrival Pool: %s, Queue: %s", state.tm, state.Arr_Pool, state.Ac_queue)

    def update_state(self, state: State):
        for (i, flight) in enumerate(state.Ac_Info):
            if flight.status != FlightStatus.FINISHED:
                traj = self.trajecs[i]
                flight.eta = traj.expected_eta(state.tm, flight)
        state.weather = self.weather_process(state.tm)

    def run(self, release_policy: ReleasePolicy, clock: Clock) -> list[FlightInfo]:

        # Initialisation
        self.event_list.clear()
        Ac_Info = []
        for datai, psi, pwi in zip(self.flight_data, self.ps_time, self.pax_weight):
            info = FlightInfo.from_flight_data(datai, psi, pwi)
            Ac_Info.append(info)
        state = State(Ac_Info)

        # Schedule pool joining
        NoA = len(self.flight_data)
        for i in range(NoA):
            traj = self.trajecs[i]
            event = Event(traj.pool_time, EventType.JOINS_POOL, i)
            self.schedule_event(event)
        # Puts flights in pool
        self.run_scheduled_events(state)

        # Start main loop of simulation
        start = time.time()
        num_iter = 0
        while state.num_remaining_flights() > 0:
            assert self.event_list
            Ac_added = release_policy.run(state)
            num_iter += 1
            if Ac_added:
                for i in Ac_added:
                    event = Event(state.tm, EventType.RELEASE_FLIGHT, i)
                    self.schedule_event(event)
                self.run_scheduled_events(state) # Need to update states immediately if flights have been marked for release

            # Calculation new simulation time
            time_elapsed = time.time() - start
            next_event = self.event_list[0].time
            new_tm = clock(state.tm, time_elapsed, num_iter, next_event)
            if new_tm != state.tm:
                num_iter = 0
                state.tm = new_tm
                self.run_scheduled_events(state)
                self.update_state(state)

        return state.Ac_Info


    def schedule_event(self, event: Event):
        heapq.heappush(self.event_list, event)

    def process_event(self, event: Event, state: State):
        assert state.tm >= event.time
        log.info("Processing %s", event)
        etype = event.type
        etime = event.time
        i = event.i
        eflight = state.Ac_Info[i]
        match event.type:
            case EventType.JOINS_POOL:
                state.flight_joins_pool(etime, i)
            case EventType.RELEASE_FLIGHT:
                state.flight_released(etime, i)
                trav_time = self.trajecs[i].travel_time
                trav_event = Event(etime + trav_time, EventType.TRAVEL_TIME_COMPLETE, i)
                self.schedule_event(trav_event)
                if len(state.Ac_queue) == 1:
                    serv_event = Event(etime, EventType.SERVICE_BEGINS, i)
                    self.schedule_event(serv_event)
            case EventType.TRAVEL_TIME_COMPLETE:
                eflight.travel_time = etime - eflight.release_time
            case EventType.SERVICE_BEGINS:
                state.flight_enters_service(etime, i)
                norm_sep = self.norm_service_time[i]
                min_sep = self.sep.sample_separation(state.prev_class,
                                                    eflight.ac_class,
                                                    eflight.weather_state,
                                                    norm_service_time=norm_sep)
                trav_time = self.trajecs[i].travel_time
                serv_fin, _ = landing_time(state.prev_completion, min_sep,
                                        eflight.release_time, trav_time)
                fin_event = Event(serv_fin, EventType.SERVICE_ENDS, i,
                                  {'service_time': min_sep})
                self.schedule_event(fin_event)
            case EventType.SERVICE_ENDS:
                serv_time = event.kwargs['service_time']
                state.flight_leaves_service(etime, i, serv_time)
                log.info("Flight %d lands. Orig. Sched. Arr.: %.2f, PS Time: %.2f, Joined pool: % .2f, Released: %.2f, Travel time: %.2f, Separation: %.2f",
                         i, eflight.orig_sched_time, eflight.ps_time, eflight.pool_time, eflight.release_time, eflight.travel_time, eflight.service_time)
                if len(state.Ac_queue) != 0:
                    flight = state.Ac_queue[0]
                    serv_event = Event(etime, EventType.SERVICE_BEGINS, flight)
                    self.schedule_event(serv_event)
            case _:
                raise RuntimeError('Invalid EventType specified')


