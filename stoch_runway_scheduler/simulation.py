from dataclasses import dataclass, field
import heapq
import logging
from enum import Enum
import time
from .utils import FlightInfo, FlightStatus, FlightData
from .trajectory import StochasticTrajectory
from .separation import StochasticSeparation, landing_time
from .weather import StochasticWeatherProcess
from .state import State
from .sim_heur import SimHeur

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
                 tau: float, release_policy: SimHeur, conv_factor: float, resolution: float):
        
        assert len(flight_data) == len(ps_time) == len(trajecs) == len(pax_weight)
        NoA = len(flight_data)
        self.weather_process = weather_process
        self.sep = sep
        self.tau = tau
        self.release_policy = release_policy
        self.conv_factor = conv_factor
        self.resolution = resolution

        Ac_Info = []
        for datai, psi, pwi in zip(flight_data, ps_time, pax_weight):
            info = FlightInfo.from_flight_data(datai, psi, pwi)
            Ac_Info.append(info)

        # Sort Ac_Info and trajecs according to order of ps_time
        Ac_Info, self.trajecs = zip(*sorted(((Ac_Infoi, traj) for Ac_Infoi, traj in zip(Ac_Info, trajecs)),
                                                 key=lambda tp: tp[0].ps_time))
        self.norm_service_time = [sep.sample_normalized_separation()
                                    for info in Ac_Info]

        self.state = State(Ac_Info)
        # self.state.initialise_pool(tau)

        self.event_list = []
        # Schedule pool joining
        for i in range(NoA):
            traj = self.trajecs[i]
            event = Event(traj.pool_time, EventType.JOINS_POOL, i)
            self.schedule_event(event)
        # Puts flights in pool
        self.run_scheduled_events()

    def run_scheduled_events(self):
        "Process all events in schedule up to current time"
        flg = False
        while self.event_list:
            if self.event_list[0].time > self.state.tm:
                break
            event = heapq.heappop(self.event_list)
            self.process_event(event)
            flg = True
        if flg: # Only log if something has happened
            log.info("Time: %.2f, Arrival Pool: %s, Queue: %s", self.state.tm, self.state.Arr_Pool, self.state.Ac_queue)

    def update_state(self):
        for (i, flight) in enumerate(self.state.Ac_Info):
            if flight.status != FlightStatus.FINISHED:
                traj = self.trajecs[i]
                flight.eta = traj.expected_eta(self.state.tm, flight)
        self.state.weather = self.weather_process(self.state.tm)

    def run(self):
        start = time.time()
        # Puts given flights in pool
        while self.state.num_remaining_flights() > 0:
            Ac_added = self.release_policy.run(self.state)
            for i in Ac_added:
                event = Event(self.state.tm, EventType.RELEASE_FLIGHT, i)
                self.schedule_event(event)
            latest_time = (time.time() - start)/self.conv_factor
            if latest_time - self.state.tm > self.resolution:
                self.state.tm = latest_time
                self.run_scheduled_events()
                self.update_state()




    def schedule_event(self, event: Event):
        heapq.heappush(self.event_list, event)

    def process_event(self, event: Event):
        assert self.state.tm >= event.time
        log.info("Processing %s", event)
        etype = event.type
        etime = event.time
        i = event.i
        eflight = self.state.Ac_Info[i]
        match event.type:
            case EventType.JOINS_POOL:
                self.state.flight_joins_pool(etime, i)
            case EventType.RELEASE_FLIGHT:
                self.state.flight_released(etime, i)
                trav_time = self.trajecs[i].travel_time
                trav_event = Event(etime + trav_time, EventType.TRAVEL_TIME_COMPLETE, i)
                self.schedule_event(trav_event)
                if len(self.state.Ac_queue) == 1:
                    serv_event = Event(etime, EventType.SERVICE_BEGINS, i)
                    self.schedule_event(serv_event)
            case EventType.TRAVEL_TIME_COMPLETE:
                eflight.travel_time = etime - eflight.release_time
            case EventType.SERVICE_BEGINS:
                self.state.flight_enters_service(etime, i)
                norm_sep = self.norm_service_time[i]
                min_sep = self.sep.sample_separation(self.state.prev_class,
                                                    eflight.ac_class,
                                                    eflight.weather_state,
                                                    norm_service_time=norm_sep)
                trav_time = self.trajecs[i].travel_time
                serv_fin, _ = landing_time(self.state.prev_completion, min_sep,
                                        eflight.release_time, trav_time)
                fin_event = Event(serv_fin, EventType.SERVICE_ENDS, i,
                                  {'service_time': min_sep})
                self.schedule_event(fin_event)
            case EventType.SERVICE_ENDS:
                serv_time = event.kwargs['service_time']
                self.state.flight_leaves_service(etime, i, serv_time)
                if len(self.state.Ac_queue) != 0:
                    flight = self.state.Ac_queue[0]
                    serv_event = Event(etime, EventType.SERVICE_BEGINS, flight)
                    self.schedule_event(serv_event)
            case _:
                raise RuntimeError('Invalid EventType specified')


