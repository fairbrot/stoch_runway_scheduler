from dataclasses import dataclass, field
from collections import deque
from .utils import FlightInfo, FlightStatus
from .weather import WeatherStatus


@dataclass
class State:
    # JF Note - FlightInfo objects should be initialised here
    Ac_Info: tuple[FlightInfo]
    "State of all flights"
    tm: float = field(default=0)
    "Current Time (minutes)"
    Arr_NotReady: set[int] = field(init=False)
    "Flights with status NOT_READY"
    Arr_Pool: set[int] = field(init=False)
    "Flights with status IN_POOL"
    Ac_queue: deque[int] = field(default_factory=deque)
    "Flights with status IN_QUEUE"
    prev_class: int = field(default=4) # initially set to 4 (dummy value)
    "Class of previous aircraft to be served in queue"
    prev_completion: float = field(default=-1e6)
    "Time of latest service completion"
    weather: WeatherStatus = field(default=WeatherStatus.GOOD)
    "Current state of weather"

    def __post_init__(self):
        self.Arr_NotReady = set(i for i, f in enumerate(self.Ac_Info) if f.status == FlightStatus.NOT_READY)
        self.Arr_Pool = set(i for i, f in enumerate(self.Ac_Info) if f.status == FlightStatus.IN_POOL)

    def flight_joins_pool(self, tm: float, i: int):
        flight = self.Ac_Info[i]
        flight.status = FlightStatus.IN_POOL
        flight.pool_time = tm
        self.Arr_NotReady.remove(i)
        self.Arr_Pool.add(i)

    def flight_released(self, tm: float, i: int):
        flight = self.Ac_Info[i]
        flight.status = FlightStatus.IN_QUEUE
        flight.weather_state = self.weather
        flight.release_time = tm
        self.Arr_Pool.remove(i)
        self.Ac_queue.append(i)

    def flight_enters_service(self, tm: float, i: int):
        flight = self.Ac_Info[i]
        flight.status = FlightStatus.IN_SERVICE
        flight.enters_service = tm

    def flight_leaves_service(self, tm: float, i: int,
                              service_time: float):
        flight = self.Ac_Info[i]
        flight.status = FlightStatus.FINISHED
        flight.service_completion_time = tm
        flight.service_time = service_time
        self.prev_class = flight.ac_class
        self.prev_completion = tm
        self.Ac_queue.popleft()

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
