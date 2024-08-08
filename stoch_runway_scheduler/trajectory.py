from typing import Protocol
import random
import math
import numpy as np
from .utils import FlightInfo, FlightStatus


def round_down(tm: float, freq: int) -> float:
    """Rounds tm down to beginning of interval defined by freq.

    For example, if freq is 100, then we round down for step size of 0.01,
    that is 0.32453 would become 0.32.
    """
    int_size = 1/freq
    i = tm // int_size
    return i*int_size

class StochasticTrajectory(Protocol):

    def expected_eta(self, tm: float, flight_info) -> float:
        """
        Calculates expected arrival time at runway threshold conditional on time `tm` and current flight status.

        If flight has already reached runway threshold then returns actual arrival time.
        """
    

    def simulate_travel_time(self, tm: float, flight_info: FlightInfo) -> float:
        """
        Simulates travel time between pool and runway threshold conditional at current time and flight status.

        If flight has already ready runway threshold returns actual travel time.
        """

    def simulate_pool_time(self, tm: float, flight_info: FlightInfo) -> float:
        """
        Simulates arrival time at pool conditional on current time and flight status.

        If flight has previously entered pool then returns actual pool arrival time.
        """


class BrownianTrajectory(StochasticTrajectory):
    
    def __init__(self, h: float, Ps_time: float, tau: int, wiener_sig: float, freq: int):
        """
        Generates a trajectory for an aircraft based on Brownian Motion.

        Parameters
        ----------
        h: time at which Brownian motion starts being used to predict ETA
        Ps_time: pre-scheduled arrival time at destination airport plus pre-tactical delay
        tau: flight is considered to have joined pool when current time >= ETA - tau
        wiener_sig: standard deviation of Brownian motion
        freq: frequency (per minute) at which trajectories are updated
        """
        brown_motion = []
        # For flights with scheduled departure less than 0 we need to initially simulate where BM would be at time 0
        if h < 0:
            ETA = random.gauss(Ps_time, math.sqrt(0-h)*wiener_sig) # Update the latest ETAs for ACs that already had their dep time before time zero
        else:
            ETA = Ps_time # ETA = pre-scheduled time
        brown_motion.append(ETA)

        # i.e. when ETA is within 30 minutes of current time (0)
        # JF note: perhaps add a check for aircraft already arriving at runway
        if 0 >= ETA-tau:
            pool_arr_time = 0
            chk = 1
            ETA = tau # JF Question - is this right?
        else:
            chk = 0

        j = 0
        while True:
            j += 1 # step forward in increments of 1/freq minutes
            if j > h*freq: # only update ETA if we've gone beyond the AC's departure time
                # JF Note: this update should be related to freq!
                ETA = random.gauss(ETA, 0.1 * wiener_sig)
            brown_motion.append(ETA)
            if j/freq >= ETA-tau and chk == 0:
                pool_arr_time = round(j/freq, 2) # pool arrival time - j/freq is the 'current time'
                chk = 1
            elif j/freq >= ETA and chk == 1:
                runway_time = round(j/freq, 2)
                # Plane arrives at runway
                travel_time = runway_time - pool_arr_time # travel time between entering pool and arriving at runway
                chk = 2
                break

        self.pool_time: float = pool_arr_time
        "Time at which flight arrives in pool"
        self.travel_time: float = travel_time
        "Time to travel from pool threshold to runway"
        self.brown_motion: list[float] = brown_motion
        "Brownian motion. Element j is the ETA at time j/freq (ignoring pool time)"
        self.freq: int = freq
        "Frequency (per minute) at which ETAs are updated"
        self.wiener_sig: float = wiener_sig
        "standard deviation of Brownian motion"
        self.tau: float = tau
        "Flight is considered to join pool once current time exceeds ETA - tau"

    def expected_eta(self, tm: float, flight_info: FlightInfo) -> float:
        status = flight_info.status
        match status:
            case FlightStatus.NOT_READY:
                return self.brown_motion[int(tm * self.freq)]
            case FlightStatus.IN_POOL:
                return tm + self.tau
            case FlightStatus.IN_QUEUE | FlightStatus.IN_SERVICE:
                if (trav_time := flight_info.travel_time): # in this case travel time to runway complete
                    return flight_info.pool_time + trav_time
                else:
                    rel_time = flight_info.release_time
                    trav_so_far = tm - rel_time # amount of time spent travelling to the runway so far
                    # JF Question: why round? Should this be round down? Is this related to freq?
                    # I think we need rounding down to the beginning of interval defined by freq
                    rounded_trav_so_far = round_down(trav_so_far, self.freq)
                    return self.brown_motion[int((flight_info.pool_time + rounded_trav_so_far) * self.freq)]


    def simulate_travel_time(self, tm: float, flight_info: FlightInfo) -> float:
        if (trav_time := flight_info.travel_time):
            return trav_time
        else:
            return np.random.wald(self.tau, (self.tau/self.wiener_sig)**2)

    def simulate_pool_time(self, tm: float, flight_info: FlightInfo) -> float:
        status = flight_info.status
        match status:
            case FlightStatus.NOT_READY:
                sched = int(round(flight_info.eta - (tm + self.tau), 1)) # JF Question: what is this?
                if sched <= 0:
                    return tm
                else:
                    return np.random.wald(sched, (sched/self.wiener_sig)**2) + tm
            case _:
                return flight_info.pool_time
