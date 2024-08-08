from typing import List, Tuple, Optional, Protocol
import numpy as np
from .weather import WeatherStatus
from .gamma import sample_cond_gamma, gamma_cond_exp

class StochasticSeparation(Protocol):

    def sample_normalized_separation(self) -> float:
        """Samples a normalized service time which can be scaled by `sample_separation`.

        This is useful for simulating with common random numbers.
        """

    def sample_separation(self, prev_class: int, cur_class: int, 
                          weather_state: WeatherStatus,
                          norm_service_time: Optional[float] = None) -> float:
        """
        Sample minimum separation time between current flight and previous one.

        This is used when flight is not already in service.

        Arguments:
        ---------
        prev_class: weight class of previous flight
        cur_class: weight class of current flight
        weather_state: code for current state of weather (0, 1 or 2)
        norm_service_time: if provided, this value is appropriately scaled and used as the service time for the flight.
                    If not provided, service time is sampled directly from distribution.
        """

    def sample_conditional_separation(self, t_elapsed: float, 
                                      prev_class: int, cur_class: int,
                                      weather_state: WeatherStatus) -> float:
        """
        Sample minimum separation time between current flight and previous one,
        conditional on a given amount of time having elapsed from previous completion.


        Arguments:
        ---------

        t_elapsed: amount of time elapsed since previous completion
        prev_class: weight class of previous flight
        cur_class: weight class of current flight
        weather_state: current state of weather
        """


    def expected_conditional_seperation(self, t_elapsed: float,
                                        prev_class: int, cur_class: int,
                                        weather_state: WeatherStatus) -> float:
        """
        Calculate expected separation time between current flight and previous one,
        conditional on a given amount of time having elapsed from previous completion.


        Arguments:
        ---------

        t_elapsed: amount of time elapsed since previous completion
        prev_class: weight class of previous flight
        cur_class: weight class of current flight
        weather_state: current state of weather
        """

    def expected_separation(self, prev_class: int, cur_class: int, 
                            weather_state: WeatherStatus) -> float:
        """
        Calculate expected minimum separation time between current flight and previous one.

        This is used when flight is not already in service.

        Arguments:
        ---------
        prev_class: weight class of previous flight
        cur_class: weight class of current flight
        weather_state: code for current state of weather (0, 1 or 2)
        """


class ErlangSeparation(StochasticSeparation):

    def __init__(self, Time_Sep: List[List[int]], k: int, w_rho: float):
        """
        Arguments
        ---------
        k: Erlang shape parameter for service time
        Time_Sep: separation times between different classes (minutes)
        w_rho: multiplier for separation time in case of bad weather

        """
        self.Time_Sep = Time_Sep
        self.k = k
        self.w_rho = w_rho

    def sample_normalized_separation(self) -> float:
        """Samples a normalized service time which can be transformed
        into an actual separation time by `sample_separation`.

        This is useful for simulating with common random numbers.
        """
        return np.random.gamma(self.k, 1)


    def sample_separation(self, prev_class: int, cur_class: int, 
                          weather_state: WeatherStatus,
                          norm_service_time: Optional[float] = None) -> float:
        rate = self.k / (self.Time_Sep[prev_class][cur_class]/60)
        if weather_state == WeatherStatus.BAD:
            rate *= 1/self.w_rho
   
        # Transformation causes norm_service_time to go from [mean k, var k] to [mean e_{ij}, var e_{ij}^2/k]
        # See page 8 of paper (Section 3.1)
        return np.random.gamma(self.k, 1/rate) if norm_service_time is None else norm_service_time/rate # service time

    def sample_conditional_separation(self, t_elapsed: float, 
                                      prev_class: int, cur_class: int,
                                      weather_state: WeatherStatus) -> float:
        rate=self.k/(self.Time_Sep[prev_class][cur_class]/60)
        if weather_state==WeatherStatus.BAD:
            rate*=1/self.w_rho

        cond_serv = sample_cond_gamma(t_elapsed, self.k, rate)
        return cond_serv

    def expected_conditional_seperation(self, t_elapsed: float,
                                        prev_class: int, cur_class: int,
                                        weather_state: WeatherStatus) -> float:
        rate = self.k/(self.Time_Sep[prev_class][cur_class]/60)
        if weather_state == WeatherStatus.BAD:
            rate *= 1/self.w_rho
        return gamma_cond_exp(t_elapsed, self.k, rate)

    def expected_separation(self, prev_class: int, cur_class: int, 
                            weather_state: WeatherStatus) -> float:
        mn = self.Time_Sep[prev_class][cur_class]/60
        if weather_state == WeatherStatus.BAD:
            mn *= self.w_rho
        return mn

def landing_time(prev_comp: float, min_sep: float,
                 rel_time: float, trav_time: float) -> Tuple[float, int]:
    """
    Calculates landing time and whether flight went straight into service on reaching runway threshold.

    Arguments
    ---------
    prev_comp: time previous aircraft finished landed
    min_sep: minimum separation time between current and previous aircraft
    rel_time: time current aircraft is released into queue
    trav_time: travel time between from leaving pool to runway threshold

    Returns:
    --------
    t_out: time flight is finished being served
    straight_into_service: indicates whether flight enters service immediately on joining queue
    """
    t1 = rel_time + trav_time
    t2 = prev_comp + min_sep

    if t1 < t2:
        straight_into_service = 0
        t_out = t2
    else:
        straight_into_service = 1
        t_out = t1

    return t_out, straight_into_service
