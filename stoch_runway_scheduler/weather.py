from dataclasses import dataclass
from typing import Protocol
from enum import Enum
import random
import numpy as np

class WeatherStatus(Enum):
    GOOD = 0
    BAD = 1

class WeatherProcess(Protocol):
    def __call__(self, tm: float) -> WeatherStatus:
        ...

class StochasticWeatherProcess(WeatherProcess, Protocol):
    def sample_process(self, tm: float) -> WeatherProcess:
        ...

    def expected_process(self, tm: float) -> WeatherProcess:
        ...

@dataclass
class SimpleWeatherProcess:
    wlb: float
    wub: float

    def __call__(self, tm: float) -> WeatherStatus:
        # if self.wlb >= self.wub: # Case not explicitly needed
        #     return WeatherStatus.GOOD
        
        if tm >= self.wlb and tm <= self.wub:
            return WeatherStatus.BAD
        else:
            return WeatherStatus.GOOD

class BrownianWeatherProcess(StochasticWeatherProcess, SimpleWeatherProcess):

    def __init__(self, wlb: int, wub: int, T: int, weather_sig: float, freq: int):
        """
        Generates brownian motion for bad weather forecasts, as well as actual times of bad weather.

        Parameters
        ----------
        wlb: prediction at time 0 for start of bad weather
        wub: prediction at time 0 for end of bad weather
        T: length of time horizon in minutes (used to pre-allocate output arrays)
        weather_sig: standard deviation of Brownian motion
        freq: number of updates for each minute
        """

        N = T * freq # Size of weather array
        weather_lb = [0] * N # Dynamically forcast for start of bad weather T_0(t)
        weather_ub = [0] * N # Dynamically forcast for end of bad weather T_1(t)

        # called U_0 and U_1 in the paper
        wlb_tm = 0 # Actual (randomly generated) time at which bad weather starts; leave this as zero
        wub_tm = 0 # Actual (randomly generated) time at which bad weather ends; leave this as zero

        weather_lb[0] = wlb # Set initial forecasted values
        weather_ub[0] = wub
        old_lb = wlb
        old_ub = wub

        j = 0
        # Brownian motion again used to get dynamic forecast prediction
        for j in range(N):
            new_lb = random.gauss(old_lb, 0.1 * weather_sig) # random.gauss(old_lb,0.01) JF Question: why is 0.1 inside here - shouldn't this be incorporated into weather_sig
            new_ub = random.gauss(old_ub, 0.1 * weather_sig)
            if new_lb > new_ub: # JF Question: Why is this needed?
                new_ub = new_lb
            if j/freq >= new_lb and wlb_tm == 0:
                wlb_tm = j/freq
            if j/freq >= new_ub and wub_tm == 0:
                wub_tm = j/freq
            weather_lb[j] = new_lb
            weather_ub[j] = new_ub
            old_lb = new_lb
            old_ub = new_ub

        self.weather_sig = weather_sig
        self.freq = freq
        self.wlb = wlb_tm # actual start time of bad weather
        self.wub = wub_tm # actual end time of bad weather
        # forecast for start of bad weather over time
        # Each element i corresponds to forecast at time j*freq.
        self.weather_lb = weather_lb 
        # forecast for end of bad weather over time.
        self.weather_ub = weather_ub

    def sample_process(self, tm: float) -> SimpleWeatherProcess:
        wlb = self.weather_lb[int(tm * self.freq)]
        wub = self.weather_ub[int(tm * self.freq)]

        chk=0
        while chk==0:
            if tm >= wlb:
                wlb_gen = wlb
            else:
                # Do wlb_gen
                sched = int(round(wlb - tm, 1))
                if sched<=0:
                    wlb_gen=tm
                else:
                    wlb_gen=np.random.wald(sched, (sched/self.weather_sig)**2) + tm
            if tm >= wub:
                wub_gen = wub
            else:
                # Do wub_gen
                sched = int(round(wub-tm,1))
                if sched <= 0:
                    wub_gen = tm
                else:
                    wub_gen = np.random.wald(sched,(sched / self.weather_sig)**2) + tm
            if wlb_gen <= wub_gen: # JF Question: this if else and chk seem pointless?
                chk = 1
            else:
                chk = 1
        return SimpleWeatherProcess(wlb_gen, wub_gen)

    def expected_process(self, tm) -> SimpleWeatherProcess:
        wlb = self.wlb if tm >= self.wlb else self.weather_lb[int(tm * self.freq)]
        wub = self.wub if tm >= self.wub else self.weather_ub[int(tm * self.freq)]
        return SimpleWeatherProcess(wlb, wub)