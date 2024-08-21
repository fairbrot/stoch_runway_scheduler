from typing import Protocol

class Clock(Protocol):
    "Function object which calculates what time clock should be moved to in simulation"

    def __call__(self, tm: float, time_elapsed: float, num_iter: int, next_event: float) -> float:
        """
        Calculates time clock should be updated to given current time `tm`, the
        computational time elapsed since the simulation began `time_elapsed`, the
        number of times the release policy has been run since last time update `num_iter`,
        and the scheduled time of the next event, `next_event`.
        """

class ComputationalClock:
    "Clock which moves simulation time forward based on how much computational time has elapsed."

    def __init__(self, conv_factor: int, resolution: float):
        """
        Initialises ComputationalClock object.

        Arguments
        ----------
        conv_factor: conversion factor between computation time and simulation time
                     More precisely, simulation time is computational time elapsed / conv_factor
        resolution: Minimal time by which clock can be incremented.
        """
        self.conv_factor = conv_factor
        self.resolution = resolution

    def __call__(self, tm: float, time_elapsed: float, num_iter: int, next_event: float) -> float:
        latest_tm = time_elapsed/self.conv_factor
        return latest_tm if (latest_tm - tm) > self.resolution else tm


class EventClock:
    "Clock which moves simulation time to next event"
    def __init__(self):
        pass

    def __call__(self, tm: float, time_elapsed: float, num_iter: int, next_event: float) -> float:
        return next_event


class ComputationalBudgetClock:
    "Clock which moves forward in fixed increments at regular numbers of iterations of simulation."
    def __init__(self, budget: int, increment: float):
        """
        Initialises ComputationalBudgetClock.

        Arguments:
        ----------
        budget: number of iterations of simulation between updating clock
        increment: increment to time when clock is updated
        """
        self.budget = budget
        self.increment = increment

    def __call__(self, tm: float, time_elapsed: float, num_iter: int, next_event: float) -> float:
        return tm + self.increment if num_iter >= self.budget else tm



    