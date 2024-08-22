from typing import Protocol
from .state import State

class ReleasePolicy(Protocol):
    """Function object which indicates which flights should be released based on current state."""

    def run(self, state: State) -> list[int]:
        "Returns a list of flight indices to be released"


class FCFS:
    "Function object which marks flights for release as and when they enter the pool"

    def run(self, state: State) -> list[int]:
        return list(state.Arr_Pool)