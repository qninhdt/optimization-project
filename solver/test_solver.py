from typing import Dict, List

import time
import random
from .base_solver import BaseSolver


class TestSolver(BaseSolver):
    """
    Test solver for the SAT problem
    """

    def __init__(self) -> None:
        pass

    def solve(self, sample: Dict) -> List[bool]:
        solution = [False] * sample["n_variables"]

        # random sleep to simulate computation
        time.sleep(random.random())

        for i in range(sample["n_variables"]):
            solution[i] = True

        return solution
