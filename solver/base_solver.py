from typing import Dict, List

from timeit import default_timer as timer
import colorful as cf


class BaseSolver:
    """
    Base class for a SAT solver
    """

    def __init__(self) -> None:
        pass

    def solve(self, sample: Dict) -> List[bool]:
        """
        Skeleton method for solving a SAT problem

        Args:
            sample (dict): Sample data

        Returns:
            List[bool]: List of boolean values representing the solution
                        index i represents the value of variable i+1
        """
        raise NotImplementedError

    def count_satisfied_clauses(self, sample: Dict, solution: List[bool]) -> int:
        """
        Count the number of satisfied clauses in a sample

        Args:
            sample (dict): Sample data
            solution (List[bool]): List of boolean values representing the solution
                                   index i represents the value of variable i+1

        Returns:
            int: Number of satisfied clauses
        """
        n_satisfied = 0
        for clause in sample["clauses"]:
            if any(
                [
                    solution[abs(v) - 1] if v > 0 else not solution[abs(v) - 1]
                    for v in clause
                ]
            ):
                n_satisfied += 1
        return n_satisfied

    def check_solution(self, sample: Dict, solution: List[bool]) -> bool:
        """
        Check if a solution is valid

        Args:
            sample (dict): Sample data
            solution (List[bool]): List of boolean values representing the solution
                                   index i represents the value of variable i+1

        Returns:
            bool: True if the solution is valid, False otherwise
        """
        for clause in sample["clauses"]:
            if not any(
                [
                    solution[abs(v) - 1] if v > 0 else not solution[abs(v) - 1]
                    for v in clause
                ]
            ):
                return False

        return True

    def solve_dataset(self, dataset: Dict) -> List[List[bool]]:
        print(f"> Solving dataset {dataset['name']}")
        solutions = []
        durations = 0

        solved = 0
        for sample in dataset["samples"]:

            start = timer()
            solution = self.solve(sample)
            end = timer()

            duration = end - start
            durations += duration

            n_valid_clauses = self.count_satisfied_clauses(sample, solution)

            if n_valid_clauses == sample["n_clauses"]:
                print(
                    cf.green
                    | f"[{sample['id']}/{dataset['n_samples']}] Solved in {duration:.4f}s ({n_valid_clauses}/{sample['n_clauses']})"
                )
                solved += 1
            else:
                print(
                    cf.red
                    | f"[{sample['id']}/{dataset['n_samples']}] Not solved in {duration:.4f}s ({n_valid_clauses}/{sample['n_clauses']})"
                )

            solutions.append(solution)

        print(
            f"> Solved {solved}/{dataset['n_samples']} samples in dataset {dataset['name']} in {durations:.4f}s, average {durations/dataset['n_samples']:.4f}s per sample"
        )
