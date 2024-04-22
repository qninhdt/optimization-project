from typing import Dict, List

import autograd.numpy as np
from autograd import grad
from autograd.misc.optimizers import adam
from . import BaseSolver


class GradientDescentSolver(BaseSolver):
    def __init__(self, weighted: bool = True, learning_rate: float = 0.1) -> None:
        self.learning_rate = learning_rate
        self.weighted = weighted

    def process_idx(self, idx, n):
        idx = np.array(idx)

        mask = idx > 0
        idx[mask] -= 1
        idx[~mask] = -idx[~mask] + n - 1

        return idx

    def solve(self, sample: Dict) -> List[bool]:
        clauses = sample["clauses"]
        m = sample["n_clauses"]
        n = sample["n_variables"]

        if self.weighted:
            # compute weight of variables and clauses
            wv = np.zeros(n)

            for clause in clauses:
                for v in clause:
                    wv[abs(v) - 1] += 1

            wv = wv / np.average(wv)

            wc = np.zeros(m)

            for i, clause in enumerate(clauses):
                for v in clause:
                    wc[i] += wv[abs(v) - 1]

            wc = wc / np.average(wc)
        else:
            wv = np.ones(n)
            wc = np.ones(m)

        np.random.seed(42)
        x = np.random.uniform(0, 1, n)
        solution = (x > 0.5).tolist()

        idx = self.process_idx(clauses, n)

        def loss(x):
            xx = np.concatenate([x, 1 - x])  # (n_variables*2,)
            y = xx[idx]  # (n_clauses*2, 3)

            loss_1 = np.sum(wc * (1 - np.minimum(1, y.sum(axis=1))))
            loss_2 = np.sum((wv * x * (1 - x)) ** 2)

            return loss_1, loss_2

        def loss_(x, _=None):
            loss_1, loss_2 = loss(x)
            return loss_1 + loss_2 * 0.1

        loss_grad = grad(loss_)

        while True:

            last_loss = loss_(x)
            for j in range(100):
                x_grad = loss_grad(x)
                alpha = last_loss / (np.linalg.norm(x_grad) ** 2)

                x -= x_grad * alpha
                solution = (x > 0.5).tolist()
                curr_loss = loss_(x)
                last_loss = curr_loss

            if self.check_solution(sample, solution):
                break

            beta = 0.7
            x = (1 - beta) * x + beta * np.random.uniform(0, 1, n)

        return (x > 0.5).tolist()
