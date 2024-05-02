from collections import defaultdict
import networkx as nx
import numpy as np

from solver import BaseSolver
from dwave.samplers import SimulatedAnnealingSampler
from dimod.binary import BinaryQuadraticModel


class SASolver(BaseSolver):
    def __init__(self, **kwargs):
        self.solver = SimulatedAnnealingSampler()
        self.q = defaultdict(int)
        self.var_map = {}
        self.var_cnt = 0
        self.offset = 0

    def solve(self, sample):
        self.q = defaultdict(int)
        self.var_map = {}
        self.var_cnt = 0
        self.offset = 0

        n_variables = sample["n_variables"]
        n_clauses = sample["n_clauses"]
        clauses = sample["clauses"]
        #
        # print(f"Number of variables: {n_variables}")
        # print(f"Number of clauses: {n_clauses}")
        # print(f"Clauses: {clauses}")

        for i in range(n_variables):
            self.var_map[("x", i + 1)] = i
            self.var_cnt += 1
        for clause in clauses:
            self.add_clause(*clause)

        bqm = BinaryQuadraticModel.from_qubo(self.q)
        response = self.solve_simulated_annealing(bqm)
        sample = response.first.sample
        energy = response.first.energy
        print("Energy: ", energy + self.offset)
        sat_solution = [bool(sample[self.var_map[("x", i + 1)]] == 1) for i in range(n_variables)]
        # print("SAT solution: ", sat_solution)
        return sat_solution

    def solve_simulated_annealing(self, bqm, method="?_", num_reads=1000):
        sampler = SimulatedAnnealingSampler()
        # beta_range = [0.1, 4]
        num_sweeps = 100
        # config = method + str(num_reads) + "-SA" + "".join(str(beta_range).split(" ")) + "s" + str(num_sweeps)
        solver_config = method + str(num_reads) + "-SA" + "s" + str(num_sweeps)
        # os.environ["SOLVER_CONFIG"] = solver_config
        # data_name = os.getenv("PHYLO_FILE")
        # output_dir = "output/" + data_name + "/"
        # if not os.path.exists(output_dir):
        #     os.makedirs(output_dir)
        # print(config)
        # start = time.time()
        response = sampler.sample(bqm,
                                  num_reads=num_reads,
                                  label=solver_config,
                                  # beta_range=beta_range,
                                  num_sweeps=num_sweeps)
        # end = time.time()
        # config_dict = {
        #     "config": solver_config,
        #     "num_vars": len(bqm.variables),
        #     "time_elapsed": end - start,
        # }
        # with open(output_dir + solver_config + "_1.json", "w") as f:
        #     json.dump(config_dict, f, indent=4)
        return response

    def add_clause(self, a, b, c, gamma=1):
        val_a = abs(a)
        val_b = abs(b)
        val_c = abs(c)
        self.var_map[("s", self.var_cnt)] = self.var_cnt
        if a > 0 and b > 0 and c > 0:
            # - c + ab - as - bs - cs
            self.q[self.var_map[("x", val_c)], self.var_map[("x", val_c)]] -= gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_b)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] += gamma
            self.offset += gamma
        elif a > 0 and b > 0 and c < 0:
            # + c + ab - as - bs - cs + s
            self.q[self.var_map[("x", val_c)], self.var_map[("x", val_c)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_b)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += gamma
        elif a > 0 and b < 0 and c > 0:
            # + b - bc - as - bs + cs + s
            self.q[self.var_map[("x", val_b)], self.var_map[("x", val_b)]] += gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("x", val_c)]] -= gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] += gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += gamma
        elif a < 0 and b > 0 and c > 0:
            # + a - ac - as - bs + cs + s
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_a)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_c)]] -= gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] += gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += gamma
        elif a > 0 and b < 0 and c < 0:
            # +bc - as - bs - cs + 2s
            self.q[self.var_map[("x", val_b)], self.var_map[("x", val_c)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += 2 * gamma
        elif a < 0 and b > 0 and c < 0:
            # +ac - as - bs - cs + 2s
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_c)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += 2 * gamma
        elif a < 0 and b < 0 and c > 0:
            # +ab - as - bs - cs + 2s
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_b)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += 2 * gamma
        elif a < 0 and b < 0 and c < 0:
            # +ab - as - bs + cs + s
            self.q[self.var_map[("x", val_a)], self.var_map[("x", val_b)]] += gamma
            self.q[self.var_map[("x", val_a)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_b)], self.var_map[("s", self.var_cnt)]] -= gamma
            self.q[self.var_map[("x", val_c)], self.var_map[("s", self.var_cnt)]] += gamma
            self.q[self.var_map[("s", self.var_cnt)], self.var_map[("s", self.var_cnt)]] += gamma
        self.var_cnt += 1
