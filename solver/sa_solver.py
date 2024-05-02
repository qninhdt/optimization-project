from collections import defaultdict
import networkx as nx
import numpy as np

from solver import BaseSolver
from dwave.samplers import SimulatedAnnealingSampler
from dimod.binary import BinaryQuadraticModel


class SASolver(BaseSolver):
    def __init__(self, **kwargs):
        self.solver = SimulatedAnnealingSampler()

    def solve(self, sample):
        n_variables = sample["n_variables"]
        n_clauses = sample["n_clauses"]
        clauses = sample["clauses"]
        #
        # print(f"Number of variables: {n_variables}")
        # print(f"Number of clauses: {n_clauses}")
        # print(f"Clauses: {clauses}")

        g = nx.Graph()
        for i in range(0, len(clauses)):
            clause = clauses[i]
            for j1 in range(0, 3):
                for j2 in range(j1 + 1, 3):
                    u = (i, j1)
                    v = (i, j2)
                    g.add_edge(u, v)
        for i1 in range(0, len(clauses)):
            for j1 in range(0, 3):
                for i2 in range(i1 + 1, len(clauses)):
                    for j2 in range(0, 3):
                        if clauses[i1][j1] == -clauses[i2][j2]:
                            u = (i1, j1)
                            v = (i2, j2)
                            g.add_edge(u, v)

        # print("Graph created with edges: ", g.edges)
        n_nodes = len(g.nodes)
        n_edges = len(g.edges)

        solver = self.solver
        var_map = {}
        var_cnt = 0
        q = defaultdict(int)
        for u in g.nodes:
            var_map[("x", u)] = var_cnt
            var_cnt += 1

        # Constraint
        # gamma = n_edges // n_nodes + 5
        gamma = 3
        for (u, v) in g.edges:
            q[(var_map[("x", u)], var_map[("x", v)])] += gamma

        # Objective Function
        gamma = 1
        for u in g.nodes:
            q[(var_map[("x", u)], var_map[("x", u)])] -= gamma

        bqm = BinaryQuadraticModel.from_qubo(q)
        response = self.solve_simulated_annealing(bqm)
        sample = response.first.sample
        energy = response.first.energy
        print("Energy: ", energy)
        sample_sum = np.sum([sample[i] for i in range(n_nodes)])
        if sample_sum >= n_clauses:
            solution = [u for u in g.nodes if sample[var_map[("x", u)]] == 1]
            # print("Solution: ", solution)
            # print("Optimal value found: ", sum)
            sat_solution_raw = [clauses[u[0]][u[1]] for u in solution]
            sat_solution = []
            for i in range(1, n_variables + 1):
                if -i in sat_solution_raw:
                    sat_solution.append(False)
                else:
                    sat_solution.append(True)
            # print("SAT solution: ", sat_solution)
            return sat_solution
        return [False] * n_variables

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
