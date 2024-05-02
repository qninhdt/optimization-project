from collections import defaultdict
import networkx as nx
import numpy as np
from dwave.embedding.chain_strength import uniform_torque_compensation
from dwave.system import EmbeddingComposite, DWaveSampler
import dwave.inspector

from solver import BaseSolver
from dimod.binary import BinaryQuadraticModel


class QASolver(BaseSolver):
    def __init__(self, **kwargs):
        pass
        # self.solver = ()

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

        # solver = self.solver
        var_map = {}
        var_cnt = 0
        q = defaultdict(int)
        for u in g.nodes:
            var_map[("x", u)] = var_cnt
            var_cnt += 1

        # Constraint
        # gamma = n_edges // n_nodes + 5
        gamma = 2
        for (u, v) in g.edges:
            q[(var_map[("x", u)], var_map[("x", v)])] += gamma

        # Objective Function
        gamma = 1
        for u in g.nodes:
            q[(var_map[("x", u)], var_map[("x", u)])] -= gamma

        bqm = BinaryQuadraticModel.from_qubo(q)
        response = self.solve_quantum_annealing(bqm)
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

    def solve_quantum_annealing(self, bqm,
                                method="?_",
                                num_reads=1000):
        chain_strength_prefactor = 0.25
        annealing_time = 200
        anneal_schedule_id = -1
        chain_strength = uniform_torque_compensation(
            bqm=bqm, prefactor=chain_strength_prefactor)
        schedules = {12: [(0.0, 0.0), (40.0, 0.4), (180.0, 0.4), (200.0, 1.0)],
                     11: [(0.0, 0.0), (40.0, 0.5), (120.0, 0.5), (200.0, 1.0)],
                     13: [(0.0, 0.0), (40.0, 0.5), (130.0, 0.5), (200.0, 1.0)],
                     14: [(0.0, 0.0), (30.0, 0.5), (160.0, 0.5), (200.0, 1.0)]}
        embed_config = method + str(num_reads) + "-" + str(chain_strength)
        solver_config = method + str(num_reads) + "-" + str(chain_strength) + "s" + str(
            anneal_schedule_id) + "_A" + str(
            annealing_time)
        # os.environ["SOLVER_CONFIG"] = solver_config
        # data_name = os.getenv("PHYLO_FILE")
        # data_name = ""
        # output_dir = "output/" + data_name + "/"
        # if not os.path.exists(output_dir):
        #     os.makedirs(output_dir)
        # get_embedding(bqm, output_dir + embed_config + ".json")
        # with open(output_dir + embed_config + ".json", "r") as f:
        #     embedding = json.load(f)
        # embedding = {int(k): v for k, v in embedding.items()}
        # try:
        #     sampler = EmbeddingComposite(DWaveSampler(), embedding=embedding)
        # except DisconnectedChainError as e:
        #     raise RuntimeError(e)
        sampler = EmbeddingComposite(DWaveSampler())

        # start = time.time()
        if anneal_schedule_id == -1:
            response = sampler.sample(bqm=bqm,
                                      chain_strength=chain_strength,
                                      num_reads=num_reads,
                                      label=solver_config,
                                      annealing_time=annealing_time)
        else:
            schedule = schedules[anneal_schedule_id]
            response = sampler.sample(bqm=bqm,
                                      chain_strength=chain_strength,
                                      num_reads=num_reads,
                                      label=solver_config,
                                      # annealing_time=annealing_time,
                                      anneal_schedule=schedule)
        # end = time.time()
        # dwave.inspector.show(response)
        # chains = response.info["embedding_context"]["embedding"].values()
        # for chain in chains:
        #     if len(chain) > 10:
        #         print(chain)

        # config_dict = {
        #     "config": solver_config,
        #     "num_vars": len(bqm.variables),
        #     "num_qubit": sum([len(chain) for chain in chains]),
        #     "time_elapsed": end - start,
        #     "best_state": {
        #         "sample": response.record.sample[0].tolist(),
        #         "energy": response.record.energy[0],
        #         "chain_break_fraction": response.record.chain_break_fraction[0],
        #     },
        #     "chain_strength_prefactor": chain_strength_prefactor,
        #     "chain_strength": chain_strength,
        #     "max_chain_length": max([len(chain) for chain in chains]),
        #     "timing_info": response.info["timing"],
        #     "embedding_info": response.info["embedding_context"]
        # }
        # with open(output_dir + solver_config + "_1.json", "w") as f:
        #     json.dump(config_dict, f, indent=4)
        # file_response_output = open(output_dir + solver_config + ".txt", "w")
        # file_response_output.write("Config: " + solver_config + "\n")
        # file_response_output.write("Number of source variables: " + str(len(bqm.variables)) + "\n")
        # file_response_output.write("Number of target variables: " + str(sum([len(chain) for chain in chains])) + "\n")
        # file_response_output.write("Time Elapsed: " + str(end - start) + "\n")
        # file_response_output.write(
        #     "Best State: " + str(response.record.sample[0]) + "\n" + str(response.record.energy[0]) + "\t" + str(
        #         response.record.chain_break_fraction[0]) + "\n")
        # file_response_output.write("ChainStr/ChainLen: " + str(chain_strength) + "/" + str(
        #     max([len(chain) for chain in chains])) + "\n")
        # file_response_output.write("Info: " + str(response.info["timing"]) + "\n")
        # file_response_output.write("Embedding Info: " + str(response.info["embedding_context"]) + "\n")

        dwave.inspector.show(response)
        return response
