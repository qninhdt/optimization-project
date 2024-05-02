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
        self.q = defaultdict(int)
        self.var_map = {}
        self.var_cnt = 0
        self.offset = 0
        # self.solver = ()

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
        response = self.solve_quantum_annealing(bqm)
        sample = response.first.sample
        energy = response.first.energy
        print("Energy: ", energy + self.offset)
        sat_solution = [bool(sample[self.var_map[("x", i + 1)]] == 1) for i in range(n_variables)]
        # print("SAT solution: ", sat_solution)
        return sat_solution

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
