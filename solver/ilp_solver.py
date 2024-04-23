import networkx as nx
from ortools.linear_solver import pywraplp
from solver.base_solver import BaseSolver
from itertools import combinations


class ILPSolver(BaseSolver):
    def __init__(self, **kwargs):
        self.solver = pywraplp.Solver.CreateSolver('SCIP')
        self.solver.EnableOutput()
        self.solver.SetNumThreads(1)
        self.solver.SetTimeLimit(60000)

    def solve(self, sample):
        n_variables = sample["n_variables"]
        n_clauses = sample["n_clauses"]
        clauses = sample["clauses"]

        print(f"Number of variables: {n_variables}")
        print(f"Number of clauses: {n_clauses}")
        print(f"Clauses: {clauses}")

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

        print("Graph created with edges: ", g.edges)
        n_nodes = len(g.nodes)
        n_edges = len(g.edges)

        solver = self.solver
        var_map = {}
        d1 = n_nodes - 1
        d2 = n_nodes - 2
        for u in g.nodes:
            var_map[("x", u)] = solver.IntVar(0, 1, "x_{}".format(u))
        for u in g.nodes:
            solver.Add(sum([var_map[("x", v)] for v in g.neighbors(u)]) <= d1 * (1 - var_map[("x", u)]))
            print("sum({}) <= d1 * (1 - {})".format([("x", v) for v in g.neighbors(u)], ("x", u)))
        # for u in g.nodes:
        #     solver.Add(sum([var_map[("x", v)] for v in g.neighbors(u)]) <= 1 + d2 * var_map[("x", u)])
        #     print("sum({}) <= 1 + d2 * {}".format([("x", v) for v in g.neighbors(u)], ("x", u)))
        solver.Maximize(sum([var_map[("x", u)] for u in g.nodes]))
        print("Maximizing: sum({})".format([("x", u) for u in g.nodes]))
        status = solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            for u in g.nodes:
                print("x_{} = ".format(u), var_map[("x", u)].solution_value())
            solution = [u for u in g.nodes if var_map[("x", u)].solution_value() == 1]
            print("Solution: ", solution)
            print("Optimal value found: ", solver.Objective().Value())
            sat_solution_raw = [clauses[u[0]][u[1]] for u in solution]
            sat_solution = []
            for i in range(1, n_variables + 1):
                if -i in sat_solution_raw:
                    sat_solution.append(-i)
                else:
                    sat_solution.append(i)
            print("SAT solution: ", sat_solution)
            return sat_solution



