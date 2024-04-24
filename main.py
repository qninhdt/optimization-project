from solver import GradientDescentSolver, ILPSolver
from dataset import load_dataset

# Load dataset
# dataset = load_dataset("uf3-3")
# dataset = load_dataset("uf20-91")
# dataset = load_dataset("uf50-218")
dataset = load_dataset("uf75-325")

# Initialize solver
solver = ILPSolver()

dataset["samples"] = dataset["samples"][:100]

# Solve dataset
solutions = solver.solve_dataset(dataset)
