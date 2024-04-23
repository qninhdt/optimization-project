from solver.ilp_solver import ILPSolver
from solver.test_solver import TestSolver
from dataset import load_dataset

# Load dataset
dataset = load_dataset("uf3-3")

# Initialize solver
solver = ILPSolver()

# Solve dataset
solutions = solver.solve_dataset(dataset)
