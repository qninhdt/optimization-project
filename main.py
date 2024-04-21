from solver.test_solver import TestSolver
from dataset import load_dataset

# Load dataset
dataset = load_dataset("uf3-3")

# Initialize solver
solver = TestSolver()

# Solve dataset
solutions = solver.solve_dataset(dataset)
