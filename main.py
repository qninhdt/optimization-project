from solver import GradientDescentSolver
from dataset import load_dataset
import json
import os

# Load dataset
DATASETS = [
    "uf20-91",
    "uf50-218",
    "uf75-325",
]

MAX_SAMPLES = 100

# create folder reports if it does not exist
os.makedirs("reports", exist_ok=True)


def run_experiment(solver, dataset_name: str):
    print(f"> Running experiment with {type(solver).__name__} on {dataset_name}")
    os.makedirs(f"reports/{type(solver).__name__}", exist_ok=True)

    dataset = load_dataset(dataset_name)

    n_samples = min(MAX_SAMPLES, len(dataset["samples"]))
    dataset["samples"] = dataset["samples"][:n_samples]
    dataset["n_samples"] = n_samples

    report = solver.solve_dataset(dataset)

    with open(f"reports/{type(solver).__name__}/{dataset_name}.json", "w") as f:
        json.dump(report, f, indent=4)


for dataset_name in DATASETS:
    solver = GradientDescentSolver()
    run_experiment(solver, dataset_name)
