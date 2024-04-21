from typing import Dict, List

import os
from utils import clean_str


def load_sample(sample_path: str) -> Dict:
    """
    Load a sample from a file

    Args:
        sample_path (str): Path to the sample file

    Returns:
        dict: Sample data

    Sample data:
    {
        "id": int,
        "n_clauses": int,
        "n_variables": int,
        "clauses": List[List[int]]
    }
    """
    with open(sample_path, "r") as f:
        data = f.read()

    # extract filename from path
    sample_name = os.path.basename(sample_path)

    # extract sample id from filename
    sample_id = int(sample_name.split("-")[1].split(".")[0])

    n_clauses = None
    n_variables = None
    clauses = []

    for line in data.split("\n"):
        line = clean_str(line)
        d = line.split(" ")

        if line.startswith("p cnf"):
            n_variables = int(d[2])
            n_clauses = int(d[3])
        elif (
            line.startswith("c")
            or line.startswith("%")
            or line.startswith("0")
            or line == ""
        ):
            continue
        else:
            v1, v2, v3, _ = d
            v1 = int(v1)
            v2 = int(v2)
            v3 = int(v3)
            clauses.append([v1, v2, v3])

    assert n_clauses is not None, "Number of clauses not found"
    assert n_variables is not None, "Number of variables not found"
    assert len(clauses) == n_clauses, "Number of clauses does not match"

    sample = {
        "id": sample_id,
        "n_clauses": n_clauses,
        "n_variables": n_variables,
        "clauses": clauses,
    }

    return sample


def load_dataset(dataset_name: str) -> Dict:
    """
    Load a dataset

    Args:
        dataset_name (str): Name of the dataset, e.g. "uf20-91"

    Returns:
        dict: Dataset data

    Dataset data:
    {
        "name": str,
        "n_clauses": int,
        "n_variables": int,
        "n_samples": int,
        "samples": List[dict]
    }
    """
    print(f"> Loading dataset {dataset_name}")

    dataset_path = os.path.join("datasets", dataset_name)
    dataset_files = os.listdir(dataset_path)
    dataset_files = [os.path.join(dataset_path, file) for file in dataset_files]

    n_clauses = int(dataset_name.split("-")[0][2:])
    n_variables = int(dataset_name.split("-")[1].split(".")[0])
    samples = [load_sample(file) for file in dataset_files]

    # sort samples by id
    samples = sorted(samples, key=lambda x: x["id"])

    dataset = {
        "n_clauses": n_clauses,
        "n_variables": n_variables,
        "n_samples": len(samples),
        "samples": samples,
        "name": dataset_name,
    }

    return dataset


def visualize_sample(sample: Dict) -> None:

    def to_sub(n: int) -> str:
        N = "₀₁₂₃₄₅₆₇₈₉"
        return "".join([N[int(i)] for i in str(n)])

    # print sample data with ∨ and ¬ symbols

    s = " ⋀ ".join(
        [
            "("
            + " ⋁ ".join([f"{'x̄' if v < 0 else 'x'}{to_sub(abs(v))}" for v in clause])
            + ")"
            for clause in sample["clauses"]
        ]
    )

    print(s)


if __name__ == "__main__":
    dataset = load_dataset("uf3-3")

    visualize_sample(dataset["samples"][0])
