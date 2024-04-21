import os
import requests
import tarfile
import shutil

DATASETS = [
    "uf20-91",
    "uf50-218",
    "uf75-325",
    "uf100-430",
    "uf125-538",
    "uf150-645",
    "uf175-753",
    "uf200-860",
    "uf225-960",
    "uf250-1065",
]

DATASET_URL = "https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/{}.tar.gz"

# remove existing datasets folder
if os.path.exists("datasets"):
    shutil.rmtree("datasets")

os.makedirs("datasets")

# download datasets and extract them
for dataset in DATASETS:
    url = DATASET_URL.format(dataset)
    filename = "datasets/{}.tar.gz".format(dataset)

    print("Downloading {}...".format(dataset))
    r = requests.get(url)
    with open(filename, "wb") as f:
        f.write(r.content)

    print("Extracting {}...".format(dataset))
    with tarfile.open(filename, "r:gz") as tar:
        tar.extractall("datasets/{}".format(dataset))

    current_dir = os.path.join("datasets", dataset)

    if len(os.listdir(current_dir)) == 1:
        first_dir = os.path.join(current_dir, os.listdir(current_dir)[0])
        while len(os.listdir(current_dir)) == 1:
            current_dir = os.path.join(current_dir, os.listdir(current_dir)[0])

        for file in os.listdir(os.path.join(current_dir)):
            shutil.move(
                os.path.join(current_dir, file),
                os.path.join("datasets/{}".format(dataset), file),
            )

        shutil.rmtree(first_dir)

    os.remove(filename)

print("Done!")
