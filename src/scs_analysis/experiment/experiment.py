import subprocess
import time
from typing import List

from cogent3 import make_tree


from scs_analysis.distance.distance import cluster_matching_distance, rooted_rf_distance

data_path = "data/superfine/"
model_trees_end = ".model_tree"
source_trees_end = ".source_trees"

SUP = "SUP"
SCS = "SCS"
MCS = "MCS"
BCD = "BCD"

SCRIPT_PATH = "scripts/method_scripts/"

SCRIPTS = {
    SUP: "run_sup.sh",
    SCS: "run_scs.sh",
    MCS: "run_mcs.sh",
    BCD: "run_bcd.sh",
}

OPTIONS = {}
DEFAULT_OPTIONS = []


def run_methods(
    source_tree_file: str, model_tree_file: str, methods: List[str], verbosity=1
):
    results = {}
    for method in methods:
        if verbosity >= 2:
            print("Running Method", method)

        command = [
            SCRIPT_PATH + SCRIPTS[method],
            source_tree_file,
            *OPTIONS.get(method, DEFAULT_OPTIONS),
        ]
        if verbosity >= 1:
            print(" ".join(command))

        start_time = time.time()
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        end_time = time.time()

        try:
            tree = make_tree(result.stdout.decode("utf-8").strip()).bifurcating()
        except:
            tree = None

        results[method] = (tree, end_time - start_time)

    with open(model_tree_file, "r") as f:
        model = make_tree(f.read().strip()).bifurcating()

    for method in methods:
        tree, time_result = results[method]
        if verbosity >= 1:
            if verbosity >= 2:
                print("Computing distances for", method)
            if tree is None:
                print(f"{method}: time={time_result:.2f}s failed: {source_tree_file}) ")
            else:
                rf = rooted_rf_distance(model, tree)
                mat = cluster_matching_distance(model, tree)
                print(f"{method}: time={time_result:.2f}s RF={rf} MAT={mat}")

    return results


def run_experiment_super_triplets(d: int, k: int, methods: List):
    assert d in [25, 50, 75]
    assert k in [10, 20, 30, 40, 50]

    number_of_experiments = 100

    for i in range(number_of_experiments):
        tree_number = i + 1
        source_file = f"data/SuperTripletsBenchmark/source-trees/d{d}/k{k}/data-d{d}-k{k}-{tree_number}_phybp-s0r.nwk.source_trees"
        model_file = f"data/SuperTripletsBenchmark/model-trees/model-{tree_number}.nwk.model_tree"
        print(f"Results for {i} ({source_file}):")
        run_methods(source_file, model_file, methods)


def run_experiment_smidgen(taxa: int, density: int, methods: List):
    assert taxa in [100, 500, 1000]
    assert density in [20, 50, 75, 100]

    number_of_experiments = 10 if taxa == 1000 else 30

    for i in range(number_of_experiments):
        file = f"data/superfine/{taxa}-taxa/{density}/sm_data.{i}"
        print(f"Results for {i} ({file}):")
        run_methods(file + ".source_trees", file + ".model_tree", methods)


def run_experiment_smidgen_og(
    taxa: int, density: int, methods: List, verbosity: int = 1
):
    assert taxa in [100, 500, 1000, 10000]
    if taxa == 10000:
        assert density == 0
    else:
        assert density in [20, 50, 75, 100]

    number_of_experiments = 30 if taxa != 10000 else 10

    for i in range(number_of_experiments):
        source_file = f"data/SMIDGenOutgrouped/{taxa}/{density}/Source_Trees/RaxML/smo.{i}.sourceTreesOG.tre"
        model_file = f"data/SMIDGenOutgrouped/{taxa}/{density}/Model_Trees/pruned/smo.{i}.modelTree.tre"
        if verbosity >= 1:
            print(f"Results for {i} ({source_file}):")
        run_methods(source_file, model_file, methods, verbosity=verbosity)


if __name__ == "__main__":
    methods = [SCS, BCD, SUP, MCS]
    taxa = 500
    density = 20
    run_experiment_smidgen_og(taxa, density, methods, verbosity=2)
