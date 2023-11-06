import subprocess
import time
import os
from typing import List, Optional
from cogent3 import make_tree
from cogent3.core.tree import TreeNode
from scs_analysis.distance.distance import matching_cluster_distance, rooted_rf_distance

SUPER_TRIPLET_D = (25, 50, 75)
SUPER_TRIPLET_K = (10, 20, 30, 40, 50)

SMIDGEN_TAXA = (100, 500, 1000)
SMIDGEN_DENSITY = (20, 50, 75, 100)

SMIDGEN_OG_NORMAL_TAXA = (100, 500, 1000, 10000)
SMIDGEN_OG_DENSITY = (20, 50, 75, 100)

DCM_TAXA = (1000,)
DCM_SUBTREE_SIZE = (100,)

SUP = "SUP"
SCS = "SCS"
MCS = "MCS"
BCD = "BCD"

SCRIPT_PATH = "scripts/method_scripts/"
RESULTS_FOLDER = "results/"

SCRIPTS = {
    SUP: "run_sup.sh",
    SCS: "run_scs.sh",
    MCS: "run_mcs.sh",
    BCD: "run_bcd.sh",
}

OPTIONS = {}
DEFAULT_OPTIONS = []


class ResultsLogger:
    def __init__(self, results_folder: str, experiment_directory: str) -> None:
        self.write_directory = results_folder + experiment_directory
        self.file_suffix = "_results.tsv"

        if not os.path.exists(self.write_directory):
            os.makedirs(self.write_directory)

    def result_already_exists(self, method: str, source_tree_file: str) -> bool:
        file_path = self.format_file_path(method)
        if not os.path.exists(file_path):
            return False

        with open(file_path, "r") as f:
            for line in f:
                mtf, stf, wall_time, tree = line.strip("\n").split("\t")
                if stf == source_tree_file:
                    return True
        return False

    def write_results(
        self,
        method: str,
        model_tree_file: Optional[str],
        source_tree_file: str,
        wall_time: float,
        tree: TreeNode,
    ) -> None:
        with open(self.format_file_path(method), "a") as f:
            f.write(
                str(model_tree_file)
                + "\t"
                + source_tree_file
                + "\t"
                + str(wall_time)
                + "\t"
                + str(tree)
                + "\n"
            )

    def format_file_path(self, method: str) -> str:
        return self.write_directory + method + self.file_suffix


def run_methods(
    source_tree_file: str,
    model_tree_file: str,
    methods: List[str],
    verbosity: int = 1,
    force_bifurcating: bool = False,
    calculate_distances: bool = True,
    logger: Optional[ResultsLogger] = None,
):
    results = {}
    for method in methods:
        if logger is not None:
            if logger.result_already_exists(method, source_tree_file):
                if verbosity >= 1:
                    print(
                        "Result already exists for",
                        method,
                        "on",
                        source_tree_file + "... skipping.",
                    )
                continue

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
            tree = make_tree(result.stdout.decode("utf-8").strip())
            if force_bifurcating:
                tree = tree.bifurcating()
        except:
            tree = None

        results[method] = (tree, end_time - start_time)

        if logger is not None:
            logger.write_results(
                method, model_tree_file, source_tree_file, end_time - start_time, tree  # type: ignore
            )

    with open(model_tree_file, "r") as f:
        model = make_tree(f.read().strip())
        if force_bifurcating:
            model = model.bifurcating()

    for method in methods:
        if method not in results:
            continue  # If it was skipped earlier
        tree, time_result = results[method]
        if verbosity >= 1:
            if verbosity >= 2:
                print("Computing distances for", method)
            if tree is None:
                print(f"{method}: time={time_result:.2f}s failed: {source_tree_file}) ")
            elif calculate_distances or verbosity >= 2:
                rf = rooted_rf_distance(model, tree)
                mc = matching_cluster_distance(model, tree)
                print(f"{method}: time={time_result:.2f}s RF={rf} MC={mc}")
            else:
                print(f"{method}: time={time_result:.2f}s")

    return results


def run_experiment_super_triplets(
    d: int,
    k: int,
    methods: List,
    verbosity: int = 1,
    calculate_distances: bool = True,
    result_logging: bool = False,
):
    assert d in SUPER_TRIPLET_D
    assert k in SUPER_TRIPLET_K

    if result_logging:
        logger = ResultsLogger(RESULTS_FOLDER, f"SuperTripletsBenchmark/d{d}/k{k}/")
    else:
        logger = None

    number_of_experiments = 100

    for i in range(number_of_experiments):
        tree_number = i + 1
        source_file = f"data/SuperTripletsBenchmark/source-trees/d{d}/k{k}/data-d{d}-k{k}-{tree_number}_phybp-s0r.nwk.source_trees"
        model_file = f"data/SuperTripletsBenchmark/model-trees/model-{tree_number}.nwk.model_tree"
        print(f"Results for {i} ({source_file}):")
        run_methods(
            source_file,
            model_file,
            methods,
            verbosity=verbosity,
            calculate_distances=calculate_distances,
            logger=logger,
        )


def run_experiment_smidgen(
    taxa: int,
    density: int,
    methods: List,
    verbosity: int = 1,
    calculate_distances: bool = True,
    result_logging: bool = False,
):
    assert taxa in SMIDGEN_TAXA
    assert density in SMIDGEN_DENSITY

    if result_logging:
        logger = ResultsLogger(RESULTS_FOLDER, f"superfine/{taxa}-taxa/{density}/")
    else:
        logger = None

    number_of_experiments = 10 if taxa == 1000 else 30

    for i in range(number_of_experiments):
        file = f"data/superfine/{taxa}-taxa/{density}/sm_data.{i}"
        print(f"Results for {i} ({file}):")
        run_methods(
            file + ".source_trees",
            file + ".model_tree",
            methods,
            verbosity=verbosity,
            calculate_distances=calculate_distances,
            logger=logger,
        )


def run_experiment_smidgen_og(
    taxa: int,
    density: int,
    methods: List,
    verbosity: int = 1,
    calculate_distances: bool = True,
    result_logging: bool = False,
):
    assert taxa in SMIDGEN_OG_NORMAL_TAXA or taxa == 10000
    if taxa == 10000:
        assert density == 0
    else:
        assert density in SMIDGEN_OG_DENSITY

    if result_logging:
        logger = ResultsLogger(RESULTS_FOLDER, f"SMIDGenOutgrouped/{taxa}/{density}/")
    else:
        logger = None

    number_of_experiments = 30 if taxa != 10000 else 10

    for i in range(number_of_experiments):
        source_file = f"data/SMIDGenOutgrouped/{taxa}/{density}/Source_Trees/RaxML/smo.{i}.sourceTreesOG.tre"
        model_file = f"data/SMIDGenOutgrouped/{taxa}/{density}/Model_Trees/pruned/smo.{i}.modelTree.tre"
        if verbosity >= 1:
            print(f"Results for {i} ({source_file}):")
        run_methods(
            source_file,
            model_file,
            methods,
            verbosity=verbosity,
            calculate_distances=calculate_distances,
            logger=logger,
        )


def run_experiment_dcm_exact(
    taxa: int,
    subtree_size: int,
    methods: List,
    verbosity: int = 1,
    calculate_distances: bool = True,
    result_logging: bool = False,
):
    assert taxa in DCM_TAXA
    assert subtree_size in DCM_SUBTREE_SIZE

    if result_logging:
        logger = ResultsLogger(
            RESULTS_FOLDER, f"birth_death/{taxa}/dcm_source_trees/{subtree_size}/"
        )
    else:
        logger = None

    number_of_experiments = 10

    for i in range(number_of_experiments):
        source_file = f"data/birth_death/{taxa}/dcm_source_trees/{subtree_size}/bd.{i}.source_trees"
        model_file = f"data/birth_death/{taxa}/model_trees/bd.{i}.model_tree"
        if verbosity >= 1:
            print(f"Results for {i} ({source_file}):")
        run_methods(
            source_file,
            model_file,
            methods,
            verbosity=verbosity,
            calculate_distances=calculate_distances,
            logger=logger,
        )


def run_experiment_dcm(
    taxa: int,
    subtree_size: int,
    methods: List,
    verbosity: int = 1,
    calculate_distances: bool = True,
    result_logging: bool = False,
):
    assert taxa in DCM_TAXA
    assert subtree_size in DCM_SUBTREE_SIZE

    if result_logging:
        logger = ResultsLogger(
            RESULTS_FOLDER, f"birth_death/{taxa}/iq_source_trees/{subtree_size}/"
        )
    else:
        logger = None

    number_of_experiments = 10

    for i in range(number_of_experiments):
        source_file = f"data/birth_death/{taxa}/iq_source_trees/{subtree_size}/bd.{i}.source_trees"
        model_file = f"data/birth_death/{taxa}/model_trees/bd.{i}.model_tree"
        if verbosity >= 1:
            print(f"Results for {i} ({source_file}):")
        run_methods(
            source_file,
            model_file,
            methods,
            verbosity=verbosity,
            calculate_distances=calculate_distances,
            logger=logger,
        )
