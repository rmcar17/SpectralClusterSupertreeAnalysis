import os
from typing import List, Optional

from ..distance.distance import matching_cluster_distance, rooted_rf_distance

from .experiment import BCD, MCS, RESULTS_FOLDER, SCS, SUP
from cogent3.core.tree import TreeNode
from cogent3 import make_tree

ORDERING = {BCD: 0, SCS: 1, SUP: 2, MCS: 3}


class DistanceLogger:
    def __init__(self, write_directory: str) -> None:
        self.write_directory = write_directory
        self.file_suffix = "_results_with_distances.tsv"

        if not os.path.exists(self.write_directory):
            os.makedirs(self.write_directory)

    def result_already_exists(self, method: str, source_tree_file: str) -> bool:
        file_path = self.format_file_path(method)
        if not os.path.exists(file_path):
            return False

        with open(file_path, "r") as f:
            for line in f:
                all_data = line.strip("\n").split("\t")
                mtf, stf, wall_time = all_data[:3]
                # b for forced bifurcation
                rf_distance, mc_distance, brf_distance, bmc_distance = all_data[3:7]
                tree = all_data[7]
                if stf == source_tree_file:
                    return True
        return False

    def write_results(
        self,
        method: str,
        model_tree_file: Optional[str],
        source_tree_file: str,
        wall_time: float,
        rf_distance: int,
        mc_distance: int,
        brf_distance: int,
        bmc_distance: int,
        tree: TreeNode,
    ) -> None:
        parts = [
            str(model_tree_file),
            source_tree_file,
            str(wall_time),
            str(rf_distance),
            str(mc_distance),
            str(brf_distance),
            str(bmc_distance),
            str(tree),
        ]
        with open(self.format_file_path(method), "a") as f:
            f.write("\t".join(parts) + "\n")

    def format_file_path(self, method: str) -> str:
        return self.write_directory + method + self.file_suffix


def calculate_distances_for_experiment(
    directory: str, result_files: List[str], verbosity: int = 1
):
    logger = DistanceLogger(directory + "/")

    result_files = sorted(
        result_files, key=lambda x: (ORDERING.get(x[:3], float("inf")), x)
    )

    methods = list(map(lambda x: x[:3], result_files))
    result_file_paths = list(map(lambda x: directory + "/" + x, result_files))

    file_objects = [
        open(result_file_path, "r") for result_file_path in result_file_paths
    ]

    next_lines = [file_object.readline() for file_object in file_objects]
    while any(next_lines):
        already_gave_stf = False
        for method, line in zip(methods, next_lines):
            if line == "":
                continue
            mtf, stf, wall_time, tree = line.strip("\n").split("\t")
            if logger.result_already_exists(method, stf):
                if verbosity >= 1:
                    print(
                        "Distance already exists for",
                        method,
                        "on",
                        stf + "... skipping.",
                    )
                continue

            if verbosity >= 1 and not already_gave_stf:
                print("Calculating distances for", stf)
                already_gave_stf = True

            all_data = line.strip("\n").split("\t")
            mtf, stf, wall_time, tree = all_data
            tree = make_tree(tree)

            with open(mtf, "r") as f:
                model_tree = make_tree(f.read().strip())

            rf = rooted_rf_distance(model_tree, tree)
            mc = matching_cluster_distance(model_tree, tree)

            b_model_tree = model_tree.bifurcating()
            b_tree = tree.bifurcating()

            brf = rooted_rf_distance(b_model_tree, b_tree)
            bmc = matching_cluster_distance(b_model_tree, b_tree)

            if verbosity >= 1:
                if brf != rf or bmc != mc:
                    print(f"{method}: RF={rf} MC={mc} BRF={brf} BMC={bmc}")
                else:
                    print(f"{method}: RF={rf} MC={mc}")
            logger.write_results(method, mtf, stf, wall_time, rf, mc, brf, bmc, tree)

        next_lines = [file_object.readline() for file_object in file_objects]

    if verbosity >= 1:
        print("Calculating Distances for")

    for file_object in file_objects:
        file_object.close()


def calculate_all_distances(verbosity=1):
    for root, subdirs, files in os.walk(RESULTS_FOLDER):
        result_files = list(filter(lambda x: x[3:] == "_results.tsv", files))
        if len(result_files) > 0:
            calculate_distances_for_experiment(root, result_files, verbosity=verbosity)
