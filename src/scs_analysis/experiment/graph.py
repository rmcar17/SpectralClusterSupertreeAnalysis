import os
from typing import List
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from scs_analysis.experiment.distance_calculator import ORDERING


def graph_experiment(image_folder: str, root: str, distance_files: List[str]):
    distance_files = sorted(
        distance_files, key=lambda x: (ORDERING.get(x[:-12], float("inf")), x)
    )
    methods = [file[:-27] for file in distance_files]
    distance_file_paths = [root + "/" + file for file in distance_files]

    header = [
        "model",
        "source",
        "walltime",
        "cputime",
        "rf",
        "mc",
        "f1",
        "brf",
        "bmc",
        "bf1",
        "tree",
    ]
    dfs = []
    for method, distance_file in zip(methods, distance_file_paths):
        df = pd.read_csv(distance_file, delimiter="\t", names=header)
        df = df.drop(columns=["model", "source", "tree"])
        df["method"] = method
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    image_directory = image_folder + "/".join(root.split("/")[1:]) + "/"

    if not os.path.exists(image_directory):
        os.makedirs(image_directory)

    for col in set(header).difference(["model", "source", "tree"]):
        plt.figure()
        sns.violinplot(df, x="method", y=col)
        plt.title(root)
        plt.savefig(image_directory + f"{col}.png")
        plt.close()


def graph_results(image_folder: str, results_folder: str):
    for root, subdirs, files in os.walk(results_folder):
        distance_files = list(
            filter(lambda x: x.endswith("_with_distances.tsv"), files)
        )
        if len(distance_files) > 0:
            graph_experiment(image_folder, root, distance_files)
