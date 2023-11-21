import os
from typing import List
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from scs_analysis.experiment.distance_calculator import ORDERING

from scs_analysis.experiment.experiment import SCS_FAST, BCDG, BCDN, MCS


METHOD_MAP = {SCS_FAST: "SCS", BCDG: "BCD (GSCM)", BCDN: "BCD (No GSCM)", MCS: "MCS"}


def graph_experiment(image_folder: str, root: str, distance_files: List[str]):
    distance_files = sorted(
        distance_files, key=lambda x: (ORDERING.get(x[:-12], float("inf")), x)
    )
    methods = [file[:-27] for file in distance_files]
    distance_file_paths = [root + "/" + file for file in distance_files]

    header = [
        "Model Tree",
        "Source Tree",
        "Wall Time",
        "CPU Time",
        "RF Distance",
        "Matching Cluster Distance",
        "F1 Score",
        "Bifurcating RF Distance",
        "Bifurcating Matching Cluster Distance",
        "Bifurcating F1 Score",
        "Tree",
    ]
    drops = ["Model Tree", "Source Tree", "Tree"]
    dfs = []
    for method, distance_file in zip(methods, distance_file_paths):
        if method not in METHOD_MAP:
            continue

        df = pd.read_csv(distance_file, delimiter="\t", names=header)
        df = df.drop(columns=drops)
        df["method"] = METHOD_MAP[method]
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    image_directory = image_folder + "/".join(root.split("/")[1:]) + "/"

    if not os.path.exists(image_directory):
        os.makedirs(image_directory)

    for col in set(header).difference(drops):
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(df, col)
            sns.boxplot(scaled, x="method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            sns.boxplot(df, x="method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(root)
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_all.png")
        plt.close()

        mcs_mask = df["method"] == METHOD_MAP[MCS]

        mcs_data = df[mcs_mask]
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(mcs_data, col)
            sns.boxplot(scaled, x="method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            sns.boxplot(mcs_data, x="method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(root)
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_mcs.png")
        plt.close()

        no_mcs_data = df[~mcs_mask & (df["method"] != METHOD_MAP[BCDN])]
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(no_mcs_data, col)
            sns.boxplot(scaled, x="method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            sns.boxplot(no_mcs_data, x="method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(root)
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_others.png")
        plt.close()


def make_time_scale_reasonable(df, column):
    df = df.copy()
    if df[column].max() < 600:
        return df, "seconds"
    df[column] /= 60
    if df[column].max() < 600:
        return df, "minutes"
    df[column] /= 60
    return df, "hours"


def graph_results(image_folder: str, results_folder: str):
    for root, subdirs, files in os.walk(results_folder):
        distance_files = list(
            filter(lambda x: x.endswith("_with_distances.tsv"), files)
        )
        print(root)
        if len(distance_files) > 0:
            print("GRAPHING", root)
            graph_experiment(image_folder, root, distance_files)
