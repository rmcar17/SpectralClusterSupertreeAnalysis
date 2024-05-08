import random
import click
from scs_analysis.data_generation.dcm_partition import dcm_precise_source_trees
from scs_analysis.data_generation.generate_model_trees import generate_model_trees
from scs_analysis.data_generation.iqtree import generate_iq_trees
from scs_analysis.data_generation.simulate_alignments import simulate_alignments
from scs_analysis.distance.distance import *
from scs_analysis.experiment.experiment import (
    BCD,
    DCM_SUBTREE_SIZE,
    DCM_TAXA,
    SCS_FAST,
    SMIDGEN_OG_DENSITY,
    SMIDGEN_OG_NORMAL_TAXA,
    MCS,
    SUPER_TRIPLET_D,
    SUPER_TRIPLET_K,
    run_experiment_dcm,
    run_experiment_dcm_exact,
    run_experiment_smidgen_og,
    run_experiment_super_triplets,
)
from scs_analysis.experiment.distance_calculator import (
    calculate_all_distances,
    calculate_experiment_distances,
)
from scs_analysis.experiment.graph import graph_results

import time

__author__ = "Robert McArthur"
__copyright__ = "Copyright 2023, Robert McArthur"
__version__ = "2024.5.7"


@click.group()
@click.version_option(__version__)  # add version option
def main():
    """
    A command line interface for the tools used to run experiments
    for the Spectral Cluster Supertree paper.
    """
    pass


# Common args
_verbose = click.option(
    "-v",
    "--verbose",
    default=0,
    show_default=True,
    help="verbosity level",
)
_seed = click.option(
    "-r",
    "--rand",
    default=time.time(),
    show_default=False,
    help="random seed; if not specified uses the system time.",
)


@main.command(no_args_is_help=True)
@click.option("-a", "--all", is_flag=True, help="include all supertree methods.")
@click.option("-b", "--bcd", is_flag=True, help="include bad clade deletion method.")
@click.option(
    "-s", "--scs", is_flag=True, help="include spectral cluster supertree method."
)
@click.option("-m", "--mcs", is_flag=True, help="include min-cut supertree method.")
@click.option("-n", "--name", type=str, help="name of method")
@click.argument("dataset-name", nargs=1, required=True, type=str)
@click.argument("dataset-params", nargs=2, required=True, type=(int, int))
@_verbose
@_seed
def run_experiment(
    all, bcd, scs, mcs, name, dataset_name, dataset_params, verbose, rand
):
    """
    Runs a supertree experiment over the given methods on a specific dataset.

    DATASET_NAME is the name of the dataset. One of [supertriplets, smidgen, smidgenog, smidgenog10000, dcmexact, dcm].

    DATASET_PARAMS is a tuple of integers referring to the first n numbers in that experiment's data folder.
    For example, `scs run-experiment smidgenog 100 20` would run the experiment on the 100 taxa 20 density dataset.
    If 0 is specified as one of the parameters, all numbers are used.

    """
    methods = []
    if all:
        methods = [BCD, SCS_FAST, MCS]
    else:
        if bcd:
            methods.append(BCD)
        if scs:
            methods.append(SCS_FAST)
        if mcs:
            methods.append(MCS)
    if name:
        methods = [name.upper()]

    dataset_name = dataset_name.lower()
    dataset_params = list(dataset_params)

    if dataset_name == "supertriplets":
        experiment = run_experiment_super_triplets
        if dataset_params[0] == 0:
            dataset_params[0] = SUPER_TRIPLET_D
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = SUPER_TRIPLET_K
        else:
            dataset_params[1] = [dataset_params[1]]
    elif dataset_name == "smidgenog":
        experiment = run_experiment_smidgen_og
        if dataset_params[0] == 0:
            dataset_params[0] = SMIDGEN_OG_NORMAL_TAXA
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = SMIDGEN_OG_DENSITY
        else:
            dataset_params[1] = [dataset_params[1]]
    elif dataset_name == "smidgenog10000":
        experiment = run_experiment_smidgen_og
        dataset_params[0] = [10000]
        dataset_params[1] = [0]
    elif dataset_name == "dcmexact":
        experiment = run_experiment_dcm_exact
        if dataset_params[0] == 0:
            dataset_params[0] = DCM_TAXA
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = DCM_SUBTREE_SIZE
        else:
            dataset_params[1] = [dataset_params[1]]
    elif dataset_name == "dcm":
        experiment = run_experiment_dcm
        if dataset_params[0] == 0:
            dataset_params[0] = DCM_TAXA
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = DCM_SUBTREE_SIZE
        else:
            dataset_params[1] = [dataset_params[1]]
    else:
        raise ValueError("Invalid Experiment")

    rng = random.Random(rand)

    for param_1 in dataset_params[0]:
        for param_2 in dataset_params[1]:
            experiment(
                param_1,
                param_2,
                methods,
                verbosity=verbose,
                calculate_distances=False,
                result_logging=True,
                rng=rng,
            )


EXPERIMENT_FOLDER_IDENTIFIERS = {
    "supertriplets": "SuperTripletsBenchmark",
    "smidgenog": "SMIDGenOutgrouped",
    "dcmexact": "dcm_source_trees",
    "dcm": "iq_source_trees",
}


@main.command(no_args_is_help=False)
@click.option(
    "-e",
    "--experiment",
    default="all",
    show_default=True,
    help="The experiment to calculate distances for. One of [all, supertriplets, smidgen, smidgenog, dcmexact, dcm].",
)
@_verbose
def calculate_distances(experiment, verbose):
    """
    Calculates the distance between the estimated and model trees for a given experiment.

    When no experiment is specified, runs on all experiments.
    """
    if experiment == "all":
        calculate_all_distances(verbosity=verbose)
    else:
        calculate_experiment_distances(
            EXPERIMENT_FOLDER_IDENTIFIERS[experiment], verbosity=verbose
        )


@main.command(no_args_is_help=False)
@_verbose
def plot(verbose):
    """
    Draws graphs for all experiments distances have been calculated for.
    """
    graph_results("images/", verbosity=verbose)


@main.command(no_args_is_help=True)
@click.option(
    "-b",
    "--birth_rate",
    default=1.0,
    show_default=True,
    help="the birth rate for the birth-death process",
)
@click.option(
    "-d",
    "--death_rate",
    default=0.2,
    show_default=True,
    help="the death rate for the birth-death process",
)
@click.option(
    "-t",
    "--target_height",
    default=1.0,
    show_default=True,
    help="the distance from the root to every tip in the generated birth-death tree -- if 0 no scaling is applied",
)
@click.option(
    "-i",
    "--initial_scaling_factor",
    default=1.0,
    show_default=True,
    help="the initial scaling factor to scale the branch length at the root",
)
@click.option(
    "-n",
    "--normal_dist_params",
    default=(0.0, 0.05),
    type=(float, float),
    show_default=True,
    help="the normal distribution that is added to the scaling factor when moving dow the tree",
)
@click.option(
    "-s",
    "--scaling_factor_bounds",
    type=(float, float),
    default=(0.05, 8.0),
    show_default=True,
    help="the minimum and maximum bounds of the scaling factor",
)
@click.argument("num_trees", type=int)
@click.argument("num_taxa", type=int)
@_verbose
@_seed
def create_bd_trees(
    birth_rate,
    death_rate,
    target_height,
    initial_scaling_factor,
    normal_dist_params,
    scaling_factor_bounds,
    num_trees,
    num_taxa,
    verbose,
    rand,
):
    """
    Creates birth-death trees with the specified parameters.
    Depending on the input options, the trees may not be ultrametric.

    NUM_TREES is the number of trees it generates.

    NUM_TAXA is the number of taxa in the resulting tree.
    """
    rng = random.Random(rand)
    generate_model_trees(
        num_taxa,
        num_trees,
        birth_rate,
        death_rate,
        target_height,
        initial_scaling_factor,
        *normal_dist_params,
        *scaling_factor_bounds,
        verbosity=verbose,
        rng=rng,
    )


@main.command(no_args_is_help=True)
@click.option(
    "-l",
    "--length",
    default=1000,
    show_default=True,
    help="the length of the sequence alignment",
)
@click.argument("num_taxa", type=int)
@_verbose
def sim_bd_seqs(length, num_taxa, verbose):
    """
    Given the generated model trees of a specific number of taxa,
    simulates sequence alignments for the taxa over each of the trees.

    This is performed under a strand-symmetric general nucleotide model.
    """
    simulate_alignments(num_taxa, length, verbosity=verbose)


@main.command(no_args_is_help=True)
@click.argument("num_taxa", type=int)
@click.argument("max_subproblem_size", type=int)
@_verbose
def dcm_source_trees(num_taxa, max_subproblem_size, verbose):
    """
    Generates DCM source trees for the given number of taxa.
    Generates source trees down to a maximum subproblem size.

    NUM_TAXA is the number of taxa.

    MAX_SUBPROBLEM_SIZE is the maximum size of any generated tree.
    """
    dcm_precise_source_trees(num_taxa, max_subproblem_size, verbosity=verbose)


@main.command(no_args_is_help=True)
@click.argument("num_taxa", type=int)
@click.argument("max_subproblem_size", type=int)
@_verbose
def iqtree(num_taxa, max_subproblem_size, verbose):
    """
    Given the taxa in each of the DCM source trees for a number of
    taxa and maximum subproblem size, and the sequence alignments over
    all the taxa in the model tree, generates source trees using IQTree2
    under a strand-symmetric model.

    NUM_TAXA is the number of taxa.

    MAX_SUBPROBLEM_SIZE is the maximum size of the DCM trees.
    """
    generate_iq_trees(num_taxa, max_subproblem_size, verbosity=verbose)


if __name__ == "__main__":
    main()
