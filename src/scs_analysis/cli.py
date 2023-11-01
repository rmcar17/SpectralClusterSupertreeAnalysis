from pathlib import Path

import click
from scs_analysis.distance.distance import *
from scs_analysis.experiment.experiment import (
    BCD,
    SCS,
    SUP,
    MCS,
    run_experiment_smidgen,
    run_experiment_smidgen_og,
    run_experiment_super_triplets,
)
from scs_analysis.experiment.distance_calculator import calculate_all_distances

from scitrack import CachingLogger


__author__ = "Robert McArthur"
__copyright__ = "Copyright 2023, Robert McArthur"
__credits__ = ["Robert McArthur"]
__license__ = "BSD"
__version__ = "2023.10.30"  # A DATE BASED VERSION
__maintainer__ = "Robert McArthur"
__email__ = "robert.mcarthur@anu.edu.au"
__status__ = "alpha"


LOGGER = CachingLogger()


@click.group()
@click.version_option(__version__)  # add version option
def main():
    """docstring explains what your cl app does"""
    pass


# note that I often define reusable options in a central place
_verbose = click.option(
    "-v",
    "--verbose",
    default=0,
    show_default=True,
    help="verbosity level",
)


# you can define custom parsers / validators
def _parse_csv_arg(*args) -> list:
    return args[-1].split(",")


_names = click.option(
    "--names",
    callback=_parse_csv_arg,
    help="converts comma separated values",
)

_outpath = click.option(
    "-o", "--outpath", type=Path, help="the input string will be cast to Path instance"
)


# the no_args_is_help=True means help is displayed if a
# user doesn't provide any arguments to a subcommand.
# Should be a click default I think!
@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "--achoice",
    type=click.Choice(["choice1", "choice2"]),
    default="choice1",
    help="make a choice",
)
@_names
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "--ensembl_account",
    envvar="ENSEMBL_ACCOUNT",
    help="shell variable with MySQL account "
    "details, e.g. export "
    "ENSEMBL_ACCOUNT='myhost.com jill jills_pass'",
)
@_verbose
def demo_log(infile, outpath, achoice, names, overwrite, verbose, ensembl_account):
    """what demo_log subcommand does"""
    # capture the local variables, at this point just provided arguments
    LOGGER.log_args()
    LOGGER.log_versions("numpy")
    LOGGER.input_file(infile)

    LOGGER.log_file_path = outpath / "some_path.log"


def _parse_csv_arg(*args):
    return args[-1].split(",")


@main.command(no_args_is_help=True)
@click.argument("message", required=True, type=str)
@click.option("-t", "--test", is_flag=True, help="test run")
@_verbose
def demo_echo(message, test, verbose):
    """what demo_echo subcommand does"""
    for _ in range(verbose):
        click.secho(message, fg="blue")


@main.command(no_args_is_help=True)
@click.option("-a", "--all", is_flag=True)
@click.option("-b", "--bcd", is_flag=True)
@click.option("-s", "--scs", is_flag=True)
@click.option("-u", "--sup", is_flag=True)
@click.option("-m", "--mcs", is_flag=True)
@click.argument("dataset", nargs=1, required=True, type=str)
@click.argument("dataset-params", nargs=2, required=True, type=(int, int))
@_verbose
def run_experiment(all, bcd, scs, sup, mcs, dataset, dataset_params, verbose):
    methods = []
    if all:
        methods = [BCD, SCS, SUP, MCS]
    else:
        if bcd:
            methods.append(BCD)
        if scs:
            methods.append(SCS)
        if sup:
            methods.append(SUP)
        if mcs:
            methods.append(MCS)

    dataset = dataset.lower()
    dataset_params = list(dataset_params)

    if dataset == "supertriplets":
        experiment = run_experiment_super_triplets
        if dataset_params[0] == 0:
            dataset_params[0] = [25, 50, 75]
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = [10, 20, 30, 40, 50]
        else:
            dataset_params[1] = [dataset_params[1]]
    elif dataset == "smidgen":
        experiment = run_experiment_smidgen
        if dataset_params[0] == 0:
            dataset_params[0] = [100, 500, 1000]
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = [20, 50, 75, 100]
        else:
            dataset_params[1] = [dataset_params[1]]
    elif dataset == "smidgenog":
        experiment = run_experiment_smidgen_og
        if dataset_params[0] == 0:
            dataset_params[0] = [100, 500, 1000]
        else:
            dataset_params[0] = [dataset_params[0]]
        if dataset_params[1] == 0:
            dataset_params[1] = [20, 50, 75, 100]
        else:
            dataset_params[1] = [dataset_params[1]]
    elif dataset == "smidgenog10000":
        experiment = run_experiment_smidgen_og
        if len(dataset) < 2:
            dataset.append([10000])
        if len(dataset) < 3:
            dataset.append([0])
    else:
        raise ValueError("Invalid Experiment")

    for param_1 in dataset_params[0]:
        for param_2 in dataset_params[1]:
            experiment(
                param_1,
                param_2,
                methods,
                verbosity=verbose,
                calculate_distances=False,
                result_logging=True,
            )


@main.command(no_args_is_help=True)
@_verbose
def calculate_distances(verbose):
    calculate_all_distances(verbosity=verbose)


if __name__ == "__main__":
    main()
