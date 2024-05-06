# Spectral Cluster Supertree Analysis

[![DOI](https://zenodo.org/badge/705986177.svg)](https://zenodo.org/badge/latestdoi/705986177)

A repository containing experiments for the Spectral Cluster Supertree algorithm.

This repository was used to document experiments for the Spectral Cluster Supertree paper. For just the Spectral Cluster Supertree algorithm, please go to the linked [repository](https://github.com/rmcar17/SpectralClusterSupertree).

## Repository layout

- `data/` - directory containing miscellaneous data.
- `images/` - figures generated from results
- `methods/` - programs which run each supertree method evaluated
- `results/` - results pertaining to each method over the datasets
- `scripts/` - scripts for running each method
- `src/` - cli interface and programs used for data generation/distance calculation, etc.
- `tests/` - tests for the matching cluster distacnce metric
- `tmp/` - temporary directory

This repository also contains an implementation of the disk-covering method Rec-I-DCM3 (Roshan et al., 2004) for the case of rooted trees.
The function exists in `src/data_generation/dcm3.py` as `dcm3(guide_tree: PhyloNode, max_problem_size: int) -> List[PhyloNode]` which recursively
decomposes a `cogent3` tree object into a collection of overlapping subtrees according to the algorithm. Each subtree is of size at most `max_problem_size`,
except in the case the tree cannot be decomposed anymore.

## Installation

### Spectral Cluster Supertree (Required)

Follow the installation instructions the Spectral Cluster Supertree [repository](https://github.com/rmcar17/SpectralClusterSupertree)

### Command line interface (Required)

In the same python environment in which Spectral Cluster Supertree was installed, run ``pip install flit`` then ``flit install -s`` under this directory.

### Bad Clade Deletion (Optional)

The JAR file from [this website](https://bio.informatik.uni-jena.de/software/bcd/) comes prepackaged in this repository. It requires `java` to be available from the command line.

### IQTree2 (Optional)

Download IQ-TREE 2 from the [official website](http://www.iqtree.org/), and ensure the command line application is available from the path.

## Download Datasets (Optional)

A script has been set up to download and unpack all datasets from Zenodo (<https://zenodo.org/doi/10.5281/zenodo.11118021>).

To download, run `python get_data.py` in the terminal. The script may take some time to run as it downloads then unpacks ~1GB of compressed data.

## Command Line Interface

### Experiment Commands

#### Running an Experiment

`scsa run-experiment [OPTIONS] DATASET_NAME DATASET_PARAMS`

Runs a supertree experiment over the given methods on a specific dataset. See the help for more information.

#### Calculating Distance Metrics

`scsa calculate-distances [OPTIONS]`

Calculates the distance between the estimated and model trees for a given experiment. See the help for more information.

#### Plotting Graphs

`scsa plot`

Draws graphs for all experiments distances have been calculated for. See the help for more information.

### Data Generation

#### Creating Birth-Death Model Trees

`scsa create-bd-trees [OPTIONS] NUM_TREES NUM_TAXA`

Creates birth-death trees over the specified parameters. Depending on the input options, the trees may not be ultrametric. See the help for more information.

#### Simulating Sequence Alignment over Model Tree

`scsa sim-bd-seqs [OPTIONS] NUM_TAXA`

Given the generated model trees of a specific number of taxa, simulates
sequence alignments for the taxa over each of the trees under a
strand-symmetric general nucleotide model. See the help for more
information.

#### Creating DCM Source Trees

`scsa dcm-source-trees [OPTIONS] NUM_TAXA MAX_SUBPROBLEM_SIZE`

Given the generated model trees of a certain number of taxa, applies
the DCM3 algorithm to generate source trees of a maximal size.

The maximal size may be violated (usually in cases where it is too small). See the help for more information.

#### Generating IQTree Source Trees

`scsa iqtree [OPTIONS] NUM_TAXA MAX_SUBPROBLEM_SIZE`

 Given the taxa in each of the DCM source trees for a number of taxa and
maximum subproblem size, and the sequence alignments over all the taxa in
the model tree, generates source trees using IQTree2 under a strand-
symmetric model. See the help for more information.
