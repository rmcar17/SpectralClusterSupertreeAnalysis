# Spectral Cluster Supertree Analysis

A repository containing experiments for the Spectral Cluster Supertree algorithm.

This repository was used to document experiments for the Spectral Cluster Supertree paper. For just the Spectral Cluster Supertree algorithm, please go to the linked [repository](https://github.com/rmcar17/SpectralClusterSupertree).

## Installation

### Spectral Cluster Supertree (Required)

Follow the installation instructions the Spectral Cluster Supertree [repository](https://github.com/rmcar17/SpectralClusterSupertree)

### Command line interface (Required)

In the same python environment in which Spectral Cluster Supertree was installed, run ``flit install -s`` under this directory.

### Bad Clade Deletion (Optional)

The JAR file from [this website](https://bio.informatik.uni-jena.de/software/bcd/) comes prepackaged in this repository. It requires `java` to be available from the command line.

### IQTree2 (Optional)

Download IQ-TREE 2 from the [official website](http://www.iqtree.org/), and ensure the command line application is available from the path. 

## Download Datasets

TODO put link to Zenodo containing data files. Perhaps add script to unpack correctly.

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


#### Stats