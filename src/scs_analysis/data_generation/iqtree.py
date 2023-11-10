import os
from typing import Set

from cogent3 import PhyloNode, Alignment, make_tree, load_aligned_seqs

from scs_analysis.data_generation.generate_model_trees import (
    BIRTH_DEATH_FOLDER,
)
import subprocess


TMP_IQ_TREE_PATH = "tmp/iqtree/"
TMP_IQ_TREE_ALN_FILE = "tmp.fasta"


def clean_tmp_iqtree_directory() -> None:
    if not os.path.exists(TMP_IQ_TREE_PATH):
        return
    files = sorted(os.listdir(TMP_IQ_TREE_PATH))
    for file in files:
        os.remove(TMP_IQ_TREE_PATH + file)


def generate_iq_tree_for(source_tree: PhyloNode, alignment: Alignment) -> PhyloNode:
    sub_aln: Alignment = alignment.take_seqs(source_tree.get_tip_names())  # type: ignore

    clean_tmp_iqtree_directory()
    if not os.path.exists(TMP_IQ_TREE_PATH):
        os.makedirs(TMP_IQ_TREE_PATH)
    sub_aln.write("tmp/iqtree/tmp.fasta")

    subprocess.run(
        [
            "iqtree2",
            "-s",
            TMP_IQ_TREE_PATH + TMP_IQ_TREE_ALN_FILE,
            "-m",
            "STRSYM",
            "--quiet",
        ]
    )
    with open(TMP_IQ_TREE_PATH + TMP_IQ_TREE_ALN_FILE + ".treefile") as f:
        iqtree = make_tree(f.read().strip())
        iqtree.length = None
    return iqtree


def generate_iq_trees(taxa: int, max_subproblem_size: int, verbosity=1):
    dcm_path = BIRTH_DEATH_FOLDER + f"{taxa}/dcm_source_trees/{max_subproblem_size}/"
    if not os.path.exists(dcm_path):
        raise IOError(
            f"Path {dcm_path} does not exist. Do you mean to generate the dcm trees for {taxa} taxa and {max_subproblem_size} max subproblem size first?"
        )

    aln_path = BIRTH_DEATH_FOLDER + f"{taxa}/sequences/"
    if not os.path.exists(aln_path):
        raise IOError(
            f"Path {aln_path} does not exist. Do you mean to generate the sequences for {taxa} taxa first?"
        )

    iq_path = BIRTH_DEATH_FOLDER + f"{taxa}/iq_source_trees/{max_subproblem_size}/"
    if not os.path.exists(iq_path):
        os.makedirs(iq_path)

    tree_files = sorted(
        list(
            filter(
                lambda x: x.startswith("bd.") and x.endswith(".source_trees"),
                os.listdir(dcm_path),
            )
        )
    )
    for i, file_name in enumerate(tree_files):
        if verbosity >= 1:
            print(f"Generating IQTrees for source trees {i+1} of {len(tree_files)}")

        tree_identifier = file_name.split(".")[1]

        alignment = load_aligned_seqs(aln_path + f"bd.{tree_identifier}.fasta")

        dcm_trees = []
        with open(dcm_path + file_name, "r") as f:
            for tree_line in f:
                dcm_trees.append(make_tree(tree_line.strip()))

        iq_trees = []
        for i, tree in enumerate(dcm_trees):
            if verbosity >= 1:
                print(f"Running IQTree on source tree {i+1} of {len(dcm_trees)}")
            iq_trees.append(generate_iq_tree_for(tree, alignment))

        with open(iq_path + f"bd.{tree_identifier}.source_trees", "w") as f:
            for tree in iq_trees:
                f.write(str(tree) + "\n")
