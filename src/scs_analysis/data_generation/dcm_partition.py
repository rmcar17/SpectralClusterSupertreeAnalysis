import os

import cogent3

from .dcm3 import dcm3
from .generate_model_trees import BIRTH_DEATH_FOLDER


def dcm_precise_source_trees(taxa: int, max_subproblem_size: int, verbosity=0):
    tree_path = BIRTH_DEATH_FOLDER + f"{taxa}/model_trees/"
    if not os.path.exists(tree_path):
        raise IOError(
            f"Path {tree_path} does not exist. Do you mean to generate model trees for {taxa} taxa first?"
        )
    dcm_path = BIRTH_DEATH_FOLDER + f"{taxa}/dcm_source_trees/{max_subproblem_size}/"
    if not os.path.exists(dcm_path):
        os.makedirs(dcm_path)

    tree_files = sorted(
        list(
            filter(
                lambda x: x.startswith("bd.") and x.endswith(".model_tree"),
                os.listdir(tree_path),
            )
        )
    )
    for i, file_name in enumerate(tree_files):
        if verbosity >= 2:
            print(f"Applying DCM3 to tree {i+1} of {len(tree_files)}")

        model_tree = cogent3.load_tree(tree_path + file_name)
        dcm_source_trees = dcm3(model_tree, max_subproblem_size)

        tree_identifier = file_name.split(".")[1]
        with open(dcm_path + f"bd.{tree_identifier}.source_trees", "w") as f:
            for tree in dcm_source_trees:
                f.write(str(tree) + "\n")
