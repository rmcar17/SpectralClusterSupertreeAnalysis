from .birth_death import birth_death_tree
from .tree_height import randomly_scale_tree_height, set_tree_height
import os

BIRTH_DEATH_FOLDER = "data/birth_death/"


def generate_model_tree(
    file_path: str,
    taxa: int,
    birth_rate: float = 1.0,
    death_rate: float = 0.2,
    target_height: float = 1.0,
    initial_scaling_factor: float = 1.0,
    adder_mean: float = 0.0,
    adder_std: float = 0.05,
    scaling_factor_lb: float = 0.05,
    scaling_factor_ub: float = 8.0,
) -> None:
    model_tree = birth_death_tree(
        birth_rate, death_rate, stopping_taxa=taxa, restart_on_fail=True, rename=False
    )
    set_tree_height(model_tree, target_height)
    randomly_scale_tree_height(
        model_tree,
        initial_scaling_factor,
        adder_mean,
        adder_std,
        scaling_factor_lb,
        scaling_factor_ub,
    )

    with open(file_path, "w") as f:
        f.write(str(model_tree))


def generate_model_trees(
    taxa: int,
    number_of_trees: int,
    birth_rate: float = 1.0,
    death_rate: float = 0.2,
    target_height: float = 1.0,
    initial_scaling_factor: float = 1.0,
    adder_mean: float = 0.0,
    adder_std: float = 0.05,
    scaling_factor_lb: float = 0.05,
    scaling_factor_ub: float = 8.0,
    verbosity=0,
) -> None:
    write_path = BIRTH_DEATH_FOLDER + f"{taxa}/model_trees/"
    if not os.path.exists(write_path):
        os.makedirs(write_path)
    for i in range(number_of_trees):
        if verbosity >= 2:
            print(f"Generating tree {i+1} of {number_of_trees}")
        model_tree_path = write_path + f"bd.{i}.model_tree"
        generate_model_tree(
            model_tree_path,
            taxa,
            birth_rate,
            death_rate,
            target_height,
            initial_scaling_factor,
            adder_mean,
            adder_std,
            scaling_factor_lb,
            scaling_factor_ub,
        )
