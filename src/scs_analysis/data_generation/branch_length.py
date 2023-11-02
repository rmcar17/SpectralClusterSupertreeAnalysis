from cogent3.core.tree import PhyloNode
import math


def is_ultrametric(tree: PhyloNode) -> bool:
    """
    Given a phylogenetic tree, check whether it is ultrametric.

    Args:
        tree (PhyloNode): A phylogenetic tree.

    Returns:
        bool: Whether the tree is ultrametric.
    """
    all_tips = tree.tips()
    ultrametric_length = all_tips[0].distance(tree)
    for i in range(1, len(all_tips)):
        if not math.isclose(all_tips[i].distance(tree), ultrametric_length):  # type: ignore
            return False
    return True


def is_tree_correct_height(tree: PhyloNode, height: float) -> bool:
    """
    Given a tree, checks whether the distance from the root to each
    tip is the given height.

    Args:
        tree (PhyloNode): A phylogenetic tree.
        height (float): The expected root-to-tip distance.

    Returns:
        bool: True if the length from each tip to the root is the height, False otherwise.
    """
    for tip in tree.tips():
        if not math.isclose(tip.distance(tree), height):  # type: ignore
            return False
    return True


def set_tree_height(tree: PhyloNode, target_height: float) -> None:
    """
    Given an ultrametric phylogenetic tree, scales the branch lengths
    by a uniform constant so that the root-to-tip distances are at the
    target height.

    Args:
        tree (PhyloNode): A phylogenetic tree to scale the branch lengths of.
        target_height (float): The distance from the root to all tips.
    """
    if not is_ultrametric(tree):
        raise ValueError("Tree is not ultrametric.")

    tip = next(tree.postorder(include_self=False))
    assert tip.is_tip()

    distance_to_root: float = tip.distance(tree)  # type: ignore

    scale_factor = target_height / distance_to_root

    for edge in tree.get_edge_vector(include_root=False):
        edge.length *= scale_factor


if __name__ == "__main__":
    from scs_analysis.data_generation.birth_death import birth_death_tree

    target_height = 0.214782135

    tree = birth_death_tree(1, 0.2, stopping_taxa=10, restart_on_fail=True, rename=True)

    print(tree, is_ultrametric(tree))

    tip = next(tree.postorder(include_self=False))
    distance_to_root: float = tip.distance(tree)  # type: ignore

    print("Tree height is", distance_to_root, "Target height is", target_height)

    set_tree_height(tree, target_height)

    tip = next(tree.postorder(include_self=False))
    new_height: float = tip.distance(tree)  # type: ignore

    print(tree)

    print("New height is", new_height)

    print("All successfully changed", is_tree_correct_height(tree, target_height))
