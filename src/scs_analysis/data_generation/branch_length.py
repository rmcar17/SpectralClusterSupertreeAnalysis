from cogent3.core.tree import PhyloNode
import math


def check_ultrametric(tree: PhyloNode) -> bool:
    """
    Given a phylogenetic tree, check whether it is ultrametric.

    Args:
        tree (PhyloNode): A phylogenetic tree

    Returns:
        bool: Whether the tree is ultrametric
    """
    all_tips = tree.tips()
    ultrametric_length = all_tips[0].distance(tree)
    for i in range(1, len(all_tips)):
        if not math.isclose(all_tips[i].distance(tree), ultrametric_length):  # type: ignore
            return False
    return True


def set_tree_height(tree: PhyloNode, target_height: float) -> None:
    """
    Given an ultrametric phylogenetic tree, scales the branch lengths
    by a uniform constant so that the root-to-tip distances are at the
    target height.

    Args:
        tree (PhyloNode): A phylogenetic tree to scale the branch lengths of.
        target_height (float): The distance from the root to all tips
    """
    if not check_ultrametric(tree):
        raise ValueError("Tree is not ultrametric.")

    tip = next(tree.postorder(include_self=False))
    assert tip.is_tip()

    distance_to_root: float = tip.distance(tree)  # type: ignore

    scale_factor = target_height / distance_to_root

    for edge in tree.get_edge_vector(include_root=False):
        edge.length *= scale_factor
