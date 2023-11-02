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
