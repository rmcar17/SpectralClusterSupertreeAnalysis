from cogent3.core.tree import PhyloNode


def check_ultrametric(tree: PhyloNode) -> bool:
    all_tips = tree.tips()
    ultrametric_length = all_tips[0].distance(tree)
    for i in range(1, len(all_tips)):
        if all_tips[i].distance(tree) != ultrametric_length:
            return False
    return True
