"""
Source:

https://lukejharmon.github.io/pcm/chapter10_birthdeath/#section-10.2-the-birth-death-model
"""

import random

from typing import Optional

from cogent3 import make_tree
from cogent3.core.tree import PhyloNode


def birth_death_tree(
    birth_rate: float,
    death_rate: float,
    stopping_time: Optional[float] = None,
    stopping_taxa: Optional[int] = None,
    restart_on_fail=False,
    rename=False,
    rng: Optional[random.Random] = None,
) -> PhyloNode:
    """
    Generates a rooted phylogenetic tree according to the birth-death model given a birth rate,
    death rate and stopping time (length of tips of generated tree from root) or stopping taxa number.

    If the birth-death model causes an extinction over the entire phylogeny, behaviors is
    dependent on restart_on_fail. If restart_on_fail is False, then raises a RuntimeException.
    Otherwise, attempts to form a birth-death tree again.

    Args:
        birth_rate (float): The birth rate
        death_rate (float): The death rate
        stopping_time (Optional[float]): Time for simulation
        stopping_taxa (Optional[int]): Total taxa to generate for
        restart_on_fail (bool, optional): If True, total extinction restarts the process. Otherwise throws a RuntimeException. Defaults to False.
        rename (bool, optional): If True, renames the leaves of the tree with formatted numbers prefixed with "t" at the end of the process. Defaults to False.

    Raises:
        RuntimeError: If the phylogeny entirely becomes extinct through the process. Can instead restart when restart_on_fail is True.

    Returns:
        PhyloNode: The resulting birth-death tree.
    """
    if stopping_time is None and stopping_taxa is None:
        raise ValueError(
            "Must choose a termination condition (time or number of taxa)."
        )
    if stopping_time is not None and stopping_time < 0:
        raise ValueError("Stopping time cannot be negative.")
    if stopping_taxa is not None and stopping_taxa < 2:
        raise ValueError("Stopping taxa must be at least 2.")

    if rng is None:
        rng = random.Random()

    total_rate = birth_rate + death_rate
    birth_probability = birth_rate / total_rate

    tree = make_tree("(t0:0.0,t1:0.0);")
    tips = list(tree.traverse(self_before=False, self_after=False))

    current_time = 0
    taxa_used = 2
    while True:
        waiting_time = rng.expovariate(len(tips) * total_rate)
        current_time += waiting_time

        if stopping_time is not None and current_time >= stopping_time:
            adder = stopping_time - (current_time - waiting_time)
            for tip in tips:
                tip.length += adder
            break

        for tip in tips:
            tip.length += waiting_time

        if len(tips) == stopping_taxa:
            break

        if rng.random() < birth_probability:
            # Handle birth event
            birth_tip = tips.pop(rng.randrange(0, len(tips)))
            old_tip_name = birth_tip.name
            assert isinstance(old_tip_name, str)
            birth_tip.name = None
            birth_tip.append(make_tree(old_tip_name + ":0.0;"))
            birth_tip.append(make_tree("t" + str(taxa_used) + ":0.0;"))
            taxa_used += 1
            tips.extend(birth_tip.children)
        else:
            # Handle death event
            extinct_tip = tips.pop(rng.randrange(0, len(tips)))

            parent: PhyloNode = extinct_tip.parent
            parent.remove_node(extinct_tip)

            assert len(parent.children) == 1
            grandparent: PhyloNode = parent.parent
            if grandparent is None:
                if restart_on_fail:
                    tree = make_tree("(t0:0.0,t1:0.0);")
                    tips = list(tree.traverse(self_before=False, self_after=False))
                    current_time = 0
                    taxa_used = 2
                    continue
                else:
                    raise RuntimeError("Extinction" + str(birth_probability))
            other_child = parent.children[0]
            other_child.length += parent.length
            if grandparent is None:
                tree = other_child
            else:
                grandparent.remove_node(parent)
                grandparent.append(other_child)

    if rename:
        digits_needed = len(str(len(tips) - 1))

        for i, tip in enumerate(tips):
            tip.name = f"t{i:0{digits_needed}}"

    return tree
