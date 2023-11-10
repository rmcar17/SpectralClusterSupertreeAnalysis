import cogent3
import os
from cogent3.core.tree import PhyloNode
from numpy.random import default_rng
from .generate_model_trees import BIRTH_DEATH_FOLDER

rng = default_rng(seed=None)  # setting seed will allow for reproducibility


def get_sim_params():
    """fits a ssGN model to the sample alignment and returns the fitted likelihood function"""
    aln: cogent3.Alignment = cogent3.load_aligned_seqs(
        "data/alignment/197113_332182_17210.nexus.gz", moltype="dna"  # type: ignore
    )
    aln = aln.omit_gap_pos(allowed_gap_frac=0)
    tree = cogent3.make_tree("(197113,(332182,17210))")
    model = cogent3.get_app(
        "model", "ssGN", optimise_motif_probs=True, tree=tree, show_progress=False
    )
    result = model(aln)
    return result.lf


def sim_alignment(tree: PhyloNode, align_length: int) -> cogent3.ArrayAlignment:
    lf = get_sim_params()
    mprobs = lf.get_motif_probs()
    rates = [p for p in lf.get_param_names() if ">" in p]
    mles = {p: lf.get_param_value(par_name=p) for p in rates}

    sm = cogent3.get_model("ssGN")
    lf = sm.make_likelihood_function(tree)
    lf.set_motif_probs(mprobs)
    for p, v in mles.items():
        lf.set_param_rule(par_name=p, value=v)
    return lf.simulate_alignment(sequence_length=align_length, random_series=rng)


def simulate_alignments(taxa: int, align_length: int = 1000, verbosity=0):
    tree_path = BIRTH_DEATH_FOLDER + f"{taxa}/model_trees/"
    if not os.path.exists(tree_path):
        raise IOError(
            f"Path {tree_path} does not exist. Do you mean to generate model trees for {taxa} taxa first?"
        )
    seq_path = BIRTH_DEATH_FOLDER + f"{taxa}/sequences/"
    if not os.path.exists(seq_path):
        os.makedirs(seq_path)

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
            print(f"Generating sequences {i+1} of {len(tree_files)}")
        model_tree = cogent3.load_tree(tree_path + file_name)
        alignment = sim_alignment(model_tree, align_length)

        tree_identifier = file_name.split(".")[1]
        alignment.write(seq_path + f"bd.{tree_identifier}.fasta")
