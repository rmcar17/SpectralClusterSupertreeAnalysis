import os

to_remove = []
to_rename_model_tree = []
to_rename_source_trees = []

for root, subfolders, subfiles in os.walk("SuperTripletsBenchmark/"):
    for file in subfiles:
        if file.endswith(".nex") or file.endswith(".txt"):
            to_remove.append(root + "/" + file)
        elif file.endswith(".nwk"):
            if "model-trees" in root:
                to_rename_model_tree.append(root + "/" + file)
            elif "source-trees" in root:
                to_rename_source_trees.append(root + "/" + file)

for file in to_remove:
    try:
        os.remove(file)
    except:
        pass  # It was previously removed

for file in to_rename_model_tree:
    os.rename(file, file + ".model_tree")

for file in to_rename_source_trees:
    os.rename(file, file + ".source_trees")
