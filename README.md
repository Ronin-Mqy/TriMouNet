# TriMouNet

TriMouNet (Trinet Multi-locus Network) reconstructs level-1 phylogenetic networks from multilocus gene-tree data by inferring best-fitting trinets (with support scores) and merging them into a network.

## Files and folders

TriMouNet.jar # main executable (network reconstruction from .tnets)
TriMouNet_User_Manual.doc # user manual
extract_outgroup_quartets.py # Step 1: quartet extraction + preprocessing
infer_trinet_types_and_generate_tnets.py # Step 2: trinet inference + .tnets generation

simulatednetwork_strong/ # simulated example (strong signal)
simulatednetwork_middle/ # simulated example (middle signal)
simulatednetwork_weak/ # simulated example (weak signal)
yeast/ # real-data example
bird/ # real-data example
Jun/ # real-data example

# Requirements

- Python 3.8+ (for scripts)
- Java 11+ (for `TriMouNet.jar`)
- Graphviz (optional, to render `.dot`)

## Quick start (example: simulatednetwork_middle)

Input: a multilocus gene-tree file (.tre) with one Newick tree per line
(all ingroup taxa + a fixed outgroup taxon O)
Step 1 — quartet extraction and preprocessing

python extract_outgroup_quartets.py
--input simulatednetwork_middle/simulatednetwork_middle.tre
--outdir simulatednetwork_middle/

Step 2 — infer trinets and generate a dense .tnets file

python infer_trinet_types_and_generate_tnets.py
--input simulatednetwork_middle/
--out simulatednetwork_middle/simulatednetwork_middle.tnets

Step 3 — reconstruct the network from .tnets

java -jar TriMouNet.jar
-o simulatednetwork_middle/myOutput.dot
simulatednetwork_middle/simulatednetwork_middle.tnets

Optional — visualize the network

dot -Tpng simulatednetwork_middle/myOutput.dot -o simulatednetwork_middle/myOutput.png
