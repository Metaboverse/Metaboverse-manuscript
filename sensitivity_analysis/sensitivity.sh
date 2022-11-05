#######################################################
# Run scripts for sensitivity analysis of Metaboverse #
#######################################################

# Generate and activate environment
conda env create -f sensitivity-requirements.yml 
conda activate sensitivity-analysis

# Generate sample datasets
python generate-sensitivity-datasets.py

# Build sensitivity analysis graphs
python generate-sensitivity-graphs-lung.py
python generate-sensitivity-graphs-mct1.py

# Perform "average" and "modreg" reaction pattern search in python and output results
### Do with or without collapsing
### This step requires node.js
##### On linux, this can be installed using `sudo apt install nodejs`
node generate-reaction-patterns.js

# Run network metrics of collapsed vs non-collapsed graphs
### Small world, connectivity, etc
### Compare results - top reactions
### Upset plots of what is missed
