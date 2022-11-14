#######################################################
# Run scripts for sensitivity analysis of Metaboverse #
#######################################################

# Generate and activate environment
echo "Loading environment..."
conda env create -f sensitivity-requirements.yml 
conda activate sensitivity-analysis

# Generate sample datasets
echo "Generating datasets for sensitivity analysis..."
python generate-sensitivity-datasets.py

# Build sensitivity analysis graphs
echo "Generating graphs from lung dataset..."
python generate-sensitivity-graphs-lung.py
echo "Generating graphs from mct1 dataset..."
python generate-sensitivity-graphs-mct1.py

# Perform "average" and "modreg" reaction pattern search in python and output results
### Do with or without collapsing
### This step requires node.js
##### On linux, this can be installed using `sudo apt install nodejs`
echo "Extracting reaction patterns..."
npm start

# Run network metrics of collapsed vs non-collapsed graphs
### Small world, connectivity, etc
echo "Generating network metrics..."
python generate-network-metrics.py
### Compare results - top reactions
### Upset plots of what is missed
echo "Generating comparison analysis metrics..."
python generate-results-metrics.py