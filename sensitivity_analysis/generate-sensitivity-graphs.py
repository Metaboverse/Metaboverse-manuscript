"""_summary_
Perform sensitivity analysis of Metaboverse using two test datasets presented in the manuscript
Using metaboverse-cli v0.10.0
"""

import os 
import sys 

__path__ = os.getcwd()

# Load resources
metaboverse =  os.path.join(__path__, "metaboverse-cli-windows.exe")
if not os.path.exists(metaboverse):
    print("\nDownloading metaboverse-cli executable for Windows...")
    exe = "https://github.com/Metaboverse/metaboverse-cli/releases/download/v0.10.0/metaboverse-cli-windows.exe"
    os.system("curl -OL {0}".format(exe))

# Create output directory
lung_output = os.path.join(__path__, "lung_graphs")
if not os.path.exists(lung_output):
    os.makedirs(lung_output)
    
# Make human template files
if not os.path.exists(os.path.join(lung_output, "HSA.mvdb")) \
or not os.path.exists(os.path.join(lung_output, "HSA.nbdb")) \
or not os.path.exists(os.path.join(lung_output, "HSA_template.mvrs")):
    print("Generating human template files...")
    human_template_cmd = """
        {0} curate \
            --output "{1}" \
            --organism_id "HSA" \
            --output_file "{2}" \
            --organism_curation_file "None" \
            --neighbor_dictionary_file "None" \
            --graph_template_file "None" \
            --transcriptomics "None" \
            --proteomics "None" \
            --metabolomics "None" \
            --database_source "reactome" \
            --experiment_type "None" \
            --experiment_name "human-templates" \
            --labels "0" \
            --blocklist "H+,H2O,CO2,e-,Mn2+,K+,Na+,Mg2+,O2,H2O2,NAD+,NADH,NADPH,NAAD,ATP,ADP,Pi,PPi,Ca2+" \
            --force_new_curation \
            --broadcast_genes \
            --broadcast_metabolites \
            --collapse_threshold "0.3" \
            --progress_log "{3}" \
            --session_data "{4}"
    """.format(
        metaboverse, 
        lung_output,
        os.path.join(lung_output, "human-sensitivity.mvrs"),
        os.path.join(lung_output, "progress_log.json"),
        os.path.join(lung_output, "session_data.json")
    )
    os.system(human_template_cmd)
    
# Get list of metabolomics datasets
lung_input = os.path.join(__path__, "lung_random_samples")
lung_files = next(os.walk(lung_input), (None, None, []))[2]
print("\nGathered {0} files for LUAD metabolomics.".format(len(lung_files)))

# Generate graph objects for lung metabolomics
print("\nGenerating graph objects...")
for f in lung_files:
    human_graph_cmd = """
        {0} curate \
            --output "{1}" \
            --organism_id "HSA" \
            --output_file "{2}" \
            --organism_curation_file "{3}" \
            --neighbor_dictionary_file "{4}" \
            --graph_template_file "{5}" \
            --transcriptomics "None" \
            --proteomics "None" \
            --metabolomics "{6}" \
            --database_source "reactome" \
            --experiment_type "None" \
            --experiment_name "{7}" \
            --labels "0" \
            --blocklist "H+,H2O,CO2,e-,Mn2+,K+,Na+,Mg2+,O2,H2O2,NAD+,NADH,NADPH,NAAD,ATP,ADP,Pi,PPi,Ca2+" \
            --broadcast_genes \
            --broadcast_metabolites \
            --collapse_threshold "0.3" \
            --progress_log "{8}" \
            --session_data "{9}"
    """.format(
        metaboverse, 
        lung_output,
        os.path.join(lung_output, str(f.split(".")[0]) + ".mvrs"),
        os.path.join(lung_output, "HSA.mvdb"),
        os.path.join(lung_output, "HSA.nbdb"),
        os.path.join(lung_output, "HSA_template.mvrs"),
        os.path.join(lung_input, f),
        str(f.split(".")[0]),
        os.path.join(lung_output, "progress_log.json"),
        os.path.join(lung_output, "session_data.json")
    )
    os.system(human_graph_cmd)
    
    



