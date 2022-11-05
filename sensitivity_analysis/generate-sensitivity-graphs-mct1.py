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
mct1_output = os.path.join(__path__, "mct1_graphs")
if not os.path.exists(mct1_output):
    os.makedirs(mct1_output)
    
# Make yeast template files
if not os.path.exists(os.path.join(mct1_output, "SCE.mvdb")) \
or not os.path.exists(os.path.join(mct1_output, "SCE.nbdb")) \
or not os.path.exists(os.path.join(mct1_output, "SCE_template.mvrs")):
    print("Generating yeast template files...")
    yeast_template_cmd = """
        {0} curate \
            --output "{1}" \
            --organism_id "SCE" \
            --output_file "{2}" \
            --organism_curation_file "None" \
            --neighbor_dictionary_file "None" \
            --graph_template_file "None" \
            --transcriptomics "None" \
            --proteomics "None" \
            --metabolomics "None" \
            --database_source "reactome" \
            --experiment_type "None" \
            --experiment_name "yeast-templates" \
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
        mct1_output,
        os.path.join(mct1_output, "yeast-sensitivity.mvrs"),
        os.path.join(mct1_output, "progress_log.json"),
        os.path.join(mct1_output, "session_data.json")
    )
    os.system(yeast_template_cmd)
    
# Get list of transcriptomics datasets
mct1_input = os.path.join(__path__, "mct1_random_samples")
mct1_transcriptomics = next(os.walk(mct1_input), (None, None, []))[2]
print("\nGathered {0} files for mct1 transcriptomics.".format(len(mct1_transcriptomics)))

# Get list of proteomics datasets
mct1_proteomics = next(os.walk(mct1_input), (None, None, []))[2]
print("\nGathered {0} files for mct1 proteomics.".format(len(mct1_proteomics)))
    
# Get list of metabolomics datasets
mct1_metabolomics = next(os.walk(mct1_input), (None, None, []))[2]
print("\nGathered {0} files for mct1 metabolomics.".format(len(mct1_metabolomics)))

# Generate graph objects for yeast metabolomics
print("\nGenerating graph objects...")
for f in mct1_metabolomics:
    for i in range(len(mct1_proteomics)):
        yeast_graph_cmd = """
            {0} curate \
                --output "{1}" \
                --organism_id "HSA" \
                --output_file "{2}" \
                --organism_curation_file "{3}" \
                --neighbor_dictionary_file "{4}" \
                --graph_template_file "{5}" \
                --transcriptomics {6} \
                --proteomics {7} \
                --metabolomics "{8}" \
                --database_source "reactome" \
                --experiment_type "None" \
                --experiment_name "{9}" \
                --labels "0" \
                --blocklist "H+,H2O,CO2,e-,Mn2+,K+,Na+,Mg2+,O2,H2O2,NAD+,NADH,NADPH,NAAD,ATP,ADP,Pi,PPi,Ca2+" \
                --broadcast_genes \
                --broadcast_metabolites \
                --collapse_threshold "0.3" \
                --progress_log "{10}" \
                --session_data "{11}"
        """.format(
            metaboverse, 
            mct1_output,
            os.path.join(mct1_output, str(f.split(".")[0]) + ".mvrs"),
            os.path.join(mct1_output, "SCE.mvdb"),
            os.path.join(mct1_output, "SCE.nbdb"),
            os.path.join(mct1_output, "SCE_template.mvrs"),
            os.path.join(mct1_input, mct1_transcriptomics[i]),
            os.path.join(mct1_input, mct1_proteomics[i]),
            os.path.join(mct1_input, f),
            str(f.split(".")[0]) + "_" + str(mct1_proteomics[i].split(".")[0]),
            os.path.join(mct1_output, "progress_log.json"),
            os.path.join(mct1_output, "session_data.json")
        )
        os.system(yeast_graph_cmd)
    
    



