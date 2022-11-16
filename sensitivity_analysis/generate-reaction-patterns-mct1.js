/*
Perform "average" and "modreg" reaction pattern search in python and output results
    Do with and without collapsing
    This step requires node.js:
        - On Linux, this can be installed by using `sudo apt install nodejs`
    This step requires npm:
        - On Linux, this can be installed by using ``
*/
//import .js scripts
// Scripts were sourced from Metaboverse v0.10.0 and untested lines of code were removed from the original scripts
const fs = require("fs");
const $ = require("jQuery");
const stat_value = 0.1; 

// Read in list of all files to analyze
console.log("Generating analysis list...")
let mct1_manifest = "mct1-graph-manifest.txt";
let mct1_file_list = fs.readFileSync(mct1_manifest).toString().split("\n");
mct1_file_list = mct1_file_list.filter(a => a !== "");
mct1_file_list = mct1_file_list.map(function(a){return a.replace(/(\r\n|\n|\r)/gm, "")});
console.log("Processing " + mct1_file_list.length + " .mvrs files from mct1 study...\n")

// Process mct1 files
let mct1_url = "mct1_graphs";
for (let f in mct1_file_list) {
    console.log("> Running " + mct1_file_list[f] + "...")
    // read in Metaboverse file
    var data = JSON.parse(fs.readFileSync(mct1_url + "/" + mct1_file_list[f]).toString());
    
    // Parse out network items
    let categories = data["categories"];
    let entity_species_reverse_dictionary = make_entity_species_r_dictionary(data);

    // Parse reaction patterns
    console.log(">>> Gathering patterns...")
    let motif_outputs = gatherMotifs(data, categories);
    let collapsed_motifs_Avg = motif_outputs[0];
    let collapsed_motifs_ModReg = motif_outputs[1];
    let motifs_Avg = motif_outputs[2];
    let motifs_ModReg = motif_outputs[3];

    // Export Average and ModReg patterns to tables for each network
    console.log(">>> Writing collapsed patterns to tables...")
    let filename_motifs_collapsed_avg = mct1_file_list[f].slice(0, -5) + "_avg_collapse.csv";
    let filename_motifs_collapsed_modreg = mct1_file_list[f].slice(0, -5) + "_mod_collapse.csv";
    exportTable(collapsed_motifs_Avg, entity_species_reverse_dictionary, filename_motifs_collapsed_avg);
    exportTable(collapsed_motifs_ModReg, entity_species_reverse_dictionary, filename_motifs_collapsed_modreg);
    
    console.log(">>> Writing non-collapsed patterns to tables...")
    let filename_motifs_avg = mct1_file_list[f].slice(0, -5) + "_avg_noCollapse.csv";
    let filename_motifs_modreg = mct1_file_list[f].slice(0, -5) + "_mod_noCollapse.csv";
    exportTable(motifs_Avg, entity_species_reverse_dictionary, filename_motifs_avg);
    exportTable(motifs_ModReg, entity_species_reverse_dictionary, filename_motifs_modreg);
}
    
alert("Processing complete. You may exit out of this window now.")