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
var fs = require("fs");
var $ = require("jQuery");
var stat_value = 0.1; 

// Read in list of all files to analyze
console.log("Generating analysis list...")
let lung_manifest = "lung-graph-manifest.txt";
let lung_file_list = fs.readFileSync(lung_manifest).toString().split("\n");
lung_file_list = lung_file_list.filter(a => a !== "");
lung_file_list = lung_file_list.map(function(a){return a.replace(/(\r\n|\n|\r)/gm, "")});
console.log("Processing " + lung_file_list.length + " .mvrs files from lung study...\n")

// Process lung files
let lung_url = "lung_graphs";
for (let f in lung_file_list) {
    console.log("> Running " + lung_file_list[f] + "...")
    // read in Metaboverse file
    var data = JSON.parse(fs.readFileSync(lung_url + "/" + lung_file_list[f]).toString());
    
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
    let filename_motifs_collapsed_avg = lung_file_list[f].slice(0, -5) + "_avg_collapse.tsv";
    let filename_motifs_collapsed_modreg = lung_file_list[f].slice(0, -5) + "_mod_collapse.tsv";
    exportTable(collapsed_motifs_Avg, entity_species_reverse_dictionary, filename_motifs_collapsed_avg);
    exportTable(collapsed_motifs_ModReg, entity_species_reverse_dictionary, filename_motifs_collapsed_modreg);
    
    console.log(">>> Writing non-collapsed patterns to tables...")
    let filename_motifs_avg = lung_file_list[f].slice(0, -5) + "_avg_noCollapse.tsv";
    let filename_motifs_modreg = lung_file_list[f].slice(0, -5) + "_mod_noCollapse.tsv";
    exportTable(motifs_Avg, entity_species_reverse_dictionary, filename_motifs_avg);
    exportTable(motifs_ModReg, entity_species_reverse_dictionary, filename_motifs_modreg);
    
    alert("Pause")
}
    
alert("Processing complete. You may exit out of this window now.")