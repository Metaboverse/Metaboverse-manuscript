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
let lung_manifest = "lung-graph-manifest.txt";
let lung_file_list = fs.readFileSync(lung_manifest).toString().split("\n");
lung_file_list = lung_file_list.filter(a => a !== "");
lung_file_list = lung_file_list.map(function(a){return a.replace(/(\r\n|\n|\r)/gm, "")});
console.log("Processing " + lung_file_list.length + " .mvrs files from lung study...\n")




/* TESTING */
lung_file_list = [lung_file_list[0]];
/* TESTING */



// Process files
let lung_url = "lung_graphs";
for (let f in lung_file_list) {
    console.log("> Running " + lung_file_list[f] + "...")
    // read in Metaboverse file
    var data = JSON.parse(fs.readFileSync(lung_url + "/" + lung_file_list[f]).toString());
    
    // Parse out network items
    let categories = data["categories"];

    // Parse reaction patterns
    let motif_outputs = gatherMotifs(data, categories);
    let global_collapsed_motifs = motif_outputs[0];
    let global_motifs = motif_outputs[1];

    // Export Average and ModReg patterns to tables for each network
    exportTable(motifs, sample_idx, exclude_idx, filename)
    

}
    



let mct1_url = "mct1_graphs";
let mct1_manifest = "mct1-graph-manifest.txt";








alert("Processing complete. You may exit out of this window now.")