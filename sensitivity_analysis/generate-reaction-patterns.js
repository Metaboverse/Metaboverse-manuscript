/*
Perform "average" and "modreg" reaction pattern search in python and output results
    Do with or without collapsing
    This step requires node.js
        - On linux, this can be installed using `sudo apt install nodejs`
*/
//import .js scripts

// Read in list of all files to analyze
let file_number;
let file_list;
let lung_url = "";
let mct1_url = "";

//

console.log("Processing " + file_number + " .mvrs files...")

// Process files
for (let f in file_list) {
    
    // read in Metaboverse file
    //

    // Parse out network items
    let data, categories;
    //

    // Parse reaction patterns
    let motif_outputs = gatherMotifs(data, categories);
    let global_collapsed_motifs = motif_outputs[0];
    let global_motifs = motif_outputs[1];

    // Export Average and ModReg patterns to tables for each network
    
    

}
    

