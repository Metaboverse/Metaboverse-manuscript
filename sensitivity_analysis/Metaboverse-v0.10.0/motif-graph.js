/*
*********************************
MODIFIED FROM METABOVERSE-V0.10.0
*********************************



Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

MIT License

Copyright (c) 2022 Metaboverse

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
var FileSaver = require('file-saver');

function exportTable(motifs, entity_species_reverse_dictionary, filename) {
    // for each item (condition) in motifs, output with column describing the condition
    var i = 0;
    var rows = [];
    let column_names = [
    "index", 
    "condition", 
    "id", 
    "name", 
    "collapsed", 
    "reactants", 
    "products", 
    "modifiers", 
    "additional_components",
    "magnitude_change",
    "source_p_value",
    "target_p_value",
    "averaged_p_value"
    ];
    rows.push(column_names);

    for (let condition in motifs) {
    for (let motif in motifs[condition]) {
        let this_entry = motifs[condition][motif];

        let these_reactants = [];
        for (let r in this_entry["reactants"]) {
        if (this_entry["reactants"][r] in entity_species_reverse_dictionary) {
            these_reactants.push(entity_species_reverse_dictionary[this_entry["reactants"][r]][0]);
        } else {
            these_reactants.push(this_entry["reactants"][r]);
        }
        }
        let these_products = [];
        for (let p in this_entry["products"]) {
        if (this_entry["products"][p] in entity_species_reverse_dictionary) {
            these_products.push(entity_species_reverse_dictionary[this_entry["products"][p]][0]);
        } else {
            these_products.push(this_entry["products"][p]);
        }
        }
        let these_modifiers = [];
        for (let m in this_entry["modifiers"]) {
        if (this_entry["modifiers"][m][0] in entity_species_reverse_dictionary) {
            these_modifiers.push([
            entity_species_reverse_dictionary[this_entry["modifiers"][m][0]][0],
            this_entry["modifiers"][m][1]
            ]);
        } else {
            these_modifiers.push(this_entry["modifiers"][m]);
        }
        }

        let entry = [
        String(i), 
        String(condition), 
        String(this_entry["id"]), 
        String(this_entry["name"]), 
        String(this_entry["collapsed"]), 
        String(these_reactants),
        String(these_products),
        String(these_modifiers), 
        String(this_entry["additional_components"]), 
        String(this_entry["magnitude_change"]), 
        String(this_entry["p_values"]["source"]),
        String(this_entry["p_values"]["target"]),
        String(this_entry["p_values"]["agg"])
        ];
        rows.push(entry);
        i += 1;
    }
    }

    // Source: https://stackoverflow.com/a/14966131/9571488
    let csvContent = "data:text/tab-separated-values;charset=utf-8,";
    rows.forEach(function(rowArray) {
      let row = rowArray.join("\t");
      csvContent += row + "\r\n";
    });

    var encodedUri = encodeURI(csvContent);
    // End code snippet

    // Source: https://stackoverflow.com/a/14966131/9571488
    var link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", filename);
    document.body.appendChild(link); // Required for FF
    link.click(); // This will download the data file named "my_data.csv".
    // End code snippet
}