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

// find global motifs at beginning, save as list, check each graph for members
function gatherMotifs(data, categories) {

  var excl_hubs = true;
  var hub_threshold = 50;

  var update_output = update_nodes_links(
    data.nodes,
    data.links
  );
  data.nodes = update_output[0];
  data.links = update_output[1];

  var dict_output = create_dictionaries(data.nodes);
  var expression_dict = dict_output[0];
  var stats_dict = dict_output[1];
  var inferred_dict = dict_output[2];

  var link_neighbors = create_link_neighbors(
    data.nodes,
    data.links
  );

  data.blocklist = data.species_blocklist;
  data.blocklist = complete_blocklist(
    data.blocklist,
    data.metadata.blocklist,
    data.nodes
  )
  var stat_type = data.metadata.stat_type;

  let threshold = 1;
  
  let collapsed_motifs_Avg = motifSearch_Avg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_ModReg = modifierReg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_Avg = motifSearch_Avg(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_ModReg = modifierReg(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  return [collapsed_motifs_Avg, collapsed_motifs_ModReg, motifs_Avg, motifs_ModReg];
}
