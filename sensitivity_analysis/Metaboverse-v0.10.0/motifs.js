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
var path = require("path");
var jStat = require("jstat");

var eval_collapsed = false;
var eval_modifiers = false;
var excl_hubs = true;
var infer_complexes = true;
var hub_threshold = 50;
var showed_alert = false;
var showed_alert_alt = false;

function collapsedChecked() {
  if (eval_collapsed === false) {
    eval_collapsed = true;
  } else {
    eval_collapsed = false;
  }
  console.log("Motif evaluation includes collapsed reactions: ", eval_collapsed)
}

function modifiersChecked() {
  if (eval_modifiers === false) {
    eval_modifiers = true;
  } else {
    eval_modifiers = false;
  }
  console.log("Motif evaluation includes modifiers: ", eval_modifiers)
}

function hubsChecked() {
  if (excl_hubs === false) {
    excl_hubs = true;
  } else {
    excl_hubs = false;
  }
  console.log("High hub exlusion (more than", hub_threshold, "perturbations): ", excl_hubs)
}

function complexChecked() {
  if (infer_complexes === false) {
    infer_complexes = true;
  } else {
    infer_complexes = false;
  }
  console.log("Infer complex values and stats: ", infer_complexes)
}

function cleanHubs(
    components,
    degree_dict,
    blocklist,
    hub_threshold) {

  let filtered_hubs = components.filter(function(x) {
    if (degree_dict[x] >= hub_threshold) {
    } else if (blocklist.includes(x)) {
    } else {
      return x
    }
  })
  return filtered_hubs;
}

function eval_ci(
    stat_array, 
    stat_value) {
  // For confidence interval array, take selected confidence interval, return 0 if no overlap, 1 if overlap 
  if (stat_array === undefined || stat_array === null) {
    return null;
  } else {
    let intervals_map = Object.assign({}, ...stat_array.map((x) => ({[x[0]]: x[1]})));
    let x1 = intervals_map[stat_value][0][0];
    let x2 = intervals_map[stat_value][0][1];
    let y1 = intervals_map[stat_value][1][0];
    let y2 = intervals_map[stat_value][1][1];
    let overlap = Math.max(x1, y1) <= Math.min(x2, y2);
    if (overlap === false) {
      return 0.00;
    } else {
      return 1.00;
    }
  }
}

function parseComponents(
    reaction,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    degree_dict,
    blocklist,
    sample_index) {

  let reactants = reaction.reactants;
  let products = reaction.products;
  let modifiers = reaction.modifiers;

  let source_expression = [];
  let target_expression = [];

  let temp_reactants = $.extend(true, [], reactants);
  let temp_products = $.extend(true, [], products);
  
  // If true, eval catalysts as products (driving more output) and inhibitors as reactants (preventing more output)
  if (eval_modifiers === true) {
    for (x in modifiers) {
      if (modifiers[x][1] === "catalyst") {
        temp_products.push(modifiers[x][0])
      } else if (modifiers[x][1] === "inhibitor") {
        temp_reactants.push(modifiers[x][0])
      } else {}
    }
  }

  // If false, consider each complex component as an individual piece in the reaction pattern 
  if (infer_complexes === false) {
    let add_reactants = [];
    let remove_reactants = [];
    for (let x in temp_reactants) {
      if (temp_reactants[x] in inferred_dict && inferred_dict[temp_reactants[x]] === "true") {
        remove_reactants.push(temp_reactants[x]);
        for (let n in link_neighbors[temp_reactants[x]]) {
          add_reactants.push(link_neighbors[temp_reactants[x]][n]);
        }
      }
    }
    temp_reactants = temp_reactants.filter(n => !remove_reactants.includes(n)).concat(add_reactants);

    let add_products = [];
    let remove_products = [];
    for (let x in temp_products) {
      if (temp_products[x] in inferred_dict && inferred_dict[temp_products[x]] === "true") {
        remove_products.push(temp_products[x]);
        for (let n in link_neighbors[temp_products[x]]) {
          add_products.push(link_neighbors[temp_products[x]][n]);
        }
      }
    }
    temp_products = temp_products.filter(n => !remove_products.includes(n)).concat(add_products);
  }

  let clean_reactants = [];
  let clean_products = [];
  if (excl_hubs === true) {
    clean_reactants = cleanHubs(
      temp_reactants,
      degree_dict,
      blocklist,
      hub_threshold)
    clean_products = cleanHubs(
      temp_products,
      degree_dict,
      blocklist,
      hub_threshold)
  } else {
    clean_reactants = temp_reactants;
    clean_products = temp_products;
  }

  clean_reactants.forEach(l => {
    let reactant_expr = expression_dict[l][sample_index];
    let reactant_stat = stats_dict[l][sample_index];
    if (stat_type === "array") {
      reactant_stat = eval_ci(reactant_stat, stat_value);
    }
    if (reactant_expr !== null && reactant_stat !== null) {
      source_expression.push([
        parseFloat(reactant_expr),
        parseFloat(reactant_stat)
      ]);
    }
  })

  clean_products.forEach(l => {
    let product_expr = expression_dict[l][sample_index];
    let product_stat = stats_dict[l][sample_index];
    if (stat_type === "array") {
      product_stat = eval_ci(product_stat, stat_value);
    }
    if (product_expr !== null && product_stat !== null) {
      target_expression.push([
        parseFloat(product_expr),
        parseFloat(product_stat),
      ]);
    }
  })

  let updated_source = source_expression.filter(function(value) {
    return !Number.isNaN(value);
  });
  let updated_target = target_expression.filter(function(value) {
    return !Number.isNaN(value);
  });

  return [updated_source, updated_target]
}

function parseComponentsMod(
    reaction,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    degree_dict,
    blocklist,
    sample_index) {

  let core = reaction.reactants.concat(reaction.products);
  let modifiers = reaction.modifiers;
  let temp_core = $.extend(true, [], core);
  let temp_mods = $.extend(true, [], modifiers.map((i) => i[0]));

  // If false, consider each complex component as an individual piece in the reaction pattern 
  if (infer_complexes === false) {
    let add_core = [];
    let remove_core = [];
    for (let x in temp_core) {
      if (temp_core[x] in inferred_dict && inferred_dict[temp_core[x]] === "true") {
        remove_core.push(temp_core[x]);
        for (let n in link_neighbors[temp_core[x]]) {
          add_core.push(link_neighbors[temp_core[x]][n]);
        }
      }
    }
    temp_core = temp_core.filter(n => !remove_core.includes(n)).concat(add_core);

    let add_mods = [];
    let remove_mods = [];
    for (let x in temp_mods) {
      if (temp_mods[x] in inferred_dict && inferred_dict[temp_mods[x]] === "true") {
        remove_mods.push(temp_mods[x]);
        for (let n in link_neighbors[temp_mods[x]]) {
          add_mods.push(link_neighbors[temp_mods[x]][n]);
        }
      }
    }
    temp_mods = temp_mods.filter(n => !remove_mods.includes(n)).concat(add_mods);
  }

  let clean_core = [];
  let clean_modifiers = [];
  if (excl_hubs === true) {
    clean_core = cleanHubs(
      temp_core,
      degree_dict,
      blocklist,
      hub_threshold)
    clean_modifiers = cleanHubs(
      temp_mods,
      degree_dict,
      blocklist,
      hub_threshold)
  } else {
    clean_core = temp_core;
    clean_modifiers = temp_mods;
  }

  let core_expression = [];
  let mods_expression = [];
  clean_core.forEach(l => {
    let core_expr = expression_dict[l][sample_index];
    let core_stats = stats_dict[l][sample_index];
    if (stat_type === "array") {
      core_stats = eval_ci(core_stats, stat_value);
    }
    if (core_expr !== null && core_stats !== null) {
      core_expression.push([parseFloat(core_expr), parseFloat(core_stats)]);
    }
  })

  clean_modifiers.forEach(l => {
    let mod_expr = expression_dict[l][sample_index];
    let mod_stats = stats_dict[l][sample_index];
    if (stat_type === "array") {
      mod_stats = eval_ci(mod_stats, stat_value);
    }
    if (mod_expr !== null && mod_stats !== null) {
      mods_expression.push([parseFloat(mod_expr), parseFloat(mod_stats)]);
    }
  })

  let updated_core = core_expression.filter(function(value) {
    return !Number.isNaN(value);
  });
  let updated_mods = mods_expression.filter(function(value) {
    return !Number.isNaN(value);
  });

  return [updated_core, updated_mods]
}


function parseComponentsEnzymes(
    reaction,
    nodes,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    degree_dict,
    blocklist,
    sample_index) {

  let core = reaction.reactants.concat(reaction.products);
  let modifiers = reaction.modifiers;

  let core_expression = [];
  let mods_expression = [];

  let temp_core = $.extend(true, [], core);
  let temp_mods = $.extend(true, [], modifiers);

  // If false, consider each complex component as an individual piece in the reaction pattern 
  if (infer_complexes === false) {
    let add_core = [];
    let remove_core = [];
    for (let x in temp_core) {
      if (temp_core[x] in inferred_dict && inferred_dict[temp_core[x]] === "true") {
        remove_core.push(temp_core[x]);
        for (let n in link_neighbors[temp_core[x]]) {
          add_core.push(link_neighbors[temp_core[x]][n]);
        }
      }
    }
    temp_core = temp_core.filter(n => !remove_core.includes(n)).concat(add_core);

    let add_mods = [];
    let remove_mods = [];
    for (let x in temp_mods) {
      if (temp_mods[x] in inferred_dict && inferred_dict[temp_mods[x]] === "true") {
        remove_mods.push(temp_mods[x]);
        for (let n in link_neighbors[temp_mods[x]]) {
          add_mods.push(link_neighbors[temp_mods[x]][n]);
        }
      }
    }
    temp_mods = temp_mods.filter(n => !remove_mods.includes(n)).concat(add_mods);
  }

  let clean_core = [];
  let clean_modifiers = [];
  if (excl_hubs === true) {
    clean_core = cleanHubs(
      temp_core, degree_dict, blocklist, hub_threshold)
    temp_mods = temp_mods.map(x => x[0]);
    clean_modifiers = cleanHubs(
      temp_mods, degree_dict, blocklist, hub_threshold)
  } else {
    clean_core = temp_core;
    clean_modifiers = temp_mods.map(x => x[0]);
  }

  clean_core.forEach(l => {
    if (nodes[l].sub_type !== "metabolite_component"
    && nodes[l].sub_type !== "") {

      let core_expr = expression_dict[l][sample_index];
      let core_stats = stats_dict[l][sample_index];
      if (stat_type === "array") {
        core_stats = eval_ci(core_stats, stat_value);
      }
      if (core_expr !== null && core_stats !== null) {
        core_expression.push([parseFloat(core_expr), parseFloat(core_stats)]);
      }
    }
  })
  clean_modifiers.forEach(l => {
    if (nodes[l].sub_type !== "metabolite_component"
    && nodes[l].sub_type !== "") {
      let mod_expr = expression_dict[l][sample_index];
      let mod_stats = stats_dict[l][sample_index];
      if (stat_type === "array") {
        mod_stats = eval_ci(mod_stats, stat_value);
      }
      if (mod_expr !== null && mod_stats !== null) {
        mods_expression.push([parseFloat(mod_expr), parseFloat(mod_stats)]);
      }
    }
  })

  let updated_core = core_expression.filter(function(value) {
    return !Number.isNaN(value);
  });
  let updated_mods = mods_expression.filter(function(value) {
    return !Number.isNaN(value);
  });

  return [updated_core, updated_mods]
}

function parseComponentsMetabolites(
    reaction,
    nodes,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    degree_dict,
    blocklist,
    sample_index) {

  let core = reaction.reactants.concat(reaction.products);
  let modifiers = reaction.modifiers;

  let core_expression = [];
  let mods_expression = [];

  let temp_core = $.extend(true, [], core);
  let temp_mods = $.extend(true, [], modifiers);
  
  // If false, consider each complex component as an individual piece in the reaction pattern 
  if (infer_complexes === false) {
    let add_core = [];
    let remove_core = [];
    for (let x in temp_core) {
      if (temp_core[x] in inferred_dict && inferred_dict[temp_core[x]] === "true") {
        remove_core.push(temp_core[x]);
        for (let n in link_neighbors[temp_core[x]]) {
          add_core.push(link_neighbors[temp_core[x]][n]);
        }
      }
    }
    temp_core = temp_core.filter(n => !remove_core.includes(n)).concat(add_core);

    let add_mods = [];
    let remove_mods = [];
    for (let x in temp_mods) {
      if (temp_mods[x] in inferred_dict && inferred_dict[temp_mods[x]] === "true") {
        remove_mods.push(temp_mods[x]);
        for (let n in link_neighbors[temp_mods[x]]) {
          add_mods.push(link_neighbors[temp_mods[x]][n]);
        }
      }
    }
    temp_mods = temp_mods.filter(n => !remove_mods.includes(n)).concat(add_mods);
  }

  let clean_core = [];
  let clean_modifiers = [];
  if (excl_hubs === true) {
    clean_core = cleanHubs(
      temp_core, degree_dict, blocklist, hub_threshold)
    temp_mods = temp_mods.map(x => x[0]);
    clean_modifiers = cleanHubs(
      temp_mods, degree_dict, blocklist, hub_threshold)
  } else {
    clean_core = temp_core;
    clean_modifiers = temp_mods.map(x => x[0]);
  }

  clean_core.forEach(l => {
    if (nodes[l].sub_type === "metabolite_component") {
      let core_expr = expression_dict[l][sample_index];
      let core_stats = stats_dict[l][sample_index];
      if (stat_type === "array") {
        core_stats = eval_ci(core_stats, stat_value);
      }
      if (core_expr !== null && core_stats !== null) {
        core_expression.push([parseFloat(core_expr), parseFloat(core_stats)]);
      }
    }
  })
  clean_modifiers.forEach(l => {
    if (nodes[l].sub_type === "metabolite_component") {
      let mod_expr = expression_dict[l][sample_index];
      let mod_stats = stats_dict[l][sample_index];
      if (stat_type === "array") {
        mod_stats = eval_ci(mod_stats, stat_value);
      }
      if (mod_expr !== null && mod_stats !== null) {
        mods_expression.push([parseFloat(mod_expr), parseFloat(mod_stats)]);
      }
    }
  })

  let updated_core = core_expression.filter(function(value) {
    return !Number.isNaN(value);
  });
  let updated_mods = mods_expression.filter(function(value) {
    return !Number.isNaN(value);
  });

  return [updated_core, updated_mods]
}


function parseComponentsTrans(
    reaction,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    degree_dict,
    blocklist,
    sample_index) {

  let reactants = reaction.reactants;
  let products = reaction.products;
  let modifiers = reaction.modifiers;

  let source_expression = [];
  let target_expression = [];
  let modifier_expression = [];

  let temp_reactants = $.extend(true, [], reactants);
  let temp_products = $.extend(true, [], products);
  let temp_modifiers = $.extend(true, [], modifiers);

  // If false, consider each complex component as an individual piece in the reaction pattern 
  if (infer_complexes === false) {
    let add_reactants = [];
    let remove_reactants = [];
    for (let x in temp_reactants) {
      if (temp_reactants[x] in inferred_dict && inferred_dict[temp_reactants[x]] === "true") {
        remove_reactants.push(temp_reactants[x]);
        for (let n in link_neighbors[temp_reactants[x]]) {
          add_reactants.push(link_neighbors[temp_reactants[x]][n]);
        }
      }
    }
    temp_reactants = temp_reactants.filter(n => !remove_reactants.includes(n)).concat(add_reactants);

    let add_products = [];
    let remove_products = [];
    for (let x in temp_products) {
      if (temp_products[x] in inferred_dict && inferred_dict[temp_products[x]] === "true") {
        remove_products.push(temp_products[x]);
        for (let n in link_neighbors[temp_products[x]]) {
          add_products.push(link_neighbors[temp_products[x]][n]);
        }
      }
    }
    temp_products = temp_products.filter(n => !remove_products.includes(n)).concat(add_products);

    let add_modifiers = [];
    let remove_modifiers = [];
    for (let x in temp_modifiers) {
      if (temp_modifiers[x] in inferred_dict && inferred_dict[temp_modifiers[x]] === "true") {
        remove_modifiers.push(temp_modifiers[x]);
        for (let n in link_neighbors[temp_modifiers[x]]) {
          add_modifiers.push(link_neighbors[temp_modifiers[x]][n]);
        }
      }
    }
    temp_modifiers = temp_modifiers.filter(n => !remove_modifiers.includes(n)).concat(add_modifiers);
  }

  let clean_reactants = [];
  let clean_products = [];
  let clean_modifiers = [];
  if (excl_hubs === true) {
    clean_reactants = cleanHubs(
      temp_reactants,
      degree_dict,
      blocklist,
      hub_threshold)
    clean_products = cleanHubs(
      temp_products,
      degree_dict,
      blocklist,
      hub_threshold)
    parse_modifiers = temp_modifiers.map(x => x[0]);
    clean_modifiers = cleanHubs(
      parse_modifiers,
      degree_dict,
      blocklist,
      hub_threshold)
  } else {
    clean_reactants = temp_reactants;
    clean_products = temp_products;
    clean_modifiers = temp_modifiers.map(x => x[0]);
  }

  clean_reactants.forEach(l => {
    let reactant_expr = expression_dict[l][sample_index];
    let reactant_stats = stats_dict[l][sample_index];
    if (stat_type === "array") {
      reactant_stats = eval_ci(reactant_stats, stat_value);
    }
    if (reactant_expr !== null && reactant_stats !== null) {
      source_expression.push([parseFloat(reactant_expr), parseFloat(reactant_stats)]);
    }
  })

  clean_products.forEach(l => {
    let product_expr = expression_dict[l][sample_index];
    let product_stats = stats_dict[l][sample_index];
    if (stat_type === "array") {
      product_stats = eval_ci(product_stats, stat_value);
    }
    if (product_expr !== null && product_stats !== null) {
      target_expression.push([parseFloat(product_expr), parseFloat(product_stats)]);
    }
  })

  clean_modifiers.forEach(l => {
    let modifier_expr = expression_dict[l][sample_index];
    let modifier_stats = stats_dict[l][sample_index];
    if (stat_type === "array") {
      modifier_stats = eval_ci(modifier_stats, stat_value);
    }
    if (modifier_expr !== null && modifier_stats !== null) {
      modifier_expression.push([parseFloat(modifier_expr), parseFloat(modifier_stats)]);
    }
  })

  let updated_source = source_expression.filter(function(value) {
    return !Number.isNaN(value);
  });
  let updated_target = target_expression.filter(function(value) {
    return !Number.isNaN(value);
  });
  let updated_modifier = modifier_expression.filter(function(value) {
    return !Number.isNaN(value);
  });

  return [updated_source, updated_target, updated_modifier]
}


function computeAvg(arr) {
  let arr_sum = arr[0];
  let arr_len = arr.length;
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] !== null) {
      arr_sum += arr[i];
    } else {
      arr_len -= 1;
    }
  }
  let arr_avg = arr_sum / arr_len;
  return arr_avg;
}


// search for motif 1
//let threshold = d3.select("#avg_num").node().value;
function motifSearch_Avg(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)

      let updated_source = comps[0];
      let updated_target = comps[1];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);
        
      if (source_values.length > 0 && target_values.length > 0) {
        let source_avg = computeAvg(source_values);
        let target_avg = computeAvg(target_values);

        if (Math.abs(source_avg - target_avg) >= threshold) {
          let p_source = Math.min(...source_stats);
          let p_target = Math.min(...target_stats);
          let reaction_copy = $.extend(true, {}, reaction);
          reaction_copy.p_values = {
            "source": p_source,
            "target": p_target,
            'agg': aggregate_p_values(flatten(source_stats.concat(target_stats)))
          };
          reaction_copy.magnitude_change = Math.abs(source_avg - target_avg);
          sample_motifs.add(reaction_copy);
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}

// MaxMax
//let threshold = d3.select("#maxmax_num").node().value;
function motifSearch_MaxMax(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);
      if (source_values.length > 0 && target_values.length > 0) {
        let source_max = Math.max(...source_values);
        let target_max = Math.max(...target_values);

        if (Math.abs(source_max - target_max) >= threshold) {
          let source_index = source_values.indexOf(source_max);
          let target_index = target_values.indexOf(target_max);
          let p_source = source_stats[source_index];
          let p_target = target_stats[target_index];
          let reaction_copy = $.extend(true, {}, reaction);
          reaction_copy.p_values = {
            "source": p_source,
            "target": p_target,
            'agg': aggregate_p_values([p_source, p_target])
          };
          reaction_copy.magnitude_change = Math.abs(source_max - target_max);
          sample_motifs.add(reaction_copy);
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}

// MinMin
//let threshold = d3.select("#minmin_num").node().value;
function motifSearch_MinMin(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);

      if (source_values.length > 0 && target_values.length > 0) {
        let source_min = Math.min(...source_values);
        let target_min = Math.min(...target_values);

        if (Math.abs(source_min - target_min) >= threshold) {
          let source_index = source_values.indexOf(source_min);
          let target_index = target_values.indexOf(target_min);
          let p_source = source_stats[source_index];
          let p_target = target_stats[target_index];
          let reaction_copy = $.extend(true, {}, reaction);
          reaction_copy.p_values = {
            "source": p_source,
            "target": p_target,
            'agg': aggregate_p_values([p_source, p_target])
          };
          reaction_copy.magnitude_change = Math.abs(source_min - target_min);
          sample_motifs.add(reaction_copy);
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}


//MaxMin
//let threshold = d3.select("#maxmin_num").node().value;
function motifSearch_MaxMin(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);
    
      if (updated_source.length > 0 && updated_target.length > 0) {
        let source_max = Math.max(...source_values);
        let target_min = Math.min(...target_values);
        if (Math.abs(source_max - target_min) >= threshold) {
          let source_index = source_values.indexOf(source_max);
          let target_index = target_values.indexOf(target_min);
          let p_source = source_stats[source_index];
          let p_target = target_stats[target_index];
          let reaction_copy = $.extend(true, {}, reaction);
          reaction_copy.p_values = {
            "source": p_source,
            "target": p_target,
            'agg': aggregate_p_values([p_source, p_target])
          };
          reaction_copy.magnitude_change = Math.abs(source_max - target_min);
          sample_motifs.add(reaction_copy);
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}


//MinMax
//let threshold = d3.select("#minmax_num").node().value;
function motifSearch_MinMax(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);

      if (updated_source.length > 0 && updated_target.length > 0) {
        let source_min = Math.min(...source_values);
        let target_max = Math.max(...target_values);

        if (Math.abs(source_min - target_max) >= threshold) {
          let source_index = source_values.indexOf(source_min);
          let target_index = target_values.indexOf(target_max);
          let p_source = source_stats[source_index];
          let p_target = target_stats[target_index];
          let reaction_copy = $.extend(true, {}, reaction);
          reaction_copy.p_values = {
            "source": p_source,
            "target": p_target,
            'agg': aggregate_p_values([p_source, p_target])
          };
          reaction_copy.magnitude_change = Math.abs(source_min - target_max);
          sample_motifs.add(reaction_copy);
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}

//Sustained
//let threshold = d3.select("#sustained_num").node().value;
// Will not include sustained motif if the value on both sides exactly the same
function motifSearch_Sustained(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);

      if (source_values.length > 0 && target_values.length > 0) {

        // Sustained up-regulation
        let up_in = [];
        let up_out = [];
        let p_source_up;
        let p_target_up;
        let magnitude_change_up;
        for (i in source_values) {
          if (source_values[i] >= threshold) {
            up_in.push(source_values[i]);
          }
        }
        for (j in target_values) {
          if (target_values[j] >= threshold) {
            up_out.push(target_values[j]);
          }
        }
        if (up_in.length > 0 && up_out > 0) {
          let source_up = Math.max(...up_in)
          let target_up = Math.max(...up_out)
          let source_index = source_values.indexOf(source_up);
          let target_index = target_values.indexOf(target_up);
          p_source_up = source_stats[source_index];
          p_target_up = target_stats[target_index];
          if (source_up === target_up) {
            magnitude_change_up = undefined;
          } else {
            magnitude_change_up = Math.max(Math.abs(source_up), Math.abs(target_up));
          }
        }

        // Sustained down-regulation
        let down_in = [];
        let down_out = [];
        let p_source_down;
        let p_target_down;
        let magnitude_change_down;
        for (k in source_values) {
          if (source_values[k] <= -(threshold)) {
            down_in.push(source_values[k]);
          }
        }
        for (l in target_values) {
          if (target_values[l] <= -(threshold)) {
            down_out.push(target_values[l]);
          }
        }
        if (down_in.length > 0 && down_out > 0) {
          let source_down = Math.min(...down_in)
          let target_down = Math.min(...down_out)
          let source_index = source_values.indexOf(source_down);
          let target_index = target_values.indexOf(target_down);
          p_source_down = source_stats[source_index];
          p_target_down = target_stats[target_index];

          if (source_down === target_down) {
            magnitude_change_down = undefined;
          } else {
            magnitude_change_down = Math.max(Math.abs(source_down), Math.abs(target_down));
          }
        }
      
        if (((down_in.length > 0) && (down_out.length > 0)) ||
        ((up_in.length > 0) && (up_out.length > 0))) {
          if (magnitude_change_up && magnitude_change_down) {
            if (magnitude_change_up >= magnitude_change_down) {
              let reaction_copy = $.extend(true, {}, reaction);
              reaction_copy.p_values = {
                "source": p_source_up,
                "target": p_target_up,
                'agg': aggregate_p_values([p_source_up, p_target_up])
              };
              reaction_copy.magnitude_change = magnitude_change_up;
              sample_motifs.add(reaction_copy);
            } else {
              let reaction_copy = $.extend(true, {}, reaction);
              reaction_copy.p_values = {
                "source": p_source_down,
                "target": p_target_down,
                'agg': aggregate_p_values([p_source_down, p_target_down])
              };
              reaction_copy.magnitude_change = magnitude_change_down;
              sample_motifs.add(reaction_copy);
            }
          } else if (magnitude_change_up !== undefined) {
            let reaction_copy = $.extend(true, {}, reaction);
            reaction_copy.p_values = {
              "source": p_source_up,
              "target": p_target_up,
              'agg': aggregate_p_values([p_source_up, p_target_up])
            };
            reaction_copy.magnitude_change = magnitude_change_up;
            sample_motifs.add(reaction_copy);
          } else if (magnitude_change_down !== undefined) {
            let reaction_copy = $.extend(true, {}, reaction);
            reaction_copy.p_values = {
              "source": p_source_down,
              "target": p_target_down,
              'agg': aggregate_p_values([p_source_down, p_target_down])
            };
            reaction_copy.magnitude_change = magnitude_change_down;
            sample_motifs.add(reaction_copy);
          } 
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}

function modifierReg(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  // If the net change between at least one modifier and one core component of a
  // reaction is greater than or equal to the threshold, return the reaction
  // *** Checking the "include modifiers" button will have no effect here
  let discovered_motifs = [];
  for (_idx in sample_indices) {
    let sample_motifs = new Set();
    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponentsMod(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_core = comps[0];
      let updated_modifiers = comps[1];

      if (updated_core.length > 0 && updated_modifiers.length > 0) {

        // Check each core/mod combination for a diff that meets threshold
        let same_diff = get_same_diff(updated_core, updated_modifiers, threshold);
        let big_diff = get_big_diff(updated_core, updated_modifiers, threshold);
        if (same_diff[0] !== 0 || big_diff[0] !== 0) {
          let mag_change, p_source, p_target;
          if (same_diff[0] > big_diff[0] && same_diff[0] !== 0) {
            mag_change = same_diff[0];
            p_source = same_diff[1][0][1];
            p_target = same_diff[1][1][1];
          } else if (big_diff[0] !== 0) {
            mag_change = big_diff[0];
            p_source = big_diff[1][0][1];
            p_target = big_diff[1][1][1];
          }

          let reaction_copy = $.extend(true, {}, reaction);
          reaction_copy.p_values = {
            'source': p_source,
            'target': p_target,
            'agg': aggregate_p_values([p_source, p_target])
          };
          reaction_copy.magnitude_change = mag_change;
          sample_motifs.add(reaction_copy);
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}

function modifierTransport(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    link_neighbors,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices) {
  // Highlight if modifier changed where inputs and outputs are same (minus
  //    compartment) --> regulation of transport reaction
  // If a componenet on both sides meets threshold and a modifier seperately
  // meets threshold, return
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = new Set();

    for (let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponentsTrans(
        reaction,
        expression_dict,
        stats_dict,
        stat_type,
        stat_value,
        inferred_dict,
        link_neighbors,
        degree_dict,
        blocklist,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];
      let updated_modifier = comps[2];
      let source_values = updated_source.map((i) => i[0]);
      let target_values = updated_target.map((i) => i[0]);
      let modifier_values = updated_modifier.map((i) => i[0]);
      let source_stats = updated_source.map((i) => i[1]);
      let target_stats = updated_target.map((i) => i[1]);
      let modifier_stats = updated_modifier.map((i) => i[1]);

      if (source_values.length > 0 &&
        target_values.length > 0 &&
        modifier_values.length > 0) {

        let intersect = source_values.filter(x => target_values.includes(x));
        for (x in intersect) {
          for (y in modifier_values) {
            if ((Math.abs(intersect[x] - modifier_values[y]) >= threshold)
            || (Math.abs(intersect[x]) >= threshold
                && Math.abs(modifier_values[y]) >= threshold)
            ) {
              let p_source = source_stats[x];
              let p_target = modifier_stats[y];
              let p_vals = {
                "source": p_source,
                "target": p_target,
                'agg': aggregate_p_values(flatten(source_stats.concat(modifier_stats)))
              };
              let place_value;
              if (Math.abs(intersect[x]) >= threshold
              && Math.abs(modifier_values[y]) >= threshold) {
                place_value = Math.max(
                  Math.abs(intersect[x]),
                  Math.abs(modifier_values[y])
                );
              } else {
                place_value = 0;
              }
              let magnitude_chg = Math.max(
                Math.abs(intersect[x] - modifier_values[y]),
                place_value
              );

              let reaction_copy = $.extend(true, {}, reaction);
              if (reaction in sample_motifs) {
                if (magnitude_chg > reaction_copy.magnitude_change) {
                  reaction_copy.p_values = p_vals;
                  reaction_copy.magnitude_change = magnitude_chg;
                  sample_motifs.add(reaction_copy);
                }
              } else {
                reaction_copy.p_values = p_vals;
                reaction_copy.magnitude_change = magnitude_chg;
                sample_motifs.add(reaction_copy);
              }
            }
          }
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}

// 1. Biggest difference b/t reactions passes threshold
function get_same_diff(
    base_comp,
    neighbor_comp,
    threshold) {
  let one_max = 0;
  let one_pair = [];
  let two_max = 0;
  let two_pair = [];
  for (let i in base_comp) {
    if (Math.abs(base_comp[i][0]) > one_max) {
      one_max = Math.abs(base_comp[i][0]);
      one_pair = base_comp[i];
    }
  }
  for (let i in neighbor_comp) {
    if (Math.abs(neighbor_comp[i][0]) > two_max) {
      two_max = Math.abs(neighbor_comp[i][0]);
      two_pair = neighbor_comp[i];
    }
  }
  if (one_max >= threshold
  && two_max >= threshold) {
    let same_max = Math.max(one_max, two_max);
    let same_diff = [one_pair, two_pair];
    return [same_max, same_diff];
  } else {
    return [0, 0];
  }
}

// 2. One component from each reaction passes threshold
function get_big_diff(
    base_comp,
    neighbor_comp,
    threshold) {
  
  let diff_max = 0;
  let diff_pair = [];

  for (let i in base_comp) {
    for (let j in neighbor_comp) {
      if (Math.abs(base_comp[i][0] - neighbor_comp[j][0]) > diff_max) {
        diff_max = Math.abs(base_comp[i][0] - neighbor_comp[j][0]);
        diff_pair = [base_comp[i], neighbor_comp[j]]
      }
    }
  }
  if (diff_max >= threshold) {
    return [diff_max, diff_pair];
  } else {
    return [0, 0];
  }
}

function enzymeMotif(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    neighbors_dictionary,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices,
    nodes,
    link_neighbors) {

  // If two neighboring reactions have significantly shifting enzyme conc.
  // If both change above/below a threshold, or the difference between the
  // two passes the threshold
  //
  // Consider all non-metabolites from each neighboring reaction
  // If largest difference passes, consider as a motif
  // Or if both change in the same direction and pass the threshold, return
  // Considers modifiers by default
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    // neighbors_dictionary
    let sample_motifs = new Set();

    for (let neighbor in neighbors_dictionary) {
      // check each neighbor pair for components (non-metabolite)
      if (neighbor in collapsed_reaction_dict && neighbors_dictionary[neighbor].length > 0) {
        // Get evaluated reaction info
        let reaction = collapsed_reaction_dict[neighbor];
        let comps = parseComponentsEnzymes(
          reaction,
          nodes,
          expression_dict,
          stats_dict,
          stat_type,
          stat_value,
          inferred_dict,
          link_neighbors,
          degree_dict,
          blocklist,
          _idx)
        comps = comps[0].concat(comps[1]);

        // Get neighboring reaction(s) info
        for (let n_reaction in neighbors_dictionary[neighbor]) {
          let this_neighbor = neighbors_dictionary[neighbor][n_reaction];
          if (this_neighbor in collapsed_reaction_dict) {
            let neighbor_reaction = collapsed_reaction_dict[this_neighbor];
            let neighbor_comps = parseComponentsEnzymes(
              neighbor_reaction,
              nodes,
              expression_dict,
              stats_dict,
              stat_type,
              stat_value,
              inferred_dict,
              link_neighbors,
              degree_dict,
              blocklist,
              _idx)
            neighbor_comps = neighbor_comps[0].concat(neighbor_comps[1]);

            // Check for conditions
            let same_diff = get_same_diff(comps, neighbor_comps, threshold);
            let big_diff = get_big_diff(comps, neighbor_comps, threshold);
            if (same_diff[0] !== 0 || big_diff[0] !== 0) {
              let collapsed = "false";
              if (reaction.collapsed === "true"
              || neighbor_reaction.collapsed === "true") {
                collapsed = "true";
              }
              let mag_change, p_source, p_target;
              if (same_diff[0] > big_diff[0] && same_diff[0] !== 0) {
                mag_change = same_diff[0];
                p_source = same_diff[1][0][1];
                p_target = same_diff[1][1][1];
              } else if (big_diff[0] !== 0) {
                mag_change = big_diff[0];
                p_source = big_diff[1][0][1];
                p_target = big_diff[1][1][1];
              }
            sample_motifs.add({
              'id': reaction.id + "_" + neighbor_reaction.id,
              'rxn1': $.extend(true, {}, reaction),
              'rxn2': $.extend(true, {}, neighbor_reaction),
              'collapsed': collapsed,
              'magnitude_change': mag_change,
              'p_values': {
                'source': p_source,
                'target': p_target,
                'agg': aggregate_p_values([p_source, p_target])
                }
              })
            }
          }
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['rxn1']['id']].concat(path_mapper[sample_motifs[m]['rxn2']['id']])
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}


function activityMotif(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    stats_dict,
    stat_type,
    stat_value,
    inferred_dict,
    neighbors_dictionary,
    path_mapper,
    degree_dict,
    blocklist,
    sample_indices,
    nodes,
    link_neighbors) {

  // Consider all metabolites from each neighboring reaction
  // If largest difference passes, consider as a motif
  // Or if both change in the same direction and pass the threshold, return
  // Considers modifiers by default
  let discovered_motifs = [];

  for (_idx in sample_indices) {

    // neighbors_dictionary
    let sample_motifs = new Set();

    for (let neighbor in neighbors_dictionary) {

      // check each neighbor pair for components (non-metabolite)
      if (neighbor in collapsed_reaction_dict && neighbors_dictionary[neighbor].length > 0) {

        // Get evaluated reaction info
        let reaction = collapsed_reaction_dict[neighbor];
        let comps = parseComponentsMetabolites(
          reaction,
          nodes,
          expression_dict,
          stats_dict,
          stat_type,
          stat_value,
          inferred_dict,
          link_neighbors,
          degree_dict,
          blocklist,
          _idx)
        comps = comps[0].concat(comps[1]);

        // Get neighboring reaction(s) info
        for (let n_reaction in neighbors_dictionary[neighbor]) {
          let this_neighbor = neighbors_dictionary[neighbor][n_reaction];
          if (this_neighbor in collapsed_reaction_dict) {
            let neighbor_reaction = collapsed_reaction_dict[this_neighbor];
            let neighbor_comps = parseComponentsMetabolites(
              neighbor_reaction,
              nodes,
              expression_dict,
              stats_dict,
              stat_type,
              stat_value,
              inferred_dict,
              link_neighbors,
              degree_dict,
              blocklist,
              _idx)
            neighbor_comps = neighbor_comps[0].concat(neighbor_comps[1]);

            // Check for conditions
            let same_diff = get_same_diff(comps, neighbor_comps, threshold);
            let big_diff = get_big_diff(comps, neighbor_comps, threshold);
            if (same_diff[0] !== 0 || big_diff[0] !== 0) {
              let collapsed = "false";
              if (reaction.collapsed === "true"
              || neighbor_reaction.collapsed === "true") {
                collapsed = "true";
              }
              let mag_change, p_source, p_target;
              if (same_diff[0] > big_diff[0] && same_diff[0] !== 0) {
                mag_change = same_diff[0];
                p_source = same_diff[1][0][1];
                p_target = same_diff[1][1][1];
              } else if (big_diff[0] !== 0) {
                mag_change = big_diff[0];
                p_source = big_diff[1][0][1];
                p_target = big_diff[1][1][1];
              }
              sample_motifs.add({
                'id': reaction.id + "_" + neighbor_reaction.id,
                'rxn1': $.extend(true, {}, reaction),
                'rxn2': $.extend(true, {}, neighbor_reaction),
                'collapsed': collapsed,
                'magnitude_change': mag_change,
                'p_values': {
                  'source': p_source,
                  'target': p_target,
                  'agg': aggregate_p_values([p_source, p_target])
                }
              })
            }
          }
        }
      }
    }
    sample_motifs = [...sample_motifs];
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['rxn1']['id']].concat(path_mapper[sample_motifs[m]['rxn2']['id']])
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
}


function make_neighbors_dictionary(
    data,
    degree_dict) {

  let nodes = data.nodes;
  let links = data.links;
  let neighbors_dictionary = {};

  for (let l in links) {
    let one = links[l].source;
    let two = links[l].target;

    if (excl_hubs === true
    && (degree_dict[one] >= hub_threshold
      || degree_dict[two] >= hub_threshold
    )) {
      // Skip if one of the two connected components don't match hub threshold
    } else {
      if (one in neighbors_dictionary) {
        neighbors_dictionary[one].add(two);
      } else {
        neighbors_dictionary[one] = new Set();
        neighbors_dictionary[one].add(two);
      }
      if (two in neighbors_dictionary) {
        neighbors_dictionary[two].add(one);
      } else {
        neighbors_dictionary[two] = new Set();
        neighbors_dictionary[two].add(one);
      }
    }
  }
  let neighbors = Object.keys(neighbors_dictionary);

  var reaction_neighbors_dictionary = {};
  let reactions = data.reaction_dictionary;
  for (let neighbor in neighbors) {
    if (neighbors[neighbor] in reactions) {
      let components = [...neighbors_dictionary[neighbors[neighbor]]];
      let connected_reactions = new Set()
      for (let _c in components) {
        let sub_components = [...neighbors_dictionary[components[_c]]];
        for (let _c_ in sub_components) {
          if (sub_components[_c_] in reactions &&
              neighbors[neighbor] !== sub_components[_c_]) {
            connected_reactions.add(sub_components[_c_]);
          }
        }
      }
      reaction_neighbors_dictionary[neighbors[neighbor]] = [
        ...connected_reactions
      ];
    }
  }

  var collapsed_reaction_neighbors_dictionary = {};
  let collapsed_reactions = data.collapsed_reaction_dictionary
  for (let neighbor in neighbors) {
    if (neighbors[neighbor] in collapsed_reactions) {
      let components = [...neighbors_dictionary[neighbors[neighbor]]];
      let connected_collapsed_reactions = new Set()
      for (let _c in components) {
        let sub_components = [...neighbors_dictionary[components[_c]]];
        for (let _c_ in sub_components) {
          if (sub_components[_c_] in collapsed_reactions &&
              neighbors[neighbor] !== sub_components[_c_]) {
            connected_collapsed_reactions.add(sub_components[_c_]);
          }
        }
      }
      collapsed_reaction_neighbors_dictionary[neighbors[neighbor]] = [
        ...connected_collapsed_reactions
      ];
    }
  }

  return [
    reaction_neighbors_dictionary,
    collapsed_reaction_neighbors_dictionary
  ];
}


function aggregate_p_values(p_array) {

  let agg_p = (Math.E * jStat.geomean(p_array));

  if (agg_p > 1.0) {
    agg_p = 1.0
  }

  return agg_p;
}


// Source: https://stackoverflow.com/a/39261045; CC BY-SA 3.0
function flatten(arr) {

  return arr.reduce(function(a, b) {
    return a.concat(Array.isArray(b) ? flatten(b) : b);
  }, []);

}