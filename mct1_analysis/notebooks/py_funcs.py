import os
import sys
import gzip
import json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import xpressplot as xp
import networkx as nx
from networkx.readwrite import json_graph

sns.set(font='arial')
sns.set(font_scale=float(.1))
jakes_cmap = sns.diverging_palette(212, 61, s=99, l=77, sep=1, n=16, center='dark') #Custom aesthetics


def read_table(
        url,
        sep='\t',
        index_col=0,
        header='infer',
        low_memory=False,
        compression='infer',
        comment='#'):
    """Read tab-delimited table
    """

    data = pd.read_csv(
        url,
        sep=sep,
        index_col=index_col,
        header=header,
        low_memory=low_memory,
        compression=compression,
        comment=comment)

    return data


def parse_genelist(
        url,
        gene_name,
        gene_id,
        sep='\t',
        comment='#'):
    """Parse gene list from file
    """

    data = pd.read_csv(
        url,
        sep=sep,
        comment=comment)

    gene_names = list(set(data[gene_name].tolist()))
    gene_ids = list(set(data[gene_id].tolist()))

    gene_dict = {}
    for index, row in data.iterrows():
        gene_dict[row[gene_id]] = row[gene_name]

    return gene_names, gene_ids, gene_dict


def g_decompress(
        path,
        file,
        output):

    g_open = open(os.path.join(path, output), "wb")
    with gzip.open(os.path.join(path, file), "rb") as f:
        g_data = f.read()
    g_open.write(g_data)
    g_open.close()


def make_gene_dict(__path__):

    g_decompress(
        path=os.path.join(
            __path__,
            "..",
            "data",
            "analysis_lists"),
        file="Saccharomyces_cerevisiae.R64-1-1.103.gtf.gz",
        output="Saccharomyces_cerevisiae.R64-1-1.103.gtf"
    )

    gtf = pd.read_csv(
        str(os.path.join(
            __path__,
            "..",
            "data",
            "analysis_lists",
            "Saccharomyces_cerevisiae.R64-1-1.103.gtf")),
        sep='\t',
        comment='#',
        low_memory=False,
        header=None)

    # Parse out old and new names from GTF
    gtf_genes = gtf.loc[gtf[2] == 'gene']
    gtf_genes['original'] = gtf[8].str.split("gene_id \"").str[1].str.split("\"; ").str[0]
    gtf_genes['new'] = gtf[8].str.split("gene_name \"").str[1].str.split("\"; ").str[0]
    gene_dict = pd.Series(gtf_genes['new'].values,index=gtf_genes['original']).to_dict()
    for k in gene_dict.keys():
        try:
            if gene_dict[k] is np.nan:
                gene_dict[k] = k
        except:
            gene_dict[k] = k

    os.remove(os.path.join(
        __path__,
        "..",
        "data",
        "analysis_lists",
        "Saccharomyces_cerevisiae.R64-1-1.103.gtf")
    )

    return gene_dict


def get_quantiles(results, _range):

    r = results.stack().quantile(_range)
    return r.iloc[0], r.iloc[1]


def read_corr_data(file="yeast_mct1_rnaseq_spqn_cor.tsv"):

    results = pd.read_csv(
        os.path.join(
            __path__,
            "..",
            "data",
            file),
        sep='\t',
        index_col=0,
        low_memory=False
    )

    return results


def get_graph(results, min, max, gene_dict=None):

    g = nx.Graph()

    for x in results.index.tolist():
        g.add_node(x)
        g.nodes()[x]['id'] = x
        if gene_dict != None:
            if x in gene_dict:
                g.nodes()[x]['name'] = gene_dict[x]
            else:
                g.nodes()[x]['name'] = x
        else:
            g.nodes()[x]['name'] = x
    print("Nodes:", str(len(g.nodes())))
    sys.stdout.flush()

    for x in results.index.tolist():
        for y in results.index.tolist():
            if results.at[x, y] >= max or results.at[x, y] <= min:
                g.add_edge(x, y)
                g.edges()[(x, y)]['r'] = results.at[x, y]
    print("Edges:", str(len(g.edges())))
    sys.stdout.flush()

    return g


def output_graph(g, file="mct1_rnaseq_correlations_spqn_graph.json"):

    data = json_graph.node_link_data(g)
    with open(os.path.join(
            __path__,
            "..",
            "data",
            file),
        'w') as f:
        json.dump(data, f, indent=4)  # Parse out as array for javascript
    print('Outputting graph to', os.path.join(
        __path__,
        "..",
        "data",
        file))

def read_graph():

    with open(os.path.join(
            __path__,
            "..",
            "data",
            "mct1_rnaseq_correlations_spqn_graph.json")) as f:
        js_graph = json.load(f)
    g = json_graph.node_link_graph(js_graph)

    return g


def make_subgraph(g, s, t, gene_dict=None):

    gg = nx.Graph()
    for x in nx.all_shortest_paths(g, source=s, target=t):
        for y in range(len(x)):
            try:
                if gene_dict != None:
                    gg.add_node(gene_dict[x[y]])
                    gg.add_edge(gene_dict[x[y]], gene_dict[x[y+1]])
                else:
                    gg.add_node(x[y])
                    gg.add_edge(x[y], x[y+1])
            except:
                pass

    data = json_graph.node_link_data(gg)
    
    if gene_dict != None:
        outfile = os.path.join(
            __path__,
            "..",
            "data",
            "mct1_rnaseq_correlations_spqn_graph_" + str(gene_dict[s]) + "_" + str(gene_dict[t]) + ".json")
    else:
        outfile = os.path.join(
            __path__,
            "..",
            "data",
            "mct1_rnaseq_correlations_spqn_graph_" + str(s) + "_" + str(t) + ".json")
    with open(
        outfile,
        'w') as f:
        json.dump(data, f, indent=4)
    print("Output as", outfile)
    
    return gg


def make_hist_gene(results, gene=''):

    f = results[gene].plot.hist(bins=50)
    fig = f.get_figure()
    fig.savefig(
        os.path.join(
            __path__,
            "..",
            "plots",
            "mct1_rnaseq_correlations_" + str(gene) + "_spqn_hist.png"
        )
    )


def make_clustermap(results):
    # Plot clustered heatmap of r values
    # generate graph, maybe more strict, make sure it was selective in the code above.

    ax = sns.clustermap(results,
        center=0,
        metric='euclidean',
        method='centroid',
        xticklabels=True,
        yticklabels=True,
        linewidths=0,
        linecolor='#DCDCDC',
        cmap=jakes_cmap,
        col_cluster=True,
        row_cluster=True,
        figsize=(100,100))
    plt.savefig(os.path.join(
            __path__,
            "..",
            "data",
            "mct1_rnaseq_correlations_clustered.pdf"), bbox_inches='tight')
