import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json
import pandas as pd
from Bio import Phylo
from datetime import datetime
import networkx as nx
from networkx.readwrite import json_graph


def _read_json(json_file):
    with open(json_file) as fin:
        return json.load(fin)

def _write_json(data, outfile):
    with open(outfile, "w") as fout:
        json.dump(data, fout)

def read_raw_graph(distance_json, meta=None):
    """Load nodes and edges from distances json file.

    :param distance_json: file path.  eg. {"000005825-000005845": 0.833528, "000005825-000006605": 0.870156, ...}
    """

    print("{} start to read raw graph.".format(datetime.now()))
    if meta:
        metadata = pd.read_csv(meta, sep='\t', dtype={'GCF': str}, index_col=None)
        GCFs = tuple(metadata['GCF'])
        print("Size of meta: {}".format(metadata.shape))

    fg = nx.Graph()
    distances = _read_json(distance_json)
    print("Length of distance: {}".format(len(distances)))

    # fg.add_weighted_edges_from((npair.split('-')[1], npair.split('-')[1], dis)
    #                            for npair, dis in distances.items())

    for npair, dis in distances.items():
        n1, n2 = npair.split('-')
        if meta and (not n1 in GCFs or not n2 in GCFs):
            continue
        fg.add_edge(n1, n2, weight=dis)
        fg.add_edge(n2, n1, weight=dis)

    return fg


def main():

    # methods = ['dc', 'dice', 'sokalsneath', 'jaccard', 'lance', 'ochiai']
    # for methods in ['jaccard', 'dc', 'do']:
    for methods in ['jaccard', 'do']:
        result_out = "/home/wbq/1/Bacteria/single_node/Spirochaetes_order/network_out/combined/{}".format(methods)

        for defined_by in ['vector_freq', 'vector', 'module', 'module_freq']:
        # for defined_by in ['vector_freq', 'module_freq']:

            distances_json = result_out + "/{}_genomes_{}_distances.json".format(defined_by, methods)
            dump_to = result_out + "/edgelist/{}_bacterial_{}_network.edgelist".format(defined_by, methods)
            mst_dump_to = result_out + "/edgelist/{}_bacterial_{}_MST.edgelist".format(defined_by, methods)

            load_dumped = False
            compute_mst = True

            if load_dumped:
                fg = nx.read_weighted_edgelist(dump_to)
            else:
                fg = read_raw_graph(distance_json=distances_json) #meta="/var/zjl/Seagate8T/9-workdir/pfam_scan/faa_microbes/clean_module_defined_microbes.tsv"
                nx.write_weighted_edgelist(fg, dump_to)
                print("Dump netwotk to: {}".format(dump_to))

            print("\n+++Netwok info+++")
            print("number_of_nodes: ", fg.number_of_nodes())
            print("number_of_edges: ", fg.number_of_edges())
            print("Graph is connected? -> {}".format(nx.is_connected(fg)))

            if compute_mst:
                print("{} Start computing MST...".format(datetime.now()))
                mst = nx.minimum_spanning_tree(fg)
                # nx.neighbors()
                print("{} Finish computing MST.".format(datetime.now()))
                print("MST number_of_nodes: ", mst.number_of_nodes())
                print("MST number_of_edges: ", mst.number_of_edges())
                print("MST Graph is connected? -> {}".format(nx.is_connected(mst)))
                nx.write_weighted_edgelist(mst, mst_dump_to)
                print("Save MST to {}".format(mst_dump_to))

if __name__ == '__main__':
    main()