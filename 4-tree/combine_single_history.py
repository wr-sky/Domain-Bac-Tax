import json

for define_by in ['vector_freq', 'vector', 'module', 'module_freq']:
    for dis_method in ['jaccard', 'do']:
        input_dir = "/home/wbq/1/Bacteria/single_node/Spirochaetes_order/network_out/{}/{}_genomes_{}_distances.json".format(dis_method, define_by, dis_method)
        input_dir_history = "/home/wbq/1/network_out/{}/{}_genomes_{}_distances.json".format(dis_method, define_by, dis_method)
        input_dir_history_2 = "/home/wbq/1/Bacteria/single_node/network_out/{}/{}_genomes_{}_distances.json".format(dis_method, define_by, dis_method)
        output = "/home/wbq/1/Bacteria/single_node/network_out/combined/{}/{}_genomes_{}_distances.json".format(dis_method, define_by, dis_method)
        distances = dict()
        for i in {input_dir, input_dir_history, input_dir_history_2}:
            with open(i, "r") as fin:
                distances.update(json.load(fin))
        with open(output, "w") as fout:
            json.dump(distances, fout)