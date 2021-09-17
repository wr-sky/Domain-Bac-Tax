import json

for define_by in ['vector_freq', 'module_freq']:
    for dis_method in ['jaccard', 'do']:
        base_dir = "/home/wbq/1/Bacteria/single_node/Spirochaetes_order/network_out/{}".format(dis_method)
        output = base_dir + "/{}_genomes_{}_distances.json".format(define_by, dis_method)
        process_num = 16
        distances = dict()
        for i in range(process_num):
            input_json = output.split('.')[0] + '_' + str(i) + '.json'
            with open(input_json, "r") as fin:
                distances.update(json.load(fin))
        with open(output, "w") as fout:
            json.dump(distances, fout)