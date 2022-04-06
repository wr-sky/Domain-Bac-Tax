from collections import defaultdict


def find_mst():
    phylum_dict = dict()
    ref = "D:\\1-研究课题\\1-基于功能模块迁移的细菌重组网络研究\\2-过程结果\\1-taxonomy_out\\bacterial_assembly_refseq_complete_with_6_archaea.csv"
    with open(ref, "r") as rin:
        for line in rin.readlines():
            # 对于输入文件1，可以自行选择一列作为关键字对MST图进行分类
            gcf = line.split('\t')[-1].split('\n')[0]
            phylum = line.split('\t')[-7]
            phylum_dict[gcf] = phylum

    for method in ('jaccard', 'do', 'dc'):
        for model in ('module', 'vector', 'module_freq', 'vector_freq'):
            print('################## ' + method + ' & ' + model + ' ###################')
            # 支持全连接和edgelist格式
            file = "D:\\1-研究课题\\1-基于功能模块迁移的细菌重组网络研究\\3-写文参考\\paper数据\\{}\\edgelist\\{}_bacterial_{}_MST.edgelist".format(method, model, method)
            connection = defaultdict(list)
            with open(file, "r") as fin:
                for line in fin.readlines():
                    a = line.split()[0]
                    b = line.split()[1]
                    connection[a].append(b)
                    connection[b].append(a)

            collection = defaultdict(list)
            used = list()
            key_set = connection.keys()
            for key in key_set:
                if key not in used:
                    used.append(key)
                    phylum_ref = phylum_dict[key]
                    count = 1
                    tmp = connection[key]
                    while len(tmp) != 0:
                        tmp_2 = list()
                        for tmp_iter in tmp:
                            if phylum_dict[tmp_iter] == phylum_ref:
                                count += 1
                                used.append(tmp_iter)
                                for tmp_2_iter in connection[tmp_iter]:
                                    if tmp_2_iter not in used:
                                        tmp_2.append(tmp_2_iter)
                        tmp = tmp_2
                    collection[phylum_ref].append(count)

            iso_num = 0
            phy_num = 0
            weighted_per = 0
            total = 0
            for i in sorted(collection.keys()):
                # print(i, end='; ')
                # print(str(len(collection[i])), end='; ')
                # print(collection[i])
                max_value = max(collection[i])
                phy_total = sum(collection[i])
                iso_num = phy_total - max_value
                total += iso_num
                weighted_per += round(iso_num / phy_total, 3)
                if len(collection[i]) > 1:
                    phy_num += 1
                # print('ios_num is: ' + str(iso_num))
            print('phy_num is: ' + str(phy_num))
            print('total is: ' + str(total))
            print('weighted_per is: ' + str(round(weighted_per / 31, 3)))


if __name__ == '__main__':
    find_mst()
