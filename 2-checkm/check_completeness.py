import os
import sys

faa_complete_path = '/home/wbq/1/Bacteria/single_node/Spirochaetes_order/faa'
faa_non_complete_path = '/home/wbq/1/Bacteria/single_node/Spirochaetes_order/non-complete_faa'
count = 0
"""
for i in range(14):
    if i == 13:
        path = '/home/wbq/1/Bacteria/Ref/checkm_out/storage/bin_stats_ext.tsv'
    else:
        path = os.path.join('/home/wbq/1/Bacteria/Ref/checkm_out', 'batch{batch_id}'.format(batch_id=i*800), 'storage', 'bin_stats_ext.tsv')
    try:
        if os.path.exists(path):
            file = open(path, 'r')
            print('file: batch{batch_id}.'.format(batch_id=i*800))
            for line in file.readlines():
                tmp1 = line.split('\'Completeness\': ')[1]
                completeness = tmp1.split(', \'Contamination\': ')[0]
                Contamination = tmp1.split(', \'Contamination\': ')[1].split(', \'GC\'')[0]
                genome_id = line[0:9]
                if float(completeness) - 5 * float(Contamination) < 65:
                    print(genome_id)
                    count = count + 1
                    file_non_complete = os.path.join(faa_complete_path, genome_id+'.faa')
                    if os.path.exists(file_non_complete):
                        cmd = 'mv %s %s' % (file_non_complete, faa_non_complete_path)
                        os.system(cmd)
            file.close()
    except Exception as error:
        raise error
"""

path = '/home/wbq/1/Bacteria/single_node/Spirochaetes_order/checkm_out/storage/bin_stats_ext.tsv'
try:
    if os.path.exists(path):
        file = open(path, 'r')
        for line in file.readlines():
            tmp1 = line.split('\'Completeness\': ')[1]
            completeness = tmp1.split(', \'Contamination\': ')[0]
            Contamination = tmp1.split(', \'Contamination\': ')[1].split(', \'GC\'')[0]
            genome_id = line[0:9]
            if float(completeness) - 5 * float(Contamination) < 65:
                print(genome_id)
                count = count + 1
                file_non_complete = os.path.join(faa_complete_path, genome_id + '.faa')
                if os.path.exists(file_non_complete):
                    cmd = 'mv %s %s' % (file_non_complete, faa_non_complete_path)
                    os.system(cmd)
        file.close()
except Exception as error:
    raise error
print(count)
