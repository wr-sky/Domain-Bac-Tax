import os
import sys

fna_complete_path = '/home/wbq/1/Bacteria/fna/complete_fna'
fna_non_complete_path = '/home/wbq/1/Bacteria/fna/non-complete_fna'
path = '/home/wbq/1/Bacteria/fna/checkm_out/storage/bin_stats_ext.tsv'

count = 0
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
                file_non_complete = os.path.join(fna_complete_path, genome_id + '.fna')
                if os.path.exists(file_non_complete):
                    cmd = 'mv %s %s' % (file_non_complete, fna_non_complete_path)
                    os.system(cmd)
        file.close()
except Exception as error:
    raise error
print(count)
