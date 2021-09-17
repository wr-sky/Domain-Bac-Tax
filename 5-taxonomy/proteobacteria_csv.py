import pandas as pd
import os

base_path = '/home/wbq/1/Bacteria/Ref/taxonomy_out/phyla'
ref_microbes = os.path.join(base_path, 'Proteobacteria.csv')

microbes = pd.read_csv(ref_microbes, sep='\t', index_col=False, dtype={"GCF": object})
class_list = set(microbes['NCBI_Class'])
for class_ex in class_list:
    microbes_phyla = microbes[(microbes['NCBI_Class'] == class_ex)]
    outfile = os.path.join(base_path, 'proteobacteria', class_ex + '.csv')
    with open(outfile, 'w') as fout:
        microbes_phyla.to_csv(fout, sep='\t', index=False)