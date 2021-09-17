import pandas as pd
import os

base_path = '/home/wbq/1/Bacteria/Ref/taxonomy_out/'
ref_microbes = os.path.join(base_path, 'bacterial_assembly_refseq_complete.csv.new')

microbes = pd.read_csv(ref_microbes, sep='\t', index_col=False, dtype={"GCF": object})
phyla_list = set(microbes['NCBI_Phylum'])
for phyla in phyla_list:
    microbes_phyla = microbes[(microbes['NCBI_Phylum'] == phyla)]
    outfile = os.path.join(base_path, 'phyla', phyla + '.csv')
    with open(outfile, 'w') as fout:
        microbes_phyla.to_csv(fout, sep='\t', index=False)