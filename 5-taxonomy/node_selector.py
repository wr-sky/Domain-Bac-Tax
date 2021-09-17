import pandas as pd

input = '/home/wbq/1/Bacteria/taxonomy_out_all/bacterial_assembly_all.csv.new'


output_Spirochaetes = '/home/wbq/1/Bacteria/taxonomy_out_all/bacterial_assembly_Spirochaetes_order.csv'


microbes = pd.read_csv(input, sep='\t', index_col=False, dtype={"GCF": object})

# microbes_Spirochaetes = microbes[(microbes['NCBI_Order'] == 'Leptospirales') & (microbes['genome_rep'] == 'Full') & (microbes['assembly_level'] == 'Complete Genome') & (microbes['version_status'] == 'latest')]
microbes_order = microbes[(microbes['NCBI_Order'] == 'Brachyspirales') | (microbes['NCBI_Order'] == 'Brevinematales')]

#microbes_Planctomycetes.to_csv(output_Planctomycetes, sep='\t', index=False)
microbes_order.to_csv(output_Spirochaetes, sep='\t', index=False)
#microbes_Proteobacteria.to_csv(output_Proteobacteria, sep='\t', index=False)