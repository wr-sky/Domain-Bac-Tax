import pandas as pd
import os

# Select data from metadata
metainfo = "/mnt/hgfs/share/1/bacterial_assembly_summary.csv"
ftp_path = 'fpt_path'
metadata = pd.read_csv(metainfo, sep='\t', dtype={ftp_path: str})
metadata = metadata[(metadata['assembly_level'] == 'Complete Genome') & (metadata['genome_rep'] == 'Full') & (metadata['version_status'] == 'latest') & (metadata['organism_name'].astype(str).str.contains('sp.')==False)]

# Download data from website according to the metadata information
wget_base = "wget {url} -O {outfile}"
out_dir = "/mnt/hgfs/share/1/Bacteria"

for seq_fmt in ["fna", "faa"]:
    for file_path in metadata['ftp_path']:
        str_split = file_path.split('/')
        url = os.path.join(file_path, "{ac}_genomic.{fmt}.gz".format(ac=str_split[9], fmt=seq_fmt))
        outfile = os.path.join(out_dir, seq_fmt, "{gid}.{fmt}.gz".format(gid=str_split[9], fmt=seq_fmt))
        cmd = wget_base.format(url=url, outfile=outfile)
        os.system(cmd)

