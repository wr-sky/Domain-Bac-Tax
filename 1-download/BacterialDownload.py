import pandas as pd
import os
import glob
import urllib.request


# Select data from metadata
out_dir = "/home/wbq/1/Bacteria/Ref"
metainfo = "/home/wbq/1/Bacteria/bacterial_assembly_summary.csv"
ftp_path = 'ftp_path'
metadata = pd.read_csv(metainfo, sep='\t', dtype={ftp_path: str})
ref_rep_category = ('representative genome', 'reference genome')
metadata = metadata[(metadata['assembly_level'] == 'Complete Genome') & (metadata['genome_rep'] == 'Full') & (metadata['version_status'] == 'latest') & (metadata['organism_name'].astype(str).str.contains('sp.')==False) & (metadata['refseq_category'].isin(ref_rep_category))]

# Download data from website according to the metadata information
wget_base = "wget {url} -O {outfile} -e \"https_proxy\"=http://127.0.0.1:8889"

j = 0
for seq_fmt in ["fna"]:
    for file_path in metadata['ftp_path']:
        str_split = file_path.split('/')
        file_path = 'https://' + file_path.split('//')[1]
        url = os.path.join(file_path, "{ac}_genomic.{fmt}.gz".format(ac=str_split[9], fmt=seq_fmt))
        outfile = os.path.join(out_dir, seq_fmt, "{gid}.{fmt}.gz".format(gid=str_split[9], fmt=seq_fmt))
        j += 1
        if j % 100 == 0:
            print("=================== Downloading " + str(j) + "/2587 =========================")
        if not os.path.exists(outfile):
            cmd = wget_base.format(url=url, outfile=outfile)
            os.system(cmd)

print("fna Completed!")

i = 0
for seq_fmt in ["faa"]:
    for file_path in metadata['ftp_path']:
        str_split = file_path.split('/')
        file_path = 'https://' + file_path.split('//')[1]
        url = os.path.join(file_path, "{ac}_protein.{fmt}.gz".format(ac=str_split[9], fmt=seq_fmt))
        outfile = os.path.join(out_dir, seq_fmt, "{gid}.{fmt}.gz".format(gid=str_split[9], fmt=seq_fmt))
        i += 1
        if i % 100 == 0:
            print("=================== Downloading " + str(i) + "/2587 =========================")
        if not os.path.exists(outfile):
            cmd = wget_base.format(url=url, outfile=outfile)
            os.system(cmd)

print("faa Completed!")


