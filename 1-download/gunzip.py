import pandas as pd
import os
import glob

#for target in ['fna', 'faa']:
for target in ['faa']:
    out_dir = "/home/wbq/1/Bacteria/Ref"
    faa_files = sorted(glob.glob(os.path.join(out_dir, target, '*')))
    for i in faa_files:
        print(i)
        cmd = "gunzip -c " + i + " > " + out_dir + "/" + target + "_gunzip/" + i.split('/')[7].split('_')[1].split('.')[0] + "." + target
        os.system(cmd)
#001941465