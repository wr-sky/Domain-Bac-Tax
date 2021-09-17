# To find out all the GCF name belonging to the same phyla

# To find out all phylum in the file

phyla_key = []
with open('/home/wbq/1/Bacteria/Ref/taxonomy_out/bacterial_assembly_refseq_complete.csv.new', 'r') as new:
    count = 0
    for i in new.readlines():
        count += 1
        if (count != 1):
            key = i.split('\t')[-7]
            if key not in phyla_key:
                phyla_key.append(key)
    print(phyla_key)


    for phyla in phyla_key:
        file_name = '/home/wbq/1/Bacteria/Ref/taxonomy_out/phyla/' + phyla + '.txt'
        accession_list = []
        try:
            with open(file_name, 'w') as phyla_file:
                count = 0
                with open('/home/wbq/1/Bacteria/Ref/taxonomy_out/bacterial_assembly_refseq_complete.csv.new', 'r') as new:
                    for i in new.readlines():
                        count += 1
                        str_tmp = i.split('\t')[-7]
                        if (count != 1 and str_tmp == phyla):
                            accession_list.append(i.split('\t')[-1])
                    phyla_file.writelines(accession_list)

        except IOError:
            print('File opens error')


