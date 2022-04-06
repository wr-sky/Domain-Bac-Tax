1-download
    -- BacteriaDonwload.py: Providing bulk download for the selected bacterial genetic sequences in fna format.
    -- gunzip.py: Providing bulk decompress for the downloaded gunzip packages and rename decompressed files with GCF numbers. 

2-checkm
    -- run_checkm.sh: Providing bulk check for the fna files with checkm application
    -- check_completeness.py: Checking the completeness of fna files according to outputs of run_checkm.sh process

3-pfam
    -- pfam_scan.py: Recognizing and collecting the domain information of each species according to pfam database and the faa file.

4-tree
    -- microbes.py:  Calculating the distance of each pair of species according to different statistic and mathematic models
    -- bactrial_mst.py: Building up MST according to the distance network provided by microbes.py

5-taxonomy
    -- ncbi_taxonomy.py: Matching the GCF number with taxonomy categories according to NCBI taxdump database
    -- phyla_itol.py: Creating the file that connects phyla that we focus with different colors so that the corresponding nodes can be colored in Cytoscape. 

6-clustering
    -- MST_counter_NCBI.py：Clustering MST according to NCBI database and collect the clustering result.
    -- MST_counter_GTDB.py：Clustering MST according to GTDB database and collect the clustering result.

7-data
    -- pfam_tophit: The output in '3-pfam' process with each file containing the protein domain information.
    -- json_edgelist: Including two kinds of file format, '*.json' and '*.edgelist'. The json file contains documents recording the distance of each pair of species; The edgelist file contains documents recording the MST structures.
    -- gcf_taxonomy: A document recording the connection between GCF and taxonomy, which can be utilized by Cytoscape to color corresponding nodes. 