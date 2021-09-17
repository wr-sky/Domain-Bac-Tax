import glob
import pandas as pd
import numpy as np
import os
from datetime import datetime as dt

def update_taxdump_files(taxdump_path):
    pass

class NCBItax(object):
    """Mapping NCBI taxid into full lineage information."""

    def __init__(self,
                 taxdump_path):
        """
        :param taxdump_path: the path of NCBI Taxonomy database dump files
        """

        self.taxdump_path = taxdump_path
        self.TAXONOMIC_RANKS = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
        self.nodes = None
        self.names = None
        self.merged = None
        self.delnodes = None

        self._load_names()
        self._load_nodes()
        self._load_merged()
        self._load_delnodes()

    def _check_dump_file(self, file_tp):
        """
        Check whether the target file exist or not.If not, raise a error.
        :param file_tp: which dmp file, eg: nodes, names
        :return: target file
        """

        tfile = os.path.join(self.taxdump_path, file_tp + ".dmp")
        if os.path.isfile(tfile) and os.path.exists(tfile):
            return tfile
        else:
            raise ValueError("[{}] not exist.".format(tfile))

    def _load_merged(self):
        """Load merged.dmp file into a DataFrame"""

        merged_file = self._check_dump_file('merged')
        self.merged = pd.read_csv(merged_file, sep='|', header=None, index_col=False,
                         names=[
                             'old_tax_id',
                             'new_tax_id'
                         ])
        self.merged.set_index('old_tax_id', inplace=True)

    def _load_delnodes(self):
        "Load deleted nodes (nodes that existed but were deleted) file into a DataFrame"

        delnodes_file = self._check_dump_file('delnodes')
        delnodes = pd.read_csv(delnodes_file, sep='|', header=None, index_col=False,
                                 names=[
                                     'tax_id'
                                 ])
        self.delnodes = np.array(delnodes)

    def _load_names(self):
        """load names.dmp into a DataFrame"""

        names_file = self._check_dump_file("names")
        df = pd.read_csv(names_file, sep='|', header=None, index_col=False,
                                 names=[
                                     'tax_id',
                                     'name_txt',
                                     'unique_name',
                                     'name_class'
                                 ],
                                 dtype={
                                     'tax_id': int,
                                     'name_txt': object,
                                     'unique_name': object,
                                     'name_class': object
                                 })
        df['name_txt'] = df['name_txt'].apply(str.strip)
        df['unique_name'] = df['unique_name'].apply(str.strip)
        df['name_class'] = df['name_class'].apply(str.strip)

        self.names = df[df['name_class'] == 'scientific name']
        self.names.reset_index(drop=True, inplace=True)
        self.names.set_index('tax_id', inplace=True)

    def _load_nodes(self):
        """Load nodes.dmp into a DataFrame"""

        nodes_file = self._check_dump_file('nodes')
        df = pd.read_csv(nodes_file, sep='|', header=None, index_col=False,
                         names=[
                             'tax_id',
                             'parent_tax_id',
                             'rank',
                             'embl_code',
                             'division_id',
                             'inherited_div_flag',
                             'genetic_code_id',
                             'inherited_GC__flag',
                             'mitochondrial_genetic_code_id',
                             'inherited_MGC_flag',
                             'GenBank_hidden_flag',
                             'hidden_subtree_root_flag',
                             'comments'
                        ],
                         dtype={
                             'rank': object,
                             'embl_code':object,
                             'comments': object
                         })
        df['rank'] = df['rank'].apply(str.strip)
        df['embl_code'] = df['embl_code'].apply(str.strip)
        df['comments'] = df['comments'].apply(str.strip)

        self.nodes = df
        self.nodes.set_index('tax_id', inplace=True)

    def get_lineage_by_id(self, tax_id):
        """Get the genome's full lineage information by its tax_id."""

        lineage = dict()
        while tax_id != 1:
            try:
                parent_tax_id = self.nodes.loc[tax_id]['parent_tax_id']
            except KeyError:
                if tax_id in self.delnodes:
                    print("node {} was deleted.".format(tax_id))
                    break
                exist = False
                while tax_id in self.merged.index:
                    exist = True
                    tax_id = self.merged.loc[tax_id]['new_tax_id']
                if not exist:
                    print("node {} isn't existed.".format(tax_id))
                    break
                else:
                    parent_tax_id = self.nodes.loc[tax_id]['parent_tax_id']
            if parent_tax_id == 1: #taxonomy tree root
                break

            rank = self.nodes.loc[tax_id]['rank']
            if rank in self.TAXONOMIC_RANKS:
                lineage[rank] = self.names.loc[tax_id]['name_txt']

            tax_id = self.nodes.loc[tax_id]['parent_tax_id']

        standardized = []
        for rank in self.TAXONOMIC_RANKS:
            if rank == 'superkingdom':
                prefix = 'd'
            else:
                prefix = rank[0]
            standardized.append('{}__{}'.format(prefix, lineage.get(rank, '')))

        return ';'.join(standardized)

    def get_full_lineage(self, tax_ids, out_file):
        """
        Get full lineage information of the input taxids, and save result to csv file.
        :param tax_ids: tax_id list
        :param out_file: lineage information file path.
        :return:
        """

        cols = ('taxid', 'ncbi_taxonomy')
        result = pd.DataFrame(columns=cols)
        got_ids = set()
        for i, tax_id in enumerate(tax_ids):
            if i % 200 == 0:
                print("==={}% ({}/{})....".format(i*100/len(tax_ids), i, len(tax_ids)))

            if tax_id in got_ids:
                continue

            got_ids.add(tax_id)

            lineage = self.get_lineage_by_id(tax_id)
            result = result.append(pd.Series((tax_id, lineage), index=cols), ignore_index=True)

        result.to_csv(out_file, sep='\t', index=False)
        print("[Save taxid2lineage information to] {}".format(out_file))
        return result


def deal_NCBI_taxonomy4draft_genome():
    """
    Add NCBI taxonomy information for all RefSeq draft genomes.
    :return:
    """

    # assembly_summary = '/var/zjl/Seagate8T/NCBI_GCA/metainfo/refseq/refseq_bacteria_assembly_summary_20190308.tsv'
    # taxid2lineage = '/var/zjl/Seagate8T/NCBI_GCA/metainfo/refseq/draft_genomes_taxid2lineage.tsv'

    taxonomy_out = '/home/wbq/1/Bacteria/taxonomy_out'

    assembly_summary = '/home/wbq/1/Bacteria/bacterial_assembly_summary.csv'
    # non_complete_faa = '/home/wbq/1/Bacteria/Ref/non-complete_faa'

    taxid2lineage = os.path.join(taxonomy_out, 'bacterial_assembly_taxid2lineage.tsv')
    taxdump_path = "/home/wbq/1/Taxonomy/taxdmp"

    # non_complete_faa_GCF = sorted(pd.split('.')[0].split('/')[-1] for pd in glob.glob(os.path.join(non_complete_faa, "*")))
    # print(non_complete_faa_GCF)

    metadata = pd.read_csv(assembly_summary, sep='\t', index_col=False)

    # for name in non_complete_faa_GCF:
    #    for i in ['1', '2', '3', '4']:
    #       metadata = metadata[metadata['assembly_accession'] != ('GCF_' + name + '.' + i)]
    # metadata = metadata[(metadata['assembly_level'] != 'Complete Genome')
    #                     & (metadata['genome_rep'] == 'Full')
    #                     & (metadata['version_status'] == 'latest')]
    metadata.to_csv(taxonomy_out + '/bacterial_assembly_refseq_complete.csv', sep='\t', index=False)
    metadata.reset_index(drop=True, inplace=True)

    reRun = True
    if os.path.exists(taxid2lineage):
        lineages = pd.read_csv(taxid2lineage, sep='\t', index_col=False)
        reRun = len(set(metadata['taxid'])) != lineages.shape[0]
    if reRun:
        ncbi_tax = NCBItax(taxdump_path)
        lineages = ncbi_tax.get_full_lineage(set(metadata['taxid']), taxid2lineage)

    assert lineages.shape[0] == len(set(metadata['taxid'])), "Row Count of Lineages and metadata must be eqaul."
    metadata = pd.merge(metadata, lineages, on='taxid')

    cols = ['assembly_accession', 'ncbi_taxonomy', 'NCBI_Domain', 'NCBI_Phylum', 'NCBI_Class', 'NCBI_Order', 'NCBI_Family', 'NCBI_Genus',
            'NCBI_Species', 'GCF']
    notation = pd.DataFrame(columns=cols)

    for i, row in metadata.iterrows():
        if i % 2000 == 0:
            print("\t[Debug] finishing {} ({}/{}) {}".format(i/metadata.shape[0], i, metadata.shape[0], dt.now()))
        nt = row['ncbi_taxonomy']
        acc = row['assembly_accession']
        for ch in "()-/:,'":
            nt = nt.replace(ch, '_')
        d, p, c, o, f, g, sp = [tx.split('__')[-1].replace(' ', '_') for tx in nt.split(';')]

        notation = notation.append(pd.Series((acc,
                                              nt,
                                              d or 'norank',
                                              p or 'norank',
                                              c or 'norank',
                                              o or 'norank',
                                              f or 'norank',
                                              g or 'norank',
                                              sp or 'norank',
                                              acc.split('_')[-1].split('.')[0]), index=cols),
                                   ignore_index=True)

    print("Draft genome metadata shape(Before merge ncbi_taxonomy): ", metadata.shape)
    metadata = pd.merge(metadata, notation, on='assembly_accession')
    print("Draft genome metadata shape(After merge ncbi_taxonomy): ", metadata.shape)
    metadata.to_csv(taxonomy_out + '/bacterial_assembly_refseq_complete.csv.new', sep='\t', index=False)
    print("[Save] draft genomes assembly summary to: {}".format(assembly_summary + '.new'))

if __name__ == '__main__':
    deal_NCBI_taxonomy4draft_genome()