import sys
import os
import glob
import json
import numpy as np
import pandas as pd
import datetime as dt
from typing import Iterable
from Bio import Phylo
from math import sqrt
from collections import defaultdict
import config as Config
import multiprocessing as mp
import math

DEBUG = True
MODULES = ('Family', 'Domain', 'Disordered', ' Coiled-Coil')

class Microbes(object):
    """Identify and define microbes."""

    def __init__(self,
                 pfam_scan_out_dir,
                 output_dir,
                 protein_file_suffix,
                 pfam_suffix,
                 pfam_top_hit_suffix,
                 checksum_suffix,
                 ):
        """Initialize"""
        self.pfam_scan_out_dir = pfam_scan_out_dir
        self.output_dir = output_dir
        self.protein_file_suffix = protein_file_suffix
        self.pfam_suffix = pfam_suffix
        self.pfam_top_hit_suffix = pfam_top_hit_suffix
        self.checksum_suffix = checksum_suffix

    def get_genome_ids(self):
        return sorted(pd.split('/')[-1] for pd in  glob.glob(os.path.join(self.pfam_scan_out_dir, "*")))

    def _top_hits(self, pfam_file):
        """Determine top hits to PFAMs, select by <bit score> value."""

        top_hits = defaultdict(dict)
        with open(pfam_file) as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                line_split  = line.split()
                gene_id     = line_split[0]
                hmm_id      = line_split[5]
                hmm_type    = line_split[7]
                bitscore    = float(line_split[11])

                if gene_id in top_hits:
                    if hmm_id in top_hits[gene_id]:
                        if bitscore > top_hits[gene_id][hmm_id][1]:
                            top_hits[gene_id][hmm_id] = (hmm_type, bitscore)
                    else:
                        top_hits[gene_id][hmm_id] = (hmm_type, bitscore)
                else:
                    top_hits[gene_id][hmm_id] = (hmm_type, bitscore)

        return top_hits

    def _get_hits(self, pfam_file, coverage=0.6, mask=True):
        """Get query hits to PFAMs, saving alignment position and coverage information.

        There may present a below hit result(GCF_000005825):
            --------------------------------------------------------------------------------------------------------------------------------
            WP_012957027.1     10    472     10    472 PF00478.24  IMPDH             Domain     1   345   345    558.1    6e-168   1 CL0036
            WP_012957027.1     97    142     82    146 PF00571.27  CBS               Domain     7    52    57     33.8   3.3e-08   1 No_clan
            WP_012957027.1    153    202    151    207 PF00571.27  CBS               Domain     3    51    57     32.1   1.1e-07   1 No_clan
            --------------------------------------------------------------------------------------------------------------------------------
            WP_012957130.1      6    331      4    331 PF04997.11  RNA_pol_Rpb1_1    Domain     3   312   312    293.4     2e-87   1 No_clan
            WP_012957130.1    333    387    333    396 PF00623.19  RNA_pol_Rpb1_2    Domain     1    55   166     38.2   1.5e-09   1 No_clan
            WP_012957130.1    379    474    375    475 PF00623.19  RNA_pol_Rpb1_2    Domain    70   165   166    100.5   1.1e-28   1 No_clan
            WP_012957130.1    478    663    478    663 PF04983.17  RNA_pol_Rpb1_3    Domain     1   158   158    116.9   7.4e-34   1 No_clan
            WP_012957130.1    692    760    692    776 PF05000.16  RNA_pol_Rpb1_4    Domain     1    85   108     53.4   1.9e-14   1 No_clan
            WP_012957130.1    770   1008    770   1012 PF04998.16  RNA_pol_Rpb1_5    Domain     1   158   267    177.9   2.4e-52   1 No_clan
            WP_012957130.1   1011   1127   1001   1129 PF04998.16  RNA_pol_Rpb1_5    Domain   187   265   267     52.4   4.5e-14   1 No_clan
            --------------------------------------------------------------------------------------------------------------------------------
        Where one protein has multiple Pfam Hits and one or more of these hits alignment
        cover several others' alignment.So
            (1) we order every protein's Pfam hits according to their alignment start po-
            sition and preferentially select those Pfam hits with longer alignment region
            and bit_score next. As a result, the example above will preserve PF00478.24
            as the last Pfam Hit and discard the another two Pfam hits.
            (2) Generally, Pfam does not allow any amino-acid to match more than one Pfam-A
            family, unless the overlapping families are part of the same clan. under this
            situations, we compared the aligned length l1 and l2 of the two overlapped Pfam
            families. If abs(l2-l1)/max(l1, l2) <= threshold(0.15), then we preserve the
            Pfam family with lower bit_score, or preserve the Pfam family with longer aligned
            length.

        Parameter:
        ----------
        coverage: Pfam entry alignment coverage threshold. Default is o.6, any alignment coverage
        below this value will be discard.
        mask: Whether discard the covered Pfam hit entries. Default is True.
        """

        # <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
        hits = defaultdict(dict)
        with open(pfam_file) as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                line_split = line.split()
                gene_id = line_split[0]
                alignment_start = int(line_split[1])
                alignment_end = int(line_split[2])
                hmm_id = line_split[5]
                hmm_type = line_split[7]
                hmm_start = int(line_split[8])
                hmm_end = int(line_split[9])
                hmm_length = int(line_split[10])
                bitscore = float(line_split[11])

                cover = round((hmm_end-hmm_start+1)/hmm_length, 6)
                if cover < coverage:
                    continue

                hits[gene_id][(hmm_id, alignment_start, alignment_end)] = (hmm_type, hmm_length, cover, bitscore)

        if mask:
            masked_hits = defaultdict(dict)
            for gene_id, g_hits in hits.items():
                last_hmm_id = ''
                last_start = -1
                last_end = -1
                last_btscore = 0.0

                for hmm_alignment in sorted(g_hits, key=lambda a: a[1]):
                    info = g_hits[hmm_alignment]
                    hmm_id, a_start, a_end = hmm_alignment

                    if (last_start == -1) or (a_start == last_start and a_end == last_end and info[-1] > last_btscore):
                        last_hmm_id = hmm_id
                        last_start = a_start
                        last_end = a_end
                        last_btscore = info[-1]
                    elif last_end >= a_start >= last_start:
                        if a_end < last_end:
                            continue
                        else:
                            l1 = last_end - last_start + 1
                            l2 = a_end - a_start + 1

                            if ((abs(l2 - l1)/max(l1, l2) <= 0.15) and info[-1] > last_btscore) or \
                                    ((l2 - l1)/l2 > 0.15):
                                last_hmm_id = hmm_id
                                last_start = a_start
                                last_end = a_end
                                last_btscore = info[-1]
                            else:
                                continue
                    else:
                        hkey = (last_hmm_id, last_start, last_end)
                        if  hkey not in masked_hits[gene_id]:
                            masked_hits[gene_id][hkey] = hits[gene_id][hkey]
                        last_hmm_id = hmm_id
                        last_start = a_start
                        last_end = a_end
                        last_btscore = info[-1]
                else:
                    hkey = (last_hmm_id, last_start, last_end)
                    if last_hmm_id and hkey not in masked_hits[gene_id]:
                        masked_hits[gene_id][hkey] = hits[gene_id][hkey]
            else:
                hits = masked_hits

        return masked_hits

    def _define_microbe_by_module_set(self, pfam_file, module="Domain"):
        """Define a microbe using `module` characteristics.

        A CDS/gene in bacterial genome may hit multiple Pfam entries at the `module` class.
        Hence the `module` entry with highest bit_score will be selected. THe microbe will
        be defined as a set of functional `module` entries.

        Parameters
        ----------
        pfam_file: pfam_scan.pl output.
        module: Pfam entry level, which is one of {Family, Domain, Coiled-Coil, Disordered}; Repeat, Motifs}
        """

        microbe = set()
        MODULES = ('Family', 'Domain', 'Disordered', ' Coiled-Coil')

        good_hits = self._get_hits(pfam_file, coverage=0.6)
        for gene_id, hits in good_hits.items():

            for hmm_alignment in sorted(hits.keys(), key=lambda k: k[1]):
                hmm_id, a_start, a_end = hmm_alignment
                hmm_type = hits[hmm_alignment][0]
                if hmm_type in MODULES:
                    microbe.add(hmm_id)

        return microbe

    def _define_microbe_by_module_freq_set(self, pfam_file):

        microbe = []
        MODULES = ('Family', 'Domain', 'Disordered', ' Coiled-Coil')

        good_hits = self._get_hits(pfam_file, coverage=0.6)
        for gene_id, hits in good_hits.items():

            for hmm_alignment in sorted(hits.keys(), key=lambda k: k[1]):
                hmm_id, a_start, a_end = hmm_alignment
                hmm_type = hits[hmm_alignment][0]
                if hmm_type in MODULES:
                    microbe.append(hmm_id)

        return microbe

    """
    def _define_microbe_by_module_freq_set(self, pfam_file):

        microbe = defaultdict(int)
        modules = ('Family', 'Domain')

        good_hits = self._get_hits(pfam_file, coverage=0.6)
        for gene_id, hits in good_hits.items():

            for hmm_alignment in sorted(hits.keys(), key=lambda k: k[1]):
                hmm_id, a_start, a_end = hmm_alignment
                hmm_type = hits[hmm_alignment][0]
                if hmm_type in modules:
                    microbe[hmm_id] += 1

        return ('{}_{}'.format(pf, freq) for pf, freq in microbe.items())
    """

    def _define_microbes_by_vector_set(self, pfam_file):
        """Define a microbe by its genome's protein vector set. If a protein(P)'s sequence has two Pfam entries
        hits PF1 and PF2, then we use <PF1, PF2> indicate P. The set of all protein vectors is used to represent
        this microbe."""

        microbe = set()

        good_hits = self._get_hits(pfam_file, coverage=0.6)
        for gene_id, hits in good_hits.items():
            protein = []

            for hmm_alignment in sorted(hits.keys(), key=lambda k: k[1]):
                hmm_id, a_start, a_end = hmm_alignment
                hmm_type = hits[hmm_alignment][0]
                if hmm_type in MODULES:
                    protein.append(hmm_id)
            if protein:
                microbe.add('_'.join(protein))

        return microbe

    def _define_microbes_by_freq_vector_set(self, pfam_file):
        """Define a microbe by its genome's protein vector set. If a protein(P)'s sequence has two Pfam entries
        hits PF1 and PF2, then we use <PF1, PF2> indicate P. The set of all protein vectors is used to represent
        this microbe."""

        microbe = []

        good_hits = self._get_hits(pfam_file, coverage=0.6)
        for gene_id, hits in good_hits.items():
            protein = []

            for hmm_alignment in sorted(hits.keys(), key=lambda k: k[1]):
                hmm_id, a_start, a_end = hmm_alignment
                hmm_type = hits[hmm_alignment][0]
                if hmm_type in MODULES:
                    protein.append(hmm_id)
            if protein:
                microbe.append('_'.join(protein))

        return microbe

    """
    def _define_microbes_by_freq_vector_set(self, pfam_file):
        # Define a microbe by its genome's protein functional module vector and its frequency. If protein <PF1, PF2>
        # appears 3 times, then add e entry `<PF1, PF2, 3>` to its genome set. Then the set of al protein and their freq
        # vectors is used to represent this genome.


        microbe = defaultdict(int)

        good_hits = self._get_hits(pfam_file, coverage=0.6)
        for gene_id, hits in good_hits.items():
            protein = []

            for hmm_alignment in sorted(hits.keys(), key=lambda k: k[1]):
                hmm_id, a_start, a_end = hmm_alignment
                hmm_type = hits[hmm_alignment][0]
                if hmm_type in MODULES:
                    protein.append(hmm_id)
            if protein:
                microbe['_'.join(protein)] += 1

        return ('{}_{}'.format(pv, freq) for pv, freq in microbe.items())
    """

    def define_microbes(self, by="module", output=None, selected_ids=None):
        """According to functional `module`, define all microbes under `pfam_scan_out_dir`(Pfam_scan output
        directory) and save them into a json file.

        Parameter
        ---------
        by: How to define a microbe. One of {'module', 'module_freq', 'vector', 'vector_freq'}. Default is `module` set.
        selected_ids: A list of GCF ids which was selected to define and analysis.
        """

        genome_ids = self.get_genome_ids()
        if selected_ids and isinstance(selected_ids, (list, tuple, set)):
            genome_ids = tuple(set(genome_ids).intersection(set(selected_ids)))

        # genome_ids = genome_ids[:80]

        microbes = dict()
        for i, genome_id in enumerate(genome_ids, 1):
            if DEBUG and (i % 100 == 0):
                print("\t[Debug] Finished processing {} of {} ({}%) genomes".format(i, len(genome_ids),
                                                                                    round(i*100/len(genome_ids), 2)))

            pfam_file = os.path.join(self.pfam_scan_out_dir, genome_id, genome_id+self.pfam_suffix)
            if by == 'module':
                microbe = tuple(self._define_microbe_by_module_set(pfam_file))
            elif by == 'module_freq':
                microbe = tuple(self._define_microbe_by_module_freq_set(pfam_file))
            elif by == 'vector':
                microbe = tuple(self._define_microbes_by_vector_set(pfam_file))
            elif by == 'vector_freq':
                microbe = tuple(self._define_microbes_by_freq_vector_set(pfam_file))
            else:
                raise ValueError("Parameter `by` must be one of {'module', 'vector'}.")
            microbes[genome_id] = microbe

        outfile = output or os.path.join(self.output_dir, "{}_defined_microbes.tsv".format(by))
        with open(outfile, "w") as fout:
            fout.write("{}\t{}\t{}\n".format("GCF", "Size", "Modules"))
            for genome_id, microbe in microbes.items():
                fout.write("{}\t{}\t{}\n".format(genome_id,
                                                 len(microbe),
                                                 ','.join(microbe)))
        print("[Debug] Save defined microbes to: {}".format(outfile))
        return outfile

    def list_intersection(self, listA, listB):

        setA = set(listA)
        setB = set(listB)

        setAandB = setA.intersection(setB)
        count = 0
        for i in setAandB:
            count += min(listA.count(i), listB.count(i))
        return count

        """
        number = 0
        for i in listA:
            count = 0
            for j in listB:
                if i == j:
                    number += 1
                    del listB[count]
                    break
                else:
                    count += 1
        return number
        """

    def compute_distance(self, defined_microbes, method, process_num, p_rank, output):
        microbes = pd.read_csv(defined_microbes, sep='\t', index_col=False, dtype={"GCF": object})
        distances = dict()
        for i in range(microbes.shape[0]):
            if i % process_num == p_rank:
                if DEBUG:
                    print("[Debug]i = {}/{} compute Jaccard distance...".format(i, microbes.shape[0]))
                gcfA = microbes.iloc[i]['GCF']
                sizeA = microbes.iloc[i]['Size']
                microbeA = list(microbes.iloc[i]['Modules'].split(','))
                assert sizeA == len(microbeA), "[A]microbe {} genome size conflict.".format(gcfA)

                for j in range(i + 1, microbes.shape[0]):
                    gcfB = microbes.iloc[j]['GCF']
                    sizeB = microbes.iloc[j]['Size']
                    microbeB = list(microbes.iloc[j]['Modules'].split(','))
                    assert sizeB == len(microbeB), "[B]microbe {} genome size conflict.".format(gcfB)

                    sizeOfAandB = self.list_intersection(microbeA, microbeB)
                    sizeOfAorB = sizeA + sizeB - sizeOfAandB
                    if method == 'jaccard':
                        jd = round((sizeOfAorB - sizeOfAandB) / sizeOfAorB, 6)
                    elif method == 'do':
                        jd = round(sqrt(math.log(sizeOfAandB / sizeA) * math.log(sizeOfAandB / sizeB)), 6)
                    else: # dc
                        smaller = sizeA if sizeA < sizeB else sizeB
                        jd = round((smaller - sizeOfAandB)/smaller, 6)
                    distances['-'.join((gcfA, gcfB))] = jd

        output = output.split('.')[0] + '_' + str(p_rank) + '.json'
        with open(output, "w") as fout:
            json.dump(distances, fout)


    def compute_jaccard_distances(self, by, defined_microbes, method, output, redefine_microbes=False, sep='-'):
        """Compute genomes' Jaccard distances and save them into a json file.
        Jaccard similarity: J(A, B) = |A and B| / |A or B|.
        Hence Jaccard distance jd(A, B) = 1 - J(A, B) = (|A or B| - |A and B|) / (|A or B|).

        Parameter
        ---------
        defined_microbes: The file path of defined microbes by functional modules.
        redefine_microbes:  Whether rerun self.define_microbes to define all genomes or not. Default is False.
        output: distance information output target.
        """

        if not os.path.exists(defined_microbes) or redefine_microbes:
            defined_microbes = self.define_microbes(defined_microbes)

        distances = dict()
        if DEBUG:
            print("[Debug] Start computing {} distance...".format(method))

        if by in ['vector_freq', 'module_freq']:
            process_num = 16
            process_list = []
            for p_rank in range(process_num):
                proc = mp.Process(target=self.compute_distance, args=(defined_microbes, method, process_num, p_rank, output))
                process_list.append(proc)
                proc.start()
            for p in process_list:
                p.join()
            """
            for i in range(process_num):
                input_json = output.split('.')[0] + '_' + str(i) + '.json'
                with open(input_json, "r") as fin:
                    distances.update(json.load(fin))
            with open(output, "w") as fout:
                json.dump(distances, fout)
            """
            print("[Debug] Save Jaccard distances information to {}\t{}".format(output, dt.datetime.now()))

        else:
            microbes = pd.read_csv(defined_microbes, sep='\t', index_col=False, dtype={"GCF": object})
            for i in range(microbes.shape[0]):
                if DEBUG:
                    print("[Debug]i = {}/{} compute Jaccard distance...".format(i, microbes.shape[0]))
                gcfA = microbes.iloc[i]['GCF']
                sizeA = microbes.iloc[i]['Size']
                microbeA = set(microbes.iloc[i]['Modules'].split(','))
                assert sizeA == len(microbeA), "[A]microbe {} genome size conflict.".format(gcfA)
                for j in range(i+1, microbes.shape[0]):
                    gcfB = microbes.iloc[j]['GCF']
                    sizeB = microbes.iloc[j]['Size']
                    microbeB = set(microbes.iloc[j]['Modules'].split(','))
                    assert sizeB == len(microbeB), "[B]microbe {} genome size conflict.".format(gcfB)
                    sizeOfAandB = len(microbeA.intersection(microbeB))
                    sizeOfAorB = sizeA + sizeB - sizeOfAandB
                    if method == 'jaccard':
                        jd = round((sizeOfAorB - sizeOfAandB) / sizeOfAorB, 6)
                    elif method == 'do':
                        jd = round(sqrt(math.log(sizeOfAandB / sizeA) * math.log(sizeOfAandB / sizeB)), 6)
                    else: # dc
                        smaller = sizeA if sizeA < sizeB else sizeB
                        jd = round((smaller - sizeOfAandB)/smaller, 6)
                    distances[sep.join((gcfA, gcfB))] = jd

            with open(output, "w") as fout:
                json.dump(distances, fout)
            print("[Debug] Save Jaccard distances information to {}\t{}".format(output, dt.datetime.now()))

    def compute_pair_distance(self, defined_microbes, output, method, redefine_microbes=False, sep='-'):
        """
         OTU:       Sample B    1   0
         sample A
           1                    a   b
           0                    c   d

        method:
          1.dc (domain content) Compute genomes' proteins domain distances and save them into a json file.
                            D = A' / (A'+ AB)
            A' is number of unique superfamily folds in the smaller of two genomes, A and B,
            and AB is the number of superfamily folds they share.
            [REF: Phylogeny determined by protein domain content(2005)]
          2. dice
            s = 2a / (2a + b + c)
          3. sokalsneath
            s = a / (a + 2b + 2c)
          4. jaccard
            s = a / (a + b + c)
          5. lance
            D = (b+c) / (2a+b+c)
          6.ochiai
            s = a / sqrt((a + b)(a + c))
          7. domain organization (do)
            D = sqrt((-ln(a/(a+b)))*(-ln(a/(a+c))))

          D = sqrt(1 - s)

        Parameter
        ---------
        defined_microbes: The file path of defined microbes by functional modules.
        redefine_microbes:  Whether rerun self.define_microbes to define all genomes or not. Default is False.
        output: distance information output target.
        """

        if not os.path.exists(defined_microbes) or redefine_microbes:
            defined_microbes = self.define_microbes(defined_microbes)

        microbes = pd.read_csv(defined_microbes, sep='\t', index_col=False, dtype={"GCF": object})
        distances = dict()

        if DEBUG:
            print("[Debug] Start computing pair distance...")
        for i in range(microbes.shape[0]):
            if DEBUG and (i % 100 == 0):
                print("[Debug] i = {}/{} compute pair [{}] distance...".format(i, microbes.shape[0], method))

            gcfA = microbes.iloc[i]['GCF']
            sizeA = microbes.iloc[i]['Size']
            microbeA = set(microbes.iloc[i]['Modules'].split(','))
            assert sizeA == len(microbeA), "[A]microbe {} genome size conflict.".format(gcfA)
            for j in range(i + 1, microbes.shape[0]):
                gcfB = microbes.iloc[j]['GCF']
                sizeB = microbes.iloc[j]['Size']
                microbeB = set(microbes.iloc[j]['Modules'].split(','))
                assert sizeB == len(microbeB), "[B]microbe {} genome size conflict.".format(gcfB)

                if method == 'dc':
                    sizeOfAandB = len(microbeA.intersection(microbeB))
                    sizeOfA_pie = len(microbeA.difference(microbeB)) if sizeA <= sizeB else len(
                        microbeB.difference(microbeA))
                    D = round(sizeOfA_pie / (sizeOfA_pie + sizeOfAandB), 6)
                else:
                    a = len(microbeA.intersection(microbeB))
                    b = len(microbeA.difference(microbeB))
                    c = len(microbeB.difference(microbeA))
                    if method == 'ochiai':
                        sim = a / sqrt((a + b) * (a + c))
                        D = round(sqrt(1 - sim), 6)
                    elif method == 'sokalsneath':
                        sim = a / (a + 2 * b + 2 * c)
                        D = round(sqrt(1 - sim), 6)
                    elif method == 'dice':
                        sim = 2 * a / (2 * a + b + c)
                        D = round(sqrt(1 - sim), 6)
                    elif method == 'lance':
                        D = round((b + c) / (2 * a + b + c), 6)
                    else:
                        # jaccard
                        sim = a / (a + b + c)
                        D = round(sqrt(1 - sim), 6)

                distances[sep.join((gcfA, gcfB))] = D

            if i % 50 == 1:
                with open(output + ".tmp", "w") as fout:
                    json.dump(distances, fout)

        with open(output, "w") as fout:
            json.dump(distances, fout)
        print("[Debug] Save pair distances information to {}\t{}".format(output, dt.datetime.now()))


def main():
    pfam_scan_our_dir = "/home/wbq/1/faa_scan_out"
    output_dir = "/home/wbq/1/network_out"

    microbes = Microbes(pfam_scan_out_dir=pfam_scan_our_dir,
                        output_dir=output_dir,
                        protein_file_suffix=Config.PROTEIN_FILE_SUFFIX,
                        pfam_suffix=Config.PFAM_SUFFIX,
                        pfam_top_hit_suffix=Config.PFAM_TOP_HIT_SUFFIX,
                        checksum_suffix=Config.CHECKSUM_SUFFIX)

    # for define_by in ['vector_freq', 'vector', 'module', 'module_freq']:
    for define_by in ['vector', 'module']:
        print("\n\tRun microbes define_by: [{}]".format(define_by))
        defined_microbes = output_dir + "/{}_defined_microbes.tsv".format(define_by)
        if not os.path.exists(defined_microbes):
            microbes.define_microbes(by=define_by, output=defined_microbes)
        #_thread.start_new_thread(microbes.compute_jaccard_distances, (define_by, defined_microbes, jaccard_distances, False))
        for define_method in ['jaccard']:
        # for define_method in ['jaccard', 'dc', 'do']:
            jaccard_distances = output_dir + "/{}/{}_genomes_{}_distances.json".format(define_method, define_by, define_method)
            microbes.compute_jaccard_distances(by=define_by, defined_microbes=defined_microbes, method=define_method, output=jaccard_distances, redefine_microbes=False)

if __name__ == '__main__':
    main()

