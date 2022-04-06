import os
import sys
import multiprocessing as mp
import glob
import datetime as DT
from collections import defaultdict
from hashlib import sha256


class PfamSearch(object):
    """Runs pfam_search.pl over a set of genomes."""

    def __init__(self,
                 threads,
                 pfam_hmm_dir,
                 protein_file_suffix,
                 pfam_suffix,
                 pfam_top_hit_suffix,
                 checksum_suffix,
                 output_dir):
        """Initialization."""

        self.threads = threads

        self.pfam_hmm_dir = pfam_hmm_dir
        self.protein_file_suffix = protein_file_suffix
        self.pfam_suffix = pfam_suffix
        self.pfam_top_hit_suffix = pfam_top_hit_suffix
        self.checksum_suffix = checksum_suffix
        self.output_dir = output_dir

    def _topHit(self, pfam_file):
        """Determine top hits to PFAMs.

        A gene may be assigned to multiple
        PFAM families from the same clan. The
        search_pfam.pl script takes care of
        most of these issues and here the results
        are simply parsed.

        Parameters
        ----------
        pfam_file : str
            Name of file containing hits to PFAM HMMs.
        """

        assembly_dir, filename = os.path.split(pfam_file)
        genome_id = filename.replace(self.pfam_suffix, '')
        output_tophit_file = os.path.join(self.output_dir, genome_id, filename.replace(self.pfam_suffix,
                                                                                       self.pfam_top_hit_suffix))
        tophits = defaultdict(dict)
        for line in open(pfam_file):
            if line[0] == '#' or not line.strip():
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            if gene_id in tophits:
                if hmm_id in tophits[gene_id]:
                    if bitscore > tophits[gene_id][hmm_id][1]:
                        tophits[gene_id][hmm_id] = (evalue, bitscore)
                else:
                    tophits[gene_id][hmm_id] = (evalue, bitscore)
            else:
                tophits[gene_id][hmm_id] = (evalue, bitscore)

        fout = open(output_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, hits in tophits.items():
            hit_str = []
            for hmm_id, stats in hits.items():
                hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
            fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
        fout.close()

        # calculate checksum
        checksum = sha256(output_tophit_file.encode("utf-8")).hexdigest()
        fout = open(output_tophit_file + self.checksum_suffix, 'w')
        fout.write(checksum)
        fout.close()

    def _workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        try:
            while True:
                gene_file = queueIn.get(block=True, timeout=None)
                if gene_file is None:
                    break

                genome_dir, filename = os.path.split(gene_file)
                genome_id = filename.replace(self.protein_file_suffix, '')
                output_hit_file = os.path.join(self.output_dir, genome_id, filename.replace(self.protein_file_suffix,
                                                                                            self.pfam_suffix))
                if not os.path.exists(os.path.join(self.output_dir, genome_id)):
                    os.system("mkdir {}".format(os.path.join(self.output_dir, genome_id)))
                else:
                    print('[Info] Skip {} for its existence.'.format(genome_id))
                    continue

                # dir_path = os.path.dirname(os.path.realpath(__file__))
                pfam_search_script = os.path.join(self.pfam_hmm_dir, 'PfamScan', 'pfam_scan.pl')
                cmd = '%s -outfile %s -cpu %d -fasta %s -dir %s' % (pfam_search_script,
                                                                    output_hit_file,
                                                                    self.cpus_per_genome,
                                                                    gene_file,
                                                                    self.pfam_hmm_dir)
                os.system(cmd)
                # calculate checksum
                checksum = sha256(output_hit_file.encode("utf-8")).hexdigest()
                fout = open(output_hit_file + self.checksum_suffix, 'w')
                fout.write(checksum)
                fout.close()

                # identify top hit for each gene
                self._topHit(output_hit_file)

                queueOut.put(gene_file)
        except Exception as error:
            raise error

    def _writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a is None:
                break

            processedItems += 1
            time_str = DT.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            statusStr = '\n%s ==> Finished processing %d of %d (%.1f%%) genomes.\n' % (time_str,
                                                                                processedItems,
                                                                                numDataItems,
                                                                                float(
                                                                                    processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, gene_files):
        """Annotate genes with Pfam HMMs.

        Parameters
        ----------
        gene_files : iterable
            Gene files in FASTA format to process.
        """

        self.cpus_per_genome = max(1, self.threads / len(gene_files))

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in gene_files:
            workerQueue.put(f)

        for _ in range(self.threads):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self._workerThread, args=(workerQueue, writerQueue)) for _ in
                          range(self.threads)]
            writeProc = mp.Process(target=self._writerThread, args=(len(gene_files), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()
                if p.exitcode == 1:
                    raise ValueError("Pfam Error")

            writerQueue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()
            writeProc.terminate()
            raise

def main():
    threads = 16
    pfam_hmm_dir = "/home/wbq/1/Pfam"
    protein_file_suffix = ".faa"
    pfam_suffix = "_pfam.tsv"
    pfam_top_hit_suffix = "_pfam_tophit.tsv"
    checksum_suffix = ".sha256"
    output_dir = "/home/wbq/1/Bacteria/single_node/Spirochaetes_order/pfam_scan_out"

    pfam_scan = PfamSearch(threads,
                 pfam_hmm_dir,
                 protein_file_suffix,
                 pfam_suffix,
                 pfam_top_hit_suffix,
                 checksum_suffix,
                 output_dir)

    input_dir = "/home/wbq/1/Bacteria/single_node/Spirochaetes_order/faa"
    faa_files = sorted(glob.glob(os.path.join(input_dir, '*'+protein_file_suffix)))
    pfam_scan.run(faa_files)


if __name__ == '__main__':
    main()