from collections import defaultdict
from traceback import print_exc
import sys


class KmerIndexer:
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, known_strings, kmer_size=6):
        self.barcode_list = list(known_strings)
        self.k = kmer_size
        self.index = defaultdict(list)
        self._index()

    def _get_kmers(self, seq):
        if len(seq) < self.k:
            return
        kmer = seq[:self.k]
        yield kmer
        for i in range(self.k, len(seq)):
            kmer = kmer[1:] + seq[i]
            yield kmer

    def _index(self):
        for i, barcode in enumerate(self.barcode_list):
            for kmer in self._get_kmers(barcode):
                self.index[kmer].append(i)

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=10, min_kmers=2):
        barcode_counts = defaultdict(int)
        for kmer in self._get_kmers(sequence):
            for i in self.index[kmer]:
                barcode_counts[i] += 1

        result = []
        for i, count in barcode_counts.items():
            if count < min_kmers:
                continue
            result.append((self.barcode_list[i], count))

        result = sorted(result, reverse=True, key=lambda x: x[1])
        return result[:max_hits]


def main():
    barcodes = set(["AAACAAGTATCTCCCA", "AAACAATCTACTAGCA", "GTACTCCGTGCTTCTC", "CTCGGAGCAAGTTAAG", "TACTCATTCACAAACC",
                    "CTAACTTTCACGATGT", "TCTCATACACAGATTC", "TACTCATCAGCTGTTA", "AAGGAGCAGCGTGAAC", "TCCCGATAGCGCTCCA",
                    "TGGCCAGCAATGACCT", "ACCCACTAGATATGCA", "ACCCACTGTCTAGTGT", "TGCACCTTCCCAGGTG", "AAGCCGCAGTGTCTCA",
                    "CAGCAGCCAAGAAGAG", "GTCTTCGTCCTTTCGG", "ACGAGGATCCTGCTTG", "ATCCGAAGTTCTGTTT", "ATTGGACTCAACCAAC",
                    "CATGACACATTTGCCC", "ACGGCCAGTCAGCTAT", "ATCACGATCGTACGGC", "TTAGGCAAGTTAGGTA", "CTCGTCAGTTGTTTGG",
                    "CACAGGCAGGGCTCTC", "AGGCCACGTCTAGCGC", "TTAGGACAGAGTACAT", "ATCTACTGTAGAGGAA", "TGGGCGTAGTTGTAGA",
                    "AGAGCGACACTCGACG", "GATCTAGCATTGCGGC", "GTACGTAGTGAAATCA", "CCACGGAGTGCTGTAT", "AGCGGTCGTCTTGTCC",
                    "GACAGAGCAGAGTGTG", "GCACTCTCAAGACGTG", "GCTGCGAAGCTAAGAT", "GACGTTAGTCCAAGTT", "CAAGTTGCACAGGAGT",
                    "CTCAGAAAGTTACGGG", "CTAACTTCACCGCTAG", "CACCAGGCAGTTCCCT", "GCTTCCAAGTGTCCAT", "AGATCTGGTAGGGACT",
                    "CTCGTACGTCCGTGAC", "CCCATACTCAGGCGAA", "ACGGGCTTCGCACTCT", "GCTGCGATCTACCTGC", "GCTGGGTGTCCGTTAA",
                    "ACAGCTATCCGTAGGC", "TAGTGGTGTTCCCGAG", "TCCACACGTCTCTCTG", "TTTATGCGTATTAGCC", "AGGGTGATCCGTCATC",
                    "TGCGGGTCACGTCAGC", "ACTGATGAGATCCCGC"])
    barcode_indexer = KmerIndexer(barcodes, kmer_size=5)
    # known barcode
    print(barcode_indexer.get_occurrences("AAACATCTGCTAGCACTATA"))


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)