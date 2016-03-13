from partfinder._tiger import TigerDNA
from partfinder.alignment import Alignment
import numpy as np
import pytest

test_alignments = {
"small1" : """
5 10
spp1   actgactgaa
spp2   tctgtctgtt
spp3   agtgagtgaa
spp4   actaactaaa
spp5   acagacagaa
""",
"small2" : """
4 11
s1 GTGCTGGTGCT
s2 AGCCTGGTGCT
s3 AAGCTGGTGCT
s4 CGTCTGGTGCT
""",
"large1" : """
4 500
s1 CTTGAGGTTCAGAATGGTAATGAAGTGCTGGTGCTGGAAGTTCAGCAGCAGCTCGGCGGCGGTATCGTACGTACCATCGCCATGGGTTCTTCCGACGGTCTGCGTCGCGGTCTGGATGTAAAAGACCTCGAGCACCCGATCGAAGTCCCAGTTGGTAAAGCAACACTGGGTCGTATCATGAACGTACTGGGTCAGCCAGTAGACATGAAGGGCGACATCGGTGAAGAAGAGCGTTGGGCTATCCACCGTGAAGCACCATCCTATGAAGAGCTGTCAAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATGTGTCCGTTTGCGAAGGGCGGTAAAGTTGGTCTGTTCGGTGGTGCGGGTGTAGGTAAAACCGTAAACATGATGGAGCTTATTCGTAACATCGCGATCGAGCACTCCGGTTATTCTGTGTTTGCGGGCGTAGGTGAACGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACCGA
s2 CTTGAGGTACAAAATGGTAATGAGAGCCTGGTGCTGGAAGTTCAGCAGCAGCTCGGTGGTGGTATCGTACGTGCTATCGCCATGGGTTCTTCCGACGGTCTGCGTCGTGGTCTGGAAGTTAAAGACCTTGAGCACCCGATCGAAGTCCCGGTTGGTAAAGCAACGCTGGGTCGTATCATGAACGTGCTGGGTCAGCCGATCGATATGAAAGGCGACATCGGCGAAGAAGAACGTTGGGCGATTCACCGTGCAGCACCTTCCTATGAAGAGCTCTCCAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATGTGTCCGTTCGCGAAGGGCGGTAAAGTCGGTCTGTTCGGTGGTGCGGGTGTTGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGATCGAACACTCCGGTTACTCCGTGTTTGCTGGTGTTGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACCGA
s3 CTTGAGGTACAGAATAACAGCGAGAAGCTGGTGCTGGAAGTTCAGCAGCAGCTCGGCGGCGGTATCGTACGTACCATCGCAATGGGTTCTTCCGACGGTCTGCGTCGTGGTCTGGAAGTGAAAGACCTCGAGCACCCGATCGAAGTCCCGGTAGGTAAAGCGACCCTGGGTCGTATCATGAACGTGCTGGGTCAGCCAATCGATATGAAAGGCGACATCGGCGAAGAAGATCGTTGGGCGATTCACCGCGCAGCACCTTCCTATGAAGAGCTGTCCAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATTTGTCCGTTCGCTAAGGGCGGTAAAGTTGGTCTGTTCGGTGGTGCGGGCGTAGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGATCGAGCACTCCGGTTACTCCGTGTTTGCAGGCGTGGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAGATGACCGA
s4 CTCGAGGTGAAAAATGGTGATGCTCGTCTGGTGCTGGAAGTTCAGCAGCAGCTGGGTGGTGGCGTGGTTCGTACCATCGCCATGGGTACTTCTGACGGCCTGAAGCGCGGTCTGGAAGTTACCGACCTGAAAAAACCTATCCAGGTTCCGGTTGGTAAAGCAACCCTCGGCCGTATCATGAACGTATTGGGTGAGCCAATCGACATGAAAGGCGACCTGCAGAATGACGACGGCACTGTAATTCACCGTGCAGCACCTTCGTATGAAGATCAGTCTAACTCGCAGGAACTGCTGGAAACCGGCATCAAGGTTATCGACCTGATGTGTCCGTTCGCTAAGGGCGGTAAAGTCGGTCTGTTCGGTGGTGCGGGTGTAGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGGCTGAGCACTCAGGTTATTCGGTATTTGCTGGTGTGGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACTGA
"""
}

@pytest.fixture(scope="module", params=test_alignments.keys())
def alignment_text(request):
    return test_alignments[request.param] 

def test_tigger(alignment_text):
    a = Alignment()
    a.parse(alignment_text)
    bitsets, slow_rates = slow_tiger(a)
    tigger = TigerDNA()
    tigger.build_bitsets(a)
    fast_rates = tigger.calc_rates()
    assert np.allclose(slow_rates, fast_rates)
    

def slow_tiger(alm):
    denom = float(alm.sequence_length - 1)
    # Initiate a list to store the bitsets for each set partition
    set_parts = []
    # Loop through the sites in each aligment to generate bitsets
    for site in alm.data.T:
        sites = []
        # Generate a bitset for each nucleotide. First see if there are any N's,
        # if there are, you must create bitsets for each nucleotide. Note: Make
        # sure that all ambiguous characters have been converted to N's prior to
        # running, i.e. all '-' should be 'N'
        N = (site == ord('N'))
        for nuc in map(ord, list('ACGT')):
            # Add the bitset ONLY if it is non-empty
            sb = (site == nuc)
            if np.any(sb) or np.any(N):
                sites.append(sb)
        if np.any(N):
            for s in sites:
                s |= N
        set_parts.append(sites)

    # Initiate a list to store the rates
    rates = []
    # Look at each set partition in the alignment
    for i, sp_i in enumerate(set_parts):
        rate = 0.0
        # Now compare it against every other set partition in the alignment,
        # except for itself
        for j, sp_j in enumerate(set_parts):
            # Do not compare the site against itself
            if j == i:
                continue

            # Keep track of the number of comparisons done
            num = 0.0
            # 'axpi' refers to the compatibility of one bitset to another. If
            # the bitset is a subset of another, the axpi is 1, otherwise it
            # is 0
            axpi = 0.0
            # Check each partition set from site pattern j to see if it is a
            # site pattern i
            for ps in sp_j:
                for ps2 in sp_i:
                    if np.array_equal((ps & ps2), ps):
                        axpi += 1.0
                        break
                num += 1.0
            rate += axpi / num

        rate /= denom
        rates.append(rate)
    return set_parts, rates
