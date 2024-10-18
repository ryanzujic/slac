import os
import sys
import unittest

# Get the local directory
test_dir = os.path.dirname(os.path.realpath(__file__))
top_dir = os.path.dirname(test_dir)

if top_dir not in sys.path:
    sys.path.append(top_dir)

from source import slac


class TestSLAC(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_metrics_all_match_genomic_only(self):
        aligned_genomic = ("A" * 5) + ("G" * 10) + ("T" * 5)
        aligned_hit = aligned_genomic
        aln = slac.SlacView(genomic=aligned_genomic, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(100), "identity for full hit alignment is not 100%")
        self.assertEqual(aln.coverage_to_genomic, float(100), "coverage for full hit alignment is not 100%")
        self.assertEqual(aln.identity_to_cds, float(0), "identity for full hit alignment to genomic only is not 0%")
        self.assertEqual(aln.coverage_to_cds, float(0), "coverage for full hit alignment to genomic only is not 0%")

    def test_metrics_all_match_genomic_and_cds(self):
        aligned_genomic = ("A" * 5) + ("G" * 10) + ("T" * 5)
        aligned_cds = aligned_genomic
        aligned_hit = aligned_genomic
        aln = slac.SlacView(genomic=aligned_genomic, cds=aligned_cds, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(100), "identity for full hit alignment is not 100%")
        self.assertEqual(aln.coverage_to_genomic, float(100), "coverage for full hit alignment is not 100%")
        self.assertEqual(aln.identity_to_cds, float(100), "identity for full hit alignment to cds is not 100%")
        self.assertEqual(aln.coverage_to_cds, float(100), "coverage for full hit alignment to cds is not 100%")

        # Make CDS only partially cover genomic with UTRs either side
        aligned_genomic = ("A" * 5) + ("C" * 10) + ("T" * 5)
        aligned_cds = ("-" * 5) + ("C" * 10) + ("-" * 5)
        aligned_hit = aligned_cds
        aln = slac.SlacView(genomic=aligned_genomic, cds=aligned_cds, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(100), "Incorrect identity to genomic")
        self.assertEqual(aln.coverage_to_genomic, float(50), "Incorrect coverage to genomic")
        self.assertEqual(aln.identity_to_cds, float(100), "Incorrect identity to cds")
        self.assertEqual(aln.coverage_to_cds, float(100), "Incorrect coverage to cds")

    def test_metrics_all_match_sample_gene(self):
        genomic = "ACTAACTGATTGGGACATGAGCTTGACCATTGGTTCACTAAATTTCCTTGATCTTTCTTCTTTGTGATGAACAAACCTATGTTTTACACATGTAAGAATGAGAAACACCTATCTCTATTTTAAACTCGGTTCTCCCTTCTCAACTTTCTGGTGTGTTAACATCAAGGGGTGTCTCAAATCTTGATCTTTATCATCCATCCTTCGCTCAGAGAAGCTTCTAAAGAGGAATAGTATGTTCGTTGAAGATAGTGGGAAAAATCATGTCAATTTCAAGAACTCAAATGGGGTGTCTCTTGGTGCTCATACTAGTATTATTAAGCAGTTCTTGTTTAAACCATGGAGCATTAGGTGCTAGATATGGGAAAAAAGTGGAGCATACACACAAGTTTAGAGTATGGTGGAAAAGCTTTAGCCATGAAACAGCTCGAGGCTCTTATAAGACAGTTAATAGAGAAGTACCTTCCTCTCCTGACCCCTTACACAACAGGTAGGTAGTTATTTCAATTGACTGTATCCTCTCTCCTTACAGAGCTTCTCGTGGGATGAGTGACAAGTTTTCTCGATCACAATTCAAGTTTTAGAAGAACTATGCCATTGTTTACTAAGTTGTTATCCTATTTTTATTTTTCTTGTGTATAGAGTTTTCTGCATATATACAATTACTTTGTATTTCATTATCGAACACTCCTAAAGCCTCCAAAAAAGAACACTCCCCAAGTACGGCATGTACAGCAATTTTCAGATTCTTAACCTAATGAAGCTGAAAACAAGGTTACTCAAGTTACCTGCCTGGAGAACCTTTGGAAAGTGGAAGTGTTCAAATTCTGTGTGTTGATGGAAATTAACTACAATGGCATGCATTTTATATTGTTACTAATTAGGCTTTGCAATACGCTGAAACATTGCTGCTCGAAAATTACAATACCATATGCTGTTTCCTAGTCTAGTAGTGATGGATGATTAATCGTCTTCCATATAGTCATTTGAAACCATCACTTTCCTTACGATTCTC"
        cds =     "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ATGTCAATTTCAAGAACTCAAATGGGGTGTCTCTTGGTGCTCATACTAGTATTATTAAGCAGTTCTTGTTTAAACCATGGAGCATTAGGTGCTAGATATGGGAAAAAAGTGGAGCATACACACAAGTTTAGAGTATGGTGGAAAAGCTTTAGCCATGAAACAGCTCGAGGCTCTTATAAGACAGTTAATAGAGAAGTACCTTCCTCTCCTGACCCCTTACACAACAGGTAG-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        hit =     "ACTAATTGGTTGGGACATGAGCTTGACCATTGGTTCACTAAATTTCCTTGATCTTTCTTCTTAGTGATGAACAGACCTATGTTTTACACATGTAAGAATGAGAAACAACTATCTCTATTTTAAACTC--AACATTATTCTCAATTGTCTGGTGTGTTAACATCAAGGGGTGTCTCAGATCTTAAT-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        aln = slac.SlacView(genomic=genomic, cds=cds, hit=hit)
        # These were informed by deciding the logic was correct at the time and locking in an expected output.
        self.assertEqual(aln.identity_to_genomic, float(90.811), "Incorrect identity to genomic")
        self.assertEqual(aln.coverage_to_genomic, float(18.083), "Incorrect coverage to genomic")

    def test_metrics_no_hit(self):
        aligned_genomic = ("A" * 5) + ("C" * 10) + ("T" * 5)
        aligned_cds = ("-" * 5) + ("C" * 10) + ("-" * 5)
        aligned_hit = "-" * 20
        aln = slac.SlacView(genomic=aligned_genomic, cds=aligned_cds, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(0), "Incorrect identity to genomic")
        self.assertEqual(aln.coverage_to_genomic, float(0), "Incorrect coverage to genomic")
        self.assertEqual(aln.identity_to_cds, float(0), "Incorrect identity to cds")
        self.assertEqual(aln.coverage_to_cds, float(0), "Incorrect coverage to cds")

        # Repeat with no cds
        aligned_cds = "_" * 20
        aln = slac.SlacView(genomic=aligned_genomic, cds=aligned_cds, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(0), "Incorrect identity to genomic")
        self.assertEqual(aln.coverage_to_genomic, float(0), "Incorrect coverage to genomic")
        self.assertEqual(aln.identity_to_cds, float(0), "Incorrect identity to cds")
        self.assertEqual(aln.coverage_to_cds, float(0), "Incorrect coverage to cds")

    def test_hit_genomic_partial_only(self):
        aligned_genomic = ("A" * 5) + ("C" * 10) + ("T" * 5)
        aligned_cds = ("-" * 5) + ("C" * 10) + ("-" * 5)
        aligned_hit = ("A" * 5) + ("-" * 15)
        aln = slac.SlacView(genomic=aligned_genomic, cds=aligned_cds, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(100), "Incorrect identity to genomic")
        self.assertEqual(aln.coverage_to_genomic, round(float(5/20)*100, 3), "Incorrect coverage to genomic")

    def test_hit_genomic_cds_partial(self):
        aligned_genomic = ("A" * 5) + ("C" * 10) + ("T" * 5)
        aligned_cds = ("-" * 5) + ("C" * 10) + ("-" * 5)
        aligned_hit = "A" + "-" * 19
        aln = slac.SlacView(genomic=aligned_genomic, cds=aligned_cds, hit=aligned_hit)
        self.assertEqual(aln.identity_to_genomic, float(100), "Incorrect identity to genomic")
        self.assertEqual(aln.coverage_to_genomic, round(float(1/20)*100, 3), "Incorrect coverage to genomic")


if __name__ == '__main__':
    unittest.main()
#%%
