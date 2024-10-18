import math
from collections import Counter

# Converts comparison symbols to alternatives that can be used to convey the quality/frequency of that symbol within a
# simplified block of symbols.
SYMBOL_QUALITY_MAP = {"-": ["-", "~"],
                      "_": ["_", "-"],
                      "|": ["|", "!"],
                      "O": ["O", "o"],
                      "^": ["^", "`"],
                      "X": ["X", "x"],
                      ".": [".", ","],
                      "?": ["?", "?"],
                      "=": ["=", ":"],
                      "U": ["U", "u"],
                      }

DEFAULT_SIZE = 50
# What proportion of the block must be the most common character for it to be used in the short version of the alignment
# and not the alternative character.
DEFAULT_FREQ_THRESHOLD = 1.0

# Number of decimal places to round the metrics to
METRIC_DECIMALS = 3

# Default string display modes for the object.
DISPLAY_MODE_SHORT = "short"
DISPLAY_MODE_FULL = "full"
DISPLAY_MODE_FULL_HIT = "hits"

class SLAC(object):
    def __init__(self, genomic=None, cds=None, hit=None, size_limit=DEFAULT_SIZE,
                 frequency_threshold=DEFAULT_FREQ_THRESHOLD, display_mode=None, auto_align_cds_to_genomic=False):
        """
        A class designed to represent a hit sequence aligned to the query gene's genomic and coding sequence, within
        a single continuous text string. This includes an abbreviated format called miniSLAC, which can fit within an
        near arbitrary character length, at the expense of resolution.

        Args:
            genomic: (str) the genomic sequence, aligned to the hit sequence.
            cds: (str) the coding sequence. Ideally this will be supplied aligned to the other sequences, however an
            unaligned CDS can be provided and aligned to the genomic sequence if required - and a simple algorithm
            will attempt to align it to the genomic sequence.
            hit: (str) the hit sequence, aligned to the genomic sequence.
            size_limit: (int) the maximum number of characters to display in the miniSLAC output. Default is 50.
            frequency_threshold: (float) the proportion of each abbreviated block upon which the miniSLAC is built,
            which must be the most common character, in order for the full version of that character to be used instead
            of the mixed character. Default is 1.0.
            Eg: The 10 character block of |||X|||||| will be abbreviated to the mixed character '!' if the frequency
            threshold is 1.0. If it was 0.8, it would be abbreviated to the full character '|'. You'll probably want it
            on the default to not miss details.
            display_mode: (str) the default display mode for the object. Must be one of: short, full, full_hit.
            Default is short (miniSLAC).
            auto_align_cds_to_genomic: (bool) if True, will attempt to align the CDS to the genomic sequence if the CDS
            is provided but not aligned. Default is False. if False, a provided CDS must be pre-aligned to the genomic.

        Raises:
            ValueError: Conditions are not met for the input sequences.
        """
        if genomic and not isinstance(genomic, str):
            raise ValueError(f"genomic must be a string, not {type(genomic)}")
        if cds and not isinstance(cds, str):
            raise ValueError(f"cds must be a string, not {type(cds)}")
        if hit and not isinstance(hit, str):
            raise ValueError(f"hit must be a string, not {type(hit)}")

        if not any([genomic, cds, hit]):
            raise ValueError("No sequences supplied")

        if cds and not genomic:
            raise ValueError("If a cds sequence is supplied, a genomic sequence must also be supplied. "
                             "Supply cds as genomic if only cds is available as the cds is only a layer of context.")

        if display_mode and display_mode not in [DISPLAY_MODE_SHORT, DISPLAY_MODE_FULL, DISPLAY_MODE_FULL_HIT]:
            raise ValueError(f"Invalid display mode: {display_mode}. Must be one of: "
                             f"{DISPLAY_MODE_SHORT}, {DISPLAY_MODE_FULL}, {DISPLAY_MODE_FULL_HIT}")

        if not display_mode:
            display_mode = DISPLAY_MODE_SHORT

        if genomic and hit:
            if len(genomic) != len(hit):
                raise ValueError("Genomic and hit sequences must be aligned, and therefore be the same length")

        max_length = max([len(genomic), len(cds), len(hit)])

        self.genomic = genomic.upper().replace(" ", "-") if genomic else ""
        self.cds = cds.upper().replace(" ", "-") if cds else ""
        self.hit = hit.upper().replace(" ", "-")  if hit else ""
        
        # Perform auto alignment of cds to genomic if requested and required.
        if auto_align_cds_to_genomic and all([self.genomic, self.hit]) and "-" not in self.cds: 
            cds = self.align_cds_to_genomic(genomic, hit)

        # Fill in blanks to allow use in cases where a sequence type was not provided
        if not self.genomic:
            self.genomic = "-" * max_length
        if not cds:
            self.cds = "-" * max_length
        if not hit:
            self.hit = "-" * max_length

        self.size_limit = int(size_limit)
        # The full length, symbolic alignment representation
        self._full = ''
        # The full length, symbolic alignment representation, characters that match the CDS in the hit use an uppercase
        # DNA character, and those that only match non-coding regions use a lowercase DNA character.
        self._full_hit = ''
        # The shortened version of the symbolic alignment
        self._short = ''
        self.display_mode = display_mode
        self._frequency_threshold = frequency_threshold
        self.identity_to_genomic = None
        self.identity_to_cds = None
        self.coverage_to_genomic = None
        self.coverage_to_cds = None
        self.concordance_to_genomic = None
        self.concordance_to_cds = None
        self.hit_span = None
        self._generate_text()

    def __str__(self):
        if self.display_mode == DISPLAY_MODE_FULL:
            return self._full
        elif self.display_mode == DISPLAY_MODE_FULL_HIT:
            return self._full_hit
        else:
            return self._short

    def __repr__(self):
        return self._short

    def short(self):
        return self._short

    def full(self, encoded_hit=False):
        """
        Returns the full length, symbolic alignment representation.
        Args:
            encoded_hit: (bool) if True, returns the full length symbolic alignment representation with the matching
            positions using the hit characters.
        Returns:
            (str) the full length, symbolic alignment representation.
        """
        if encoded_hit:
            return self._full_hit
        return self._full

    def set_size(self, size_limit):
        self.size_limit = size_limit
        self._generate_text()

    def _generate_text(self):
        """
        Generate the text representation of the alignment.
        Requires:
            Hit and query sequences are the same length
        Returns:
            None
        """

        self._generate_full()
        self._generate_short()
        self._generate_full_hit()

    def _generate_full(self):
        """
        Assumes genomic and cds match in all positions where cds is not a gap.
        """

        """
        if all match |
        if only g is present _
        if only h is present ^
        if only g and c are present =
        if only c is present ? (unexpected, genomic should be present in g)
        if h does not match g X
        if g and h match but c is absent - O
        else U (unknown case to fix)
        TODO - Update this to match the actual logic below
        """

        coding_matches = 0
        coding_mismatches = 0
        non_coding_matches = 0
        non_coding_mismatches = 0
        coding_inserts = 0
        non_coding_inserts = 0
        coding_deletions = 0
        non_coding_deletions = 0
        uncounted = 0
        unmatched_hit_overhang = 0
        coding_overhang = 0  # Number of positions of the coding sequence outside the hit span
        non_coding_overhang = 0  # Number of positions of the non-coding sequence outside the hit span

        i = 0
        in_exon = False
        in_genomic = False
        in_hit = False
        for g, c, h in zip(self.genomic, self.cds, self.hit):

            # --- Identify larger structures within the sequences ---

            # We need to know if we're in an exon to correctly measure insertions relative to cds
            if c != "-":
                in_exon = True
            # Kick out of an exon if we can see we've left it based on the genomic sequence.
            # If we just said "if c == '-'" we could be inside an insertion and not know if it's a coding region or not.
            elif c == "-" and g != "-":
                in_exon = False

            # Track if we're in an overhanging part of the hit
            if g != "-":
                in_genomic = True
            elif self.genomic[i:len(self.genomic)] == "-" * (len(self.genomic) - i):
                in_genomic = False
                # If we've left the genomic seq we've definitely left any exons
                in_exon = False

            # Determine if we're within the hit span
            if not in_hit and h != "-":
                in_hit = True
            elif in_hit and self.hit[i:len(self.hit)] == "-" * (len(self.hit) - i):
                in_hit = False

            # --- Determine character-by-character cases relative to these structures ---

            # Match.
            if g != "-" and g == c == h:
                self._full += "|"
                coding_matches += 1
            # Only in genomic.
            elif g != "-" and c == "-" and h == "-":
                self._full += "_"
                if in_hit:
                    non_coding_deletions += 1
                else:
                    non_coding_overhang += 1
            # Only in hit.
            elif g == "-" and c == "-" and h != "-":
                self._full += "^"
                if in_exon:
                    coding_inserts += 1
                elif in_genomic:
                    non_coding_inserts += 1
                else:
                    unmatched_hit_overhang += 1
            # Only in genomic and cds.
            elif g != "-" and c != "-" and h == "-":
                self._full += "="
                if in_hit:
                    coding_deletions += 1
                else:
                    coding_overhang += 1
            # Only in cds (not expected).
            elif c != "-" and g == "-" and h == "-":
                self._full += "?"
                uncounted += 1
            # Mismatch between coding and genomic (not expected).
            elif c != "-" and g != c:
                self._full += "?"
                uncounted += 1
            # Mismatch in a region covered by the cds.
            elif h != g and c != "-":
                self._full += "X"
                coding_mismatches += 1
            # Mismatch in a region not covered by cds (ie, in a non coding region)
            elif h != g and c == "-":
                self._full += "."
                non_coding_mismatches += 1
            # Match in a region not covered by cds.
            elif g != "-" and h == g and c == "-":
                self._full += "O"
                non_coding_matches += 1
            # All gaps, unexpected.
            elif g == "-" and g == c == h:
                raise ValueError(f"All sequences have gaps at position {i}. There may be issues with the alignment.")
            # Unplanned logical error, make obvious
            else:
                self._full += "U"
                raise ValueError(f"Unexpected case in alignment: {g}, {c}, {h}")

            i += 1

        # --- Calculate metrics ---

        if uncounted:
            # print(f"There were {uncounted} uncounted positions in the alignment, applying None to metrics to avoid "
            #       f"misleading results")
            self.identity_to_genomic = None
            self.identity_to_cds = None
            self.coverage_to_genomic = None
            self.coverage_to_cds = None

        if coding_matches or non_coding_matches:
            # This represents how much of the hit sequence is the same as the genomic sequence
            self.identity_to_genomic = (float(coding_matches) + non_coding_matches) / (coding_matches +
                                                                                       coding_mismatches +
                                                                                       non_coding_matches +
                                                                                       non_coding_mismatches +
                                                                                       coding_inserts +
                                                                                       non_coding_inserts +
                                                                                       coding_deletions +
                                                                                       non_coding_deletions)

            self.identity_to_genomic = round(self.identity_to_genomic * 100, METRIC_DECIMALS)

            # This represents how much of the genomic sequence mapped to matches or mismatches. This is not the span of
            # the total alignment.
            self.coverage_to_genomic = (float(coding_matches) +
                                        coding_mismatches +
                                        non_coding_matches +
                                        non_coding_mismatches) / (coding_matches +
                                                                  coding_mismatches +
                                                                  non_coding_matches +
                                                                  non_coding_mismatches +
                                                                  coding_deletions +
                                                                  non_coding_deletions +
                                                                  non_coding_overhang +
                                                                  coding_overhang)

            self.coverage_to_genomic = round(self.coverage_to_genomic * 100, METRIC_DECIMALS)

            self.concordance_to_genomic = round((self.coverage_to_genomic * self.identity_to_genomic) / 100,
                                                METRIC_DECIMALS)
        else:
            self.identity_to_genomic = 0
            self.coverage_to_genomic = 0
            self.concordance_to_genomic = 0

        if coding_matches:
            # This represents how much of the hit sequence is the same as the cds sequence
            self.identity_to_cds = float(coding_matches) / (coding_matches +
                                                            coding_mismatches +
                                                            coding_inserts +
                                                            coding_deletions)

            self.identity_to_cds = round(self.identity_to_cds * 100, METRIC_DECIMALS)

            # This represents how much of the cds sequence mapped to matches or mismatches. This is not the span of the
            # alignment to the cds.
            self.coverage_to_cds = (float(coding_matches) + coding_mismatches) / (coding_matches +
                                                                                  coding_mismatches +
                                                                                  coding_deletions +
                                                                                  coding_overhang)

            self.coverage_to_cds = round(self.coverage_to_cds * 100, METRIC_DECIMALS)

            self.concordance_to_cds = round((self.coverage_to_cds * self.identity_to_cds) / 100, METRIC_DECIMALS)

        else:
            self.identity_to_cds = 0
            self.coverage_to_cds = 0
            self.concordance_to_cds = 0

    def _generate_full_hit(self):
        if not self._full:
            self._generate_full()

        for i, char in enumerate(self._full):
            if char == "|":
                self._full_hit += self.hit[i].upper()
            elif char == "O":
                self._full_hit += self.hit[i].lower()
            else:
                self._full_hit += self._full[i]

    def _generate_short(self):
        block_size = int(math.floor(float(len(self.genomic)) / self.size_limit))
        if block_size < 1 or len(self.genomic) <= self.size_limit:
            self._short = self._full
            return

        self._short = ""

        # Generate shortened version
        i = 0
        while i < len(self._full):

            # If our next block would be the last one...
            if len(self._short) + 1 == self.size_limit:
                # Capture all the remaining characters and increment our counter.
                block_size = len(self._full) - i

            # If we're about to run out of sequence without filling the size limit...
            elif i + block_size >= len(self._full) and len(self._short) + 1 < self.size_limit:
                # Shorten the final blocks to squeeze an extra one at the end.
                block_size = max(1, block_size - 1)

            block = self._full[i:i + block_size]

            if len(block) > 1:
                # Find the most common character in the block
                counter = Counter(block)
                common_char, count = counter.most_common(1)[0]
                block_char = common_char

                # Retain the character if it's preserved enough across this block
                if count >= (float(block_size) * self._frequency_threshold):
                    self._short += block_char

                # If the block only contains matches, whether they be to the cds + genomic, or only to the genomic
                # sequence (ie UTR/intron), use the most frequent character.
                elif all([char in ["|", "O"] for char in block]):
                    self._short += block_char

                # If the block only contains parts of the query and has no coverage of the hit, use the most frequent
                # character.
                elif all([char in ["=", "_"] for char in block]):
                    self._short += block_char

                else:
                    self._short += SYMBOL_QUALITY_MAP[block_char][1]

            else:
                self._short += block

            i += block_size

    @staticmethod
    def align_cds_to_genomic(aligned_genomic, cds, kmer_size=5):
        """
        As we can assume the coding sequence will align to the genomic sequence with perfect identity, with gaps,
        the user can provide an un-aligned cds, we can perform simple kmer matching to determine how the cds
        aligns with the genomic, and retain full functionality of the tool.

        Disclaimer: This is not a full alignment tool and has not yet been tested in edge cases, particularly if
        insertions are common.

        Args:
            aligned_genomic: (str) the genomic sequence as aligned to the hit, if one was provided.
            cds: (str) the coding sequence
            kmer_size: (int) the size of the kmers to use for matching. We check kmer windows between the sequences and
            insert gaps in the CDS until each kmer matches the genomic sequence. When a subsequent kmer does not match,
            we retry with one size down until we reach 0, at which point we insert gaps until the next full kmer matches.
            This assumes that any exon is at least the size of the kmer. Default is 5.

        Assumes:
            The genomic sequences contains the full and exact CDS sequence.
            The CDS sequence does not contain gaps.
            No exons are smaller than the kmer size.

        Returns:
            (str) the CDS sequence with gaps inserted to align it to the genomic sequence

        Raises:
            ValueError: If the CDS contains gaps or if the entire CDS cannot be aligned by the end of the genomic sequence.
        """

        aligned_genomic = aligned_genomic.replace(" ", "-")
        cds = cds.replace(" ", "-")

        # Confirm cds has no gaps, as this would suggest it's already aligned.
        if "-" in cds:
            raise ValueError("CDS should not contain gaps for auto alignment to genomic sequence")

        aligned_cds = []
        cds_index = 0

        i = 0
        while i < len(aligned_genomic):

            # A gap in the genomic is always a gap in the CDS
            if aligned_genomic[i] == "-":
                aligned_cds.append('-')
                i += 1
                continue

            match_found = False

            # Try k-mer sizes from the provided kmer_size down to 1
            for k in range(kmer_size, 0, -1):
                genomic_kmer = aligned_genomic[i:i + k]

                # Ensure the k-mer from the CDS is not longer than the remaining CDS sequence
                if cds_index + k <= len(cds):
                    cds_kmer = cds[cds_index:cds_index + k]
                else:
                    cds_kmer = cds[cds_index:]

                # If the k-mers match, append them and move the index forward
                if genomic_kmer == cds_kmer:
                    aligned_cds.append(cds_kmer)
                    cds_index += len(cds_kmer)
                    i += len(genomic_kmer)
                    match_found = True
                    break  # Exit the k-mer size reduction loop since we found a match

            if not match_found:
                # If no match is found, insert a gap in the CDS
                aligned_cds.append('-')
                i += 1  # Move to the next base in the genomic sequence

        # Failure detection: check if we have unaligned characters left in the CDS
        if cds_index < len(cds):
            raise ValueError(f"Failed CDS to genomic alignment: {len(cds) - cds_index} bases of the CDS were not aligned to the genomic sequence. "
                             f"Check the full and exact CDS is within the genomic sequence.")

        return ''.join(aligned_cds)

