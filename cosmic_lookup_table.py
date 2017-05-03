"""
Parses CosmicMutantExport TSV to create a lookup table that connects the following
data fields: gene, HGVS genomic and protein change, chrom, start, end, ref, and alt.
"""

import sys
import re
import pandas

class CosmicLookup(object):
    """
    Class that uses lookup table to return relevant information.
    """

    def __init__(self, lookup_table_file):
        self.lookup_table = pandas.read_csv(lookup_table_file, sep="\t")

    def get_entry(self, gene, hgvs_p):
        """
        Returns a dataframe of results from filtering on gene and hgvs_p
        """
        lt = self.lookup_table
        hgvs_p = "p." + hgvs_p
        return lt[(lt['gene'] == gene) & (lt['hgvs_p'] == hgvs_p)]

#
# Methods below parse CosmicMutantExport TSV and create the lookup table.
#

def parse_hgvc_c(hgvs_c):
    """
    Parses HGVS like c.4146T>A and returns a dictionary with the keys pos, type, ref, alt
    """

    if not hgvs_c.startswith("c.") or "?" in hgvs_c:
        return {}

    parts = re.split(r"([[_\.ACGT]+|[0-9]+|del|ins])", hgvs_c[2:])

    bases = "ATCG"
    pos = ctype = ref = alt = None
    if parts[4] == ">":
        # Substitution.
        pos = parts[1]
        ctype = "sub"
        ref = parts[3]
        alt = parts[5]

    return {
        "pos": pos,
        "type": ctype,
        "ref": ref,
        "alt": alt
    }

def parse_genome_pos(genome_pos):
    """
    Parse genome position in format chrom:start-end and returns the tuple chrom, start, end
    """

    if not genome_pos:
        return None, None, None
    chrom, pos = genome_pos.split(":")
    start, end = pos.split("-")
    return chrom, start, end


def print_lookup_table(input_stream):
    """
    Print COSMIC lookup table by reading from input_stream and outputing corresponding line.
    """
    print "\t".join(["gene", "hgvs_c", "hgvs_p", "build", "chrom", "start", "end", "ref", "alt", "strand"])
    for line in input_stream:
        fields = line.split("\t")
        gene = fields[0]
        hgvs_c, hgvs_p = fields[17:19]
        build, genome_pos, strand = fields[22:25]

        parse_results = parse_hgvc_c(hgvs_c)

        chrom, start, end = parse_genome_pos(genome_pos)
        if bool(parse_results) and parse_results["ref"] and chrom:
            print "\t".join([gene, hgvs_c, hgvs_p, build, chrom, start, end, parse_results["ref"], parse_results["alt"], strand])
        else:
            # TODO: Debugging.
            #print >> sys.stderr, "LINE IGNORED", hgvs_c
            pass


if __name__ == "__main__":
    print_lookup_table(sys.stdin)
