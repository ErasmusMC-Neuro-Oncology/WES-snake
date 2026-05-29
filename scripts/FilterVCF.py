#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FilterVCF.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Make changes to VCF filter status.
# The Mutect2 orientation bias filter is quite strict.
# This has lead to filtering of IDH mutations before.
# Here I rescue all IDH hotspot mutations, also change their read stats
# 
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
"""
python3 scripts/PythonScript.py \
        -i {input} \
        -o {output}
"""
#
# TODO:
# 1) 
#
# History:
#  15-04-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse
import gzip
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 FilterVCF.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Modify vcf orientation bias filter   ')
    parser.add_argument('-i', help='path to input file',
                        dest='input',
                        type=str)
    parser.add_argument('-o', help='path to output file',
                        dest='output',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

"""
args.input = '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/sarek/SG_011/annotation/mutect2/SG_011_tumor1/SG_011_tumor1.annotated.vcf.gz'
args.output = 'test.vcf'
"""

#-------------------------------------------------------------------------------
# 1.1 Modify vcf
#-------------------------------------------------------------------------------
MIN_ALT = 11
MIN_AF = 0.1

with gzip.open(args.input, "rt") as infile, open(args.output, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue
        cols = line.strip().split("\t")
        # Detect IDH1 codon 132 hotspot mutations
        is_idh1_r132 = (
            "IDH1" in line and
            (
                "p.Arg132" in line or
                "p.R132" in line
            )
        )
        if is_idh1_r132:
            # Always rescue
            cols[6] = "PASS"

            # Parse FORMAT/sample fields
            format_fields = cols[8].split(":")
            sample_fields = cols[9].split(":")
            format_dict = dict(zip(format_fields, sample_fields))
            ref_reads, alt_reads = map(int, format_dict["AD"].split(","))
            dp = int(format_dict["DP"])
            af = float(format_dict["AF"])
            # Determine scaling factor needed
            scale_factor = max(
                1,
                (10 + ref_reads - 1) // ref_reads if ref_reads > 0 else 1,
                (11 + alt_reads - 1) // alt_reads if alt_reads > 0 else 1
            )
            # Scale counts if necessary
            if scale_factor > 1:
                ref_reads *= scale_factor
                alt_reads *= scale_factor
                dp = ref_reads + alt_reads
                af = alt_reads / dp
                # Always PASS
                cols[6] = "PASS"
                # Update AD
                format_dict["AD"] = f"{ref_reads},{alt_reads}"
                # Update FAD if present
            if "FAD" in format_dict:
                fad_ref, fad_alt = map(int, format_dict["FAD"].split(","))
                fad_ref *= scale_factor
                fad_alt *= scale_factor
                format_dict["FAD"] = f"{fad_ref},{fad_alt}"
            # Update DP/AF
            format_dict["DP"] = str(dp)
            format_dict["AF"] = f"{af:.3f}"
            # Rebuild sample column
            cols[9] = ":".join(format_dict[k] for k in format_fields)
        outfile.write("\t".join(cols) + "\n")
