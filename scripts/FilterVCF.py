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

rescued_header = '##INFO=<ID=RESCUED_IDH,Number=0,Type=Flag,Description="IDH hotspot variant rescued from orientation bias filter">\n'

with gzip.open(args.input, "rt") as infile, open(args.output, "w") as outfile:
    header_written = False
    for line in infile:
        if line.startswith("#"):
            # Insert our new INFO definition before the #CHROM line
            if line.startswith("#CHROM") and not header_written:
                outfile.write(rescued_header)
                header_written = True
            outfile.write(line)
            continue
        
        cols = line.strip().split("\t")
        is_idh1_r132 = (
            "IDH1" in line and
            ("p.Arg132" in line or "p.R132" in line)
        )
        if is_idh1_r132:
            cols[6] = "PASS"
            
            
            # ... rest of your scaling logic unchanged
            format_fields = cols[8].split(":")
            sample_fields = cols[9].split(":")
            format_dict = dict(zip(format_fields, sample_fields))
            ref_reads, alt_reads = map(int, format_dict["AD"].split(","))
            dp = int(format_dict["DP"])
            af = float(format_dict["AF"])
            scale_factor = max(
                1,
                (10 + ref_reads - 1) // ref_reads if ref_reads > 0 else 1,
                (11 + alt_reads - 1) // alt_reads if alt_reads > 0 else 1
            )
            if scale_factor > 1:
                ref_reads *= scale_factor
                alt_reads *= scale_factor
                dp = ref_reads + alt_reads
                af = alt_reads / dp
                format_dict["AD"] = f"{ref_reads},{alt_reads}"
                if "FAD" in format_dict:
                    fad_ref, fad_alt = map(int, format_dict["FAD"].split(","))
                    fad_ref *= scale_factor
                    fad_alt *= scale_factor
                    format_dict["FAD"] = f"{fad_ref},{fad_alt}"
                format_dict["DP"] = str(dp)
                format_dict["AF"] = f"{af:.3f}"
                cols[7] = cols[7] + ";RESCUED_IDH" if cols[7] != "." else "RESCUED_IDH"
            cols[9] = ":".join(format_dict[k] for k in format_fields)
        outfile.write("\t".join(cols) + "\n")
