#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FilterVCF.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Make changes to VCF filter status.
# The Mutect2 orientation bias filter is quite strict.
# This has lead to filtering of IDH mutations before.
# Here I rescue pathogenic or likely pathogenic variants flagged as 'orientation' 
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
args.input = '../output/sarek/I18-1048-01/annotation/mutect2/I18-1048-01_tumor1/I18-1048-01_tumor1.annotated.vcf.gz'
args.output = 'test.vcf'
"""

#-------------------------------------------------------------------------------
# 1.1 Modify vcf
#-------------------------------------------------------------------------------
# Open files
with gzip.open(args.input, "rt") as infile, open(args.output, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue
        cols = line.strip().split("\t")
        if cols[6] == "orientation":
            if 'pathogenic' in line:
                cols[6] = "PASS"
        outfile.write("\t".join(cols) + "\n")
