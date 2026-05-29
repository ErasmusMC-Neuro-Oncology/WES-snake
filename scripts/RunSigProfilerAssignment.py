#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RunSigProfilerAssignment.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Merge per-sample trinucleotide count files into a single SBS96 matrix
# and run SigProfilerAssignment cosmic_fit.
# Designed for WES data; uses exome-renormalized COSMIC v3.5 signatures.
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: sigprofiler
# Usage:
"""
python3 scripts/RunSigProfilerAssignment.py \
        --input {input.trinuc} \
        --output {params.out_dir} \
        --genome {params.genome}
"""
#
# TODO:
# 1)
#
# History:
#  29-05-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse
import os
import pandas as pd
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 RunSigProfilerAssignment.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Merge per-sample trinucleotide count files and run SigProfilerAssignment  ')
    parser.add_argument('--input', help='space-separated list of per-sample trinucleotide count files',
                        dest='input',
                        type=str,
                        nargs='+')
    parser.add_argument('--output', help='path to output folder',
                        dest='output',
                        type=str)
    parser.add_argument('--genome', help='reference genome build (e.g. GRCh37, GRCh38)',
                        dest='genome',
                        type=str,
                        default='GRCh38')
    parser.add_argument('--cosmic-version', help='COSMIC signature version (default: 3.5)',
                        dest='cosmic_version',
                        type=float,
                        default=3.3)
    args = parser.parse_args()
    return args
args = parse_args()

#-------------------------------------------------------------------------------
# 1.1 Read and merge trinucleotide count files
#-------------------------------------------------------------------------------
dfs = []
for filepath in args.input:
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    dfs.append(df)

merged = pd.concat(dfs, axis=1)

assert merged.shape[0] == 96, (
    f"Expected 96 mutation types after merging, got {merged.shape[0]}. "
    "Check that all input files use the same MutationType labels."
)

print(f"Merged {merged.shape[1]} samples across {merged.shape[0]} mutation types")
#-------------------------------------------------------------------------------
# 1.2 Write merged matrix
#-------------------------------------------------------------------------------
os.makedirs(args.output, exist_ok=True)
merged_matrix_path = os.path.join(args.output, 'merged_sbs96_matrix.txt')
merged.to_csv(merged_matrix_path, sep='\t')

print(f"Merged matrix written to {merged_matrix_path}")
#-------------------------------------------------------------------------------
# 2.1 Run SigProfilerAssignment
#-------------------------------------------------------------------------------
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(
    samples        = merged_matrix_path,
    output         = args.output,
    input_type     = 'matrix',
    cosmic_version = args.cosmic_version,
    exome          = True,           # WES data: use exome-renormalized signatures
    genome_build   = args.genome,
    export_probabilities = True,
    make_plots     = True,
    verbose        = False
)

print(f"SigProfilerAssignment completed. Results written to {args.output}")
