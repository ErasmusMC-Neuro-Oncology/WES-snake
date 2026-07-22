#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RunSigProfilerAssignment.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run SigProfilerAssignment with COSMIC v3.3
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
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

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 RunSigProfilerAssignment.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Run SigProfilerAssignment on a folder of per-sample somatic VCF files  ')
    parser.add_argument('--input', help='path to folder containing per-sample somatic VCF files',
                        dest='input',
                        type=str)
    parser.add_argument('--output', help='path to output folder',
                        dest='output',
                        type=str)
    parser.add_argument('--genome', help='reference genome build (e.g. GRCh37, GRCh38)',
                        dest='genome',
                        type=str,
                        default='GRCh38')
    parser.add_argument('--cosmic-version', help='COSMIC signature version (default: 3.3)',
                        dest='cosmic_version',
                        type=float,
                        default=3.3)
    args = parser.parse_args()
    return args
args = parse_args()
"""
args.input = 'output/WES/SigProfilerAssignment/vcf/'
args.output = 'output/WES/SigProfilerAssignment/'
args.genome = 'GRCh38'
args.cosmic_version = 3.3
"""
#-------------------------------------------------------------------------------
# 1.1 Validate input directory
#-------------------------------------------------------------------------------
assert os.path.isdir(args.input), (
    f"Input path '{args.input}' is not a directory. "
    "Please provide a folder containing per-sample VCF files."
)

vcf_files = [f for f in os.listdir(args.input) if f.endswith('.vcf')]
assert len(vcf_files) > 0, (
    f"No .vcf files found in '{args.input}'. "
    "Check that FilterSomaticVCF completed successfully."
)

print(f"Found {len(vcf_files)} VCF files in {args.input}")
os.makedirs(args.output, exist_ok=True)

#-------------------------------------------------------------------------------
# 2.1 Install reference genome if not already present
#-------------------------------------------------------------------------------
from SigProfilerMatrixGenerator import install as genInstall

try:
    genInstall.install(args.genome, rsync=False, bash=True)
    print(f"Reference genome {args.genome} installed successfully")
except Exception as e:
    print(f"Genome installation note: {e}")
    print("Continuing — genome may already be installed")

#-------------------------------------------------------------------------------
# 2.2 Run SigProfilerAssignment
#-------------------------------------------------------------------------------
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(
    samples = args.input,
    output = args.output,
    input_type = 'vcf',
    cosmic_version = args.cosmic_version,
    exome = True,
    genome_build = args.genome,
    exclude_signature_subgroups = [
        'Immunosuppressants_signatures',   # drops SBS32, SBS87 — no plausible azathioprine/thiopurine exposure in glioma
        'UV_signatures',                   # not CNS-relevant
        'AA_signatures',                   # aristolochic acid — not CNS-relevant
        'Colibactin_signatures',           # gut-specific
        'Lymphoid_signatures',             # AID/RAG-related, not relevant outside lymphoid malignancies
    ],
    export_probabilities = True,
    export_probabilities_per_mutation = True,
    make_plots = True,
    verbose = False
)

print(f"SigProfilerAssignment completed. Results written to {args.output}")
