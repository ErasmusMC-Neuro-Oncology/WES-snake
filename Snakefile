configfile: "config.yaml"
from datetime import datetime
import pandas as pd
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]
# Fetch sample wildcards
Patients = pd.read_csv(config['all']['samplesheet'])['patient'].to_numpy()
Samples = pd.read_csv(config['all']['samplesheet'])['sample'].to_numpy()
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand(output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz", patient = Patients)

#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Run Sarek tumor-only variant calling
rule Sarek:
    input:
        config['all']['samplesheet']
    output:
        samplesheet = output_dir + 'sarek/{patient}/csv/samplesheet.csv',
        bam = output_dir + "sarek/{patient}/preprocessing/mapped/{patient}_tumor1/{patient}_tumor1.sorted.bam",
        vcf = output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz"
    threads: 2
    resources:
        mem_mb=100000
    conda:
        "envs/nextflow.yaml"
    log:
        "logs/sarek/{patient}/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".log"
    params:
        genome = 'GATK.GRCh38',
        profile = "singularity",
        tools = "mutect2,merge",
        reference = config['sarek']['reference'],
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        HMF_PON = config['sarek']['HMF_PON'],
        workdir = config['sarek']['workdir'],
        outdir=lambda wildcards: f"{output_dir}sarek/{wildcards.patient}",
    shell:
        """
        run_dir={params.workdir}/sessions/{wildcards.patient}
        mkdir -p $run_dir
        cd $run_dir

        # Create samplesheet
        awk -F',' '$1=="patient" || $1=="{wildcards.patient}"' {input} > {output.samplesheet}

        nextflow -log {log} run nf-core/sarek -r 3.8.1 \
           -profile {params.profile} \
           -work-dir {params.workdir} \
           -resume \
              --input {output.samplesheet} \
              --outdir {params.outdir} \
              --genome {params.genome} \
              --tools {params.tools} \
              --intervals {params.targets} \
              --interval_padding {params.intervals} \
              --max_cpus {threads} \
              --fastp_max_cpus {threads} \
              --bwa_max_cpus {threads} \
              --gatk_max_cpus {threads} \
              --max_memory '{resources.mem_mb} MB' \
              --vep \
              --bcftools_annotations {params.HMF_PON} \
              --save_mapped \
              --wes \
              --mutect2_extra_args "--genotype-germline-sites true --genotype-pon-sites true"
        """
#+++++++++++++++++++++++++++++++++++++++++ 2 PERFORM CNA ANALYSIS +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Run CopywriteR, Normalize with QDNAseq and export results
"""
rule CNA_analysis:
    input:
        bam= output_dir + 'sarek/{patient}/preprocessing/recalibrated/{sample}/{sample}.recal.bam'
    output:
        output_dir + 'copywriter/{sample}/CNAprofiles/read_counts.txt'
    params:
genome = 'hg38',
outdir = output_dir + 'copywriter/'
cores = config'CopyWriteR']['cores']
        binsize = config['CopyWriteR']['binsize'],

    conda:
        "envs/copywritr.yaml"
    script:
        'scripts/CNA_analysis.R'
"""
