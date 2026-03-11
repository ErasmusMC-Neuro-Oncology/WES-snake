configfile: "config.yaml"
from datetime import datetime
import pandas as pd
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]
# Fetch sample wildcards
Samples = pd.read_csv(config['all']['samplesheet'])['sample'].to_numpy()
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand("results/sarek/{sample}.done", sample = Samples[0])

#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Run Sarek
rule Sarek:
    input:
        config['all']['samplesheet']
    output:
        "results/sarek/{sample}.done"
    threads: 2
    resources:
        mem_mb=100000
    conda:
        "envs/nextflow.yaml"
    log:
        "logs/sarek/{sample}/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".log"
    params:
        profile="singularity",
        tools="mutect2,merge",
        reference = config['sarek']['reference'],
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        gnomAD = config['sarek']['gnomAD'],
        dbSNP = config['sarek']['dbSNP'],
        COSMIC = config['sarek']['COSMIC'],
        HMF_PON = config['sarek']['HMF_PON'],
        workdir = config['sarek']['workdir'],
        outdir=lambda wildcards: f"{output_dir}sarek/{wildcards.sample}",
    shell:
        """
        nextflow -log {log} run nf-core/sarek -r 3.8.1 \
            -profile {params.profile} \
            -work-dir {params.workdir} \
            --sample {wildcards.sample} \
            --input {input} \
            --outdir {params.outdir} \
            --fasta {params.reference} \
            --tools {params.tools} \
            --intervals {params.targets} \
            --interval_padding {params.intervals} \
            --germline_resource {params.gnomAD} \
            --dbsnp {params.dbSNP} \
            --max_cpus {threads} \
            --fastp_max_cpus {threads} \
            --bwa_max_cpus {threads} \
            --gatk_max_cpus {threads} \
            --max_memory '{resources.mem_mb} MB' \
            --vep \
            --vep_custom {params.COSMIC},COSMIC,vcf,exact,0,ID \
            --vep_custom {params.HMF_PON},HMF_PON,vcf,exact,0 \
            --tumor_only \
            --wes \
            --mutect2_extra_args "--genotype-germline-sites true --genotype-pon-sites true"
        -resume
        
        touch {output}
        """

