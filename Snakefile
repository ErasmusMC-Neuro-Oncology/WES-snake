configfile: "config.yaml"
from datetime import datetime
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        "results/sarek/.done"

#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Run Sarek
rule Sarek:
    input:
        "../MINT/data/samplesheets/samplesheet_WES.csv"
    output:
        "results/sarek/.done"
    threads: 2
    resources:
        mem_mb=10000
    conda:
        "envs/nextflow.yaml"
    log:
        "logs/sarek/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".log"
    params:
        genome="hg38",
        profile="singularity",
        tools="mutect2,strelka,merge",
        outdir = output_dir + 'sarek',
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        gnomAD = config['sarek']['gnomAD'],
        dbSNP = config['sarek']['dbSNP']
    shell:
        """
        nextflow -log {log} run nf-core/sarek -r 3.8.1 \
            -profile {params.profile} \
            --input {input} \
            --outdir {params.outdir} \
            --genome {params.genome} \
            --tools {params.tools} \
            --intervals {params.targets} \
            --interval_padding {params.intervals} \
            --germline_resource {params.gnomAD} \
            --dbsnp {params.dbSNP} \
            --max_cpus {threads} \
            --max_memory '{resources.mem_mb} MB'
            --somatic \
            --germline \
            --wes \
            -resume 
        
        touch {output}
        """

