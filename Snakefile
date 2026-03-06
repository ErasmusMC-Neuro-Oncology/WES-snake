configfile: "config.yaml"
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
        "samplesheet.csv"
    output:
        "results/sarek/.done"
    threads: 2
    resources:
        mem_mb=10000
    conda:
        "envs/nextflow.yaml"
    params:
        genome="GRCh38",
        profile="singularity",
        tools="mutect2",
        outdir = output_dir + 'sarek'
    shell:
        """
        nextflow run nf-core/sarek -r 3.8.1 \
            -profile {params.profile} \
            --input {input} \
            --outdir results/sarek \
            --genome {params.genome} \
            --tools {params.tools} \
            -resume \
            -max_cpus {threads} \
            -max_memory '{resources.mem_mb} MB'

        touch {output}
        """

#-------------------------------------------------------------------------------------------------------------------
# 1.2
