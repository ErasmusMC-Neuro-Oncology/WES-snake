configfile: "config.yaml"
from datetime import datetime
import pandas as pd
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]

# Fetch Patient wildcards
Patients = pd.read_csv(config['all']['samplesheet'])['patient'].to_numpy()
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand(output_dir + 'QDNAseq/{binsize}/{patient}/plots/QDNAseq_segmented_profile.pdf', patient = Patients, binsize = config['CopyWriteR']['binsizes'])

#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Run Sarek variant calling
rule Sarek:
    input:
        config['all']['samplesheet']
    output:
        samplesheet = output_dir + 'sarek/{patient}/csv/samplesheet.csv',
        recal = temp(directory(output_dir + "sarek/{patient}/preprocessing/recalibrated/")),
        md = temp(directory(output_dir + "sarek/{patient}/preprocessing/markduplicates/")),
        sorted_cram = temp(output_dir + "sarek/{patient}/preprocessing/mapped/{patient}_tumor1/{patient}_tumor1.sorted.cram"),
        sorted_bam = output_dir + "sarek/{patient}/preprocessing/mapped/{patient}_tumor1/{patient}_tumor1.sorted.bam",
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
        workdir=lambda wildcards: f"{config['sarek']['workdir']}/{wildcards.patient}",
        outdir=lambda wildcards: f"{output_dir}/sarek/{wildcards.patient}",
    shell:
        """
        export NXF_WORK={params.workdir}
        
        # Subset samplesheet
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

        # Save alignment as .bam (to be fixed with --save-output-as-bam in new sarek release)
        samtools view -b -o {output.sorted_bam} {output.sorted_cram}
        
        # Clean cache and intermediate files upon completion but keep on failure
        status=$?
        if [ $status -eq 0 ]; then
        rm -rf {params.workdir}
        else
        echo "Sarek failed"
        fi
        
        """


#+++++++++++++++++++++++++++++++++++++++++ 2 PERFORM CNA ANALYSIS +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Run CopywriteR, Normalize with QDNAseq and export results
rule CNA_analysis:
    input:
        bam= output_dir + "sarek/{patient}/preprocessing/mapped/{patient}_tumor1/{patient}_tumor1.sorted.bam"
    output:
        sample_dir = temp(directory(output_dir + "copywriter/{binsize}/{patient}/")),
        QDNAseq = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segmented.Rds',
        Segments = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segments.txt',
        Profile = output_dir + 'QDNAseq/{binsize}/{patient}/plots/QDNAseq_segmented_profile.pdf',
    params:
        genome = 'hg38',
        cores = config['CopyWriteR']['cores'],
        outdir=lambda wildcards: f"{output_dir}/copywriter/{wildcards.binsize}/",
    conda:
        "envs/copywriter.yaml"
    script:
        'scripts/CNA_analysis.R'
