configfile: "config.yaml"
from datetime import datetime
import pandas as pd
import os
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
output_dir = config["all"]["output_dir"]

# Fetch Patient wildcards
Patients = pd.read_csv('samplesheet.csv')['patient'].unique() if os.path.isfile('samplesheet.csv') else []

#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        'samplesheet.csv',
        expand(output_dir + 'QDNAseq/{binsize}/{patient}/plots/QDNAseq_segmented_profile.pdf', patient = Patients, binsize = config['CopyWriteR']['binsizes'])

#++++++++++++++++++++++++++++++++++++++++++++ 0 CREATE SAMPLESHEET ++++++++++++++++++++++++++++++++++++++++++++++++
rule Create_Samplesheet:
    params:
        data_dir = config['all']['data_dir'],
        sample_overview = '../data/MINT_db.xlsx',
    output:
        'samplesheet.csv'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_Samplesheet_MINT.R'

        
#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Download Sarek
rule Download_Sarek:
    output:
        ".nf-core-sarek/"
    params:
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/"
    conda:
        "envs/nextflow.yaml"
    shell:
        """
	export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}        
	nf-core pipelines download --outdir {output} --container-system singularity --compress none -r dev sarek
        """

# 1.2 Run
rule Sarek:
    input:
        samplesheet = 'samplesheet.csv',
        sarek = ".nf-core-sarek/"
    output:
        samplesheet = output_dir + 'sarek/{patient}/csv/samplesheet.csv',
        mapped = temp(directory(output_dir + "sarek/{patient}/preprocessing/mapped/")),
        recal = temp(directory(output_dir + "sarek/{patient}/preprocessing/recalibrated/")),
        md_cram = temp(output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.cram"),
        md_bam = temp(output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.bam"),
        vcf = temp(output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz"),
        vcf_annotated = output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz"
    threads: 8
    resources:
        mem_mb=50000,
        gpu=0,
        runtime='30h'
    conda:
        "envs/nextflow.yaml"
    log:
        "logs/sarek/{patient}/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".log"
    params:
        version = 'dev',
        genome = 'GATK.GRCh38',
        profile = "singularity",
        tools = "mutect2,merge",
        Mutect2_params = 'params/mutect2_params.json',
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        HMF_PON = config['sarek']['HMF_PON'],
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/",
        workdir=lambda wildcards: f"{config['sarek']['workdir']}/{wildcards.patient}",
        outdir=lambda wildcards: f"{output_dir}/sarek/{wildcards.patient}",
    shell:
        """
        export NXF_WORK={params.workdir}
        export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}

        # Subset samplesheet
        awk -F',' '$1=="patient" || $1=="{wildcards.patient}"' {input.samplesheet} > {output.samplesheet}

        nextflow -log {log} run .nf-core-sarek/{params.version}/ \
           -profile {params.profile} \
           -work-dir {params.workdir} \
           -c {params.Mutect2_params} \
              --input {output.samplesheet} \
              --outdir {params.outdir} \
              --genome {params.genome} \
              --tools {params.tools} \
              --intervals {params.targets} \
              --interval_padding {params.intervals} \
              --max_memory '{resources.mem_mb} MB' \
              --save_mapped  \
              --wes

        # Save alignment as .bam (to be fixed with --save-output-as-bam in new sarek release)
        samtools view -b -o {output.md_bam} {output.md_cram}

        # Add HMF PON annotation
        bcftools annotate {output.vcf} -a {params.HMF_PON} -c INFO -O z -o {output.vcf_annotated}
        bcftools index -t {output.vcf_annotated}
        
        # Clean cache and intermediate files upon completion but keep on failure
        status=$?
        if [ $status -eq 0 ]; then
        rm -rf {params.workdir}
        else
        echo "Sarek failed"
        fi
        
        """


#+++++++++++++++++++++++++++++++++++++++++ 2 PERFORM CNA ANALYSIS +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Run CopywriteR, QDNAseq, ACE, CNH, calculate stats and export results
rule CNA_analysis:
    input:
        bam= output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.bam"
    output:
        sample_dir = temp(directory(output_dir + "copywriter/{binsize}/{patient}/")),
        QDNAseq = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segmented.Rds',
        Profile = output_dir + 'QDNAseq/{binsize}/{patient}/plots/QDNAseq_segmented_profile.pdf',
        Segments = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segments.txt',
        Segments_igv = temp(output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segments.igv'),
        Called = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_calls.txt',
        CNA_stats = output_dir + 'QDNAseq/{binsize}/{patient}/data/CNA_stats.txt',
        ACE_results = output_dir + 'ACE/{binsize}/{patient}/ACE_fits.txt',
        ACE_matrix = output_dir + 'ACE/{binsize}/{patient}/ACE_matrixplot.pdf',
        CNH_results = output_dir + 'CNH/{binsize}/{patient}/CNH_results.txt',
        CNH_plot = output_dir + 'CNH/{binsize}/{patient}/CNH_plot.pdf',
        CNH_error_plot = output_dir + 'CNH/{binsize}/{patient}/CNH_errorplot.pdf',
    resources:
        mem_mb=10000,
        gpu=0,
        runtime='30h'
    params:
        genome = 'hg38',
        cores = config['CopyWriteR']['cores'],
        cytobands = config['CopyWriteR']['cytobands'],
        ACE_purity_penalty = config['ACE']['penalty'],
        ACE_ploidy_penalty = config['ACE']['penploidy'],
        CNH_path = config['CNH']['path']
    conda:
        "envs/CNA.yaml"
    script:
        'scripts/CNA_analysis.R'


#+++++++++++++++++++++++++++++++++++++++++ 3 DISTINGUISH GERMLINE-SOMATIC +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Run PureCN to call tumor purity/ploidy, classify variants and calculate CCF  
rule PureCN:
    input:
        vcf =  output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz",
        Segments = output_dir + 'QDNAseq/100kbp/{patient}/data/QDNAseq_Segments.txt'
    output:
        intervals = temp(output_dir + 'PureCN/{patient}/baits_hg19_intervals.txt'),
        PureCN_rds = output_dir + 'PureCN/{patient}/{patient}_tumor1.rds',
        variants = output_dir + 'PureCN/{patient}/{patient}_tumor1_variants.csv'
    params:
        genome = 'hg38',
        outdir = output_dir + 'PureCN/{patient}/',
        ref = config['PureCN']['ref'],
        targets = config['sarek']['targetregions']
    conda:
        "envs/purecn.yaml"
    shell:
        """
        # Find PureCN installation
        PureCN_lib=$CONDA_PREFIX/lib/R/library/PureCN/extdata

        # Create intervals file
        Rscript $PureCN_lib/IntervalFile.R \
        --in-file {params.targets} \
        --fasta {params.ref} \
        --out-file {output.intervals} \
        --off-target \
        --genome {params.genome}
        
        # Run PureCN
        Rscript $PureCN_lib/PureCN.R \
        --out {params.output_dir} \
        --sampleid {patient} \
        --segfile {input.Segments} \
        --vcf {input.vcf} \
        --intervals {output.intervals} \
        --genome {params.genome}

        # Calculate signatures/statistics
         Rscript $PureCN_lib/Dx.R \
        --rds {output.PureCN_rds} \
        --callable {params.targets} \
        --signatures
        """
