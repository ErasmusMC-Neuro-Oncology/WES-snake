configfile: "config.yaml"
from datetime import datetime
import pandas as pd
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]

# Fetch Patient wildcards
Patients = pd.read_csv(config['all']['samplesheet'])['patient'].to_numpy()
Patients = ['MINT03']
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand(output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz", patient = Patients)
        #expand(output_dir + 'QDNAseq/{binsize}/{patient}/plots/QDNAseq_segmented_profile.pdf', patient = Patients, binsize = config['CopyWriteR']['binsizes'])

#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Run Sarek variant calling
rule Sarek:
    input:
        config['all']['samplesheet']
    output:
        samplesheet = output_dir + 'sarek/{patient}/csv/samplesheet.csv',
        mapped = temp(directory(output_dir + "sarek/{patient}/preprocessing/mapped/")),
        md = temp(directory(output_dir + "sarek/{patient}/preprocessing/markduplicates/")),
        recal_cram = temp(output_dir + "sarek/{patient}/preprocessing/recalibrated/{patient}_tumor1/{patient}_tumor1.recal.cram"),
        recal_bam = output_dir + "sarek/{patient}/preprocessing/recalibrated/{patient}_tumor1/{patient}_tumor1.recal.bam",
        vcf = output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz"
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
        genome = 'GATK.GRCh38',
        profile = "singularity",
        tools = "mutect2,merge",
        Mutect2_params = 'params/mutect2_params.json',
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

        nextflow -log {log} run nf-core/sarek -r dev \
           -profile {params.profile} \
           -work-dir {params.workdir} \
           -resume \
           -c {params.Mutect2_params} \
              --input {output.samplesheet} \
              --outdir {params.outdir} \
              --genome {params.genome} \
              --tools {params.tools} \
              --intervals {params.targets} \
              --interval_padding {params.intervals} \
              --max_memory '{resources.mem_mb} MB' \
              --bcftools_annotations {params.HMF_PON} \
              --save_mapped  \
              --wes

        # Save alignment as .bam (to be fixed with --save-output-as-bam in new sarek release)
        samtools view -b -o {output.recal_bam} {output.recal_cram}
        
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
        bam= output_dir + "sarek/{patient}/preprocessing/recalibrated/{patient}_tumor1/{patient}_tumor1.recal.bam"
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
    params:
        genome = 'hg38',
        cores = config['CopyWriteR']['cores'],
        cytobands = config['CopyWriteR']['cytobands'],
        ACE_purity_penalty = config['ACE']['penalty'],
        ACE_ploidy_penalty = config['ACE']['penploidy'],
        outdir=lambda wildcards: f"{output_dir}/copywriter/{wildcards.binsize}/",
    conda:
        "envs/CNA.yaml"
    script:
        'scripts/CNA_analysis.R'


#+++++++++++++++++++++++++++++++++++++++++ 3 DISTINGUISH GERMLINE-SOMATIC +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Run PureCN to call tumor purity/ploidy, classify variants and calculate CCF         
rule PureCN:
    input:
        vcf = output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz",        
        Segments = output_dir + 'QDNAseq/100kbp/{patient}/data/QDNAseq_Segments.txt'
    output:
        intervals = temp(output_dir + 'PureCN/{patient}/baits_hg19_intervals.txt'),
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
        """
