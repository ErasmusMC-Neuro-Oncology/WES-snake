configfile: "config.yaml"
from datetime import datetime
import pandas as pd
import os

#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
output_dir = config["all"]["output_dir"]

# Fetch Patient wildcards
Patients = pd.read_csv(config['all']['samplesheet'])['patient'].unique()
Patients = [pat for pat in Patients if not pat in ['SG_004','SG_040']]


#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        output_dir + 'plots/Barplot_target_depth.pdf',
        expand(output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_variants.csv', patient = Patients, binsize = config['CopyWriteR']['binsizes'])

        
#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Download Sarek
rule Download_Sarek:
    output:
        sarek = directory(config['sarek']['workdir'] + ".nf-core-sarek/"),
        singularity_dir = directory(f"{config['sarek']['workdir']}singularity/cache/")
    conda:
        "envs/nextflow.yaml"
    shell:
        """
	export NXF_SINGULARITY_CACHEDIR={output.singularity_dir}        
	nf-core pipelines download --outdir {output.sarek} --container-system singularity --force --compress none -r dev sarek
        """

# 1.2 Run Sarek
rule Sarek:
    input:
        samplesheet = output_dir + 'samplesheets/Samplesheet_WES_{patient}.csv',
        sarek = config['sarek']['workdir'] + ".nf-core-sarek/"
    output:
        mapped = temp(directory(output_dir + "sarek/{patient}/preprocessing/mapped/")),
        recal = temp(directory(output_dir + "sarek/{patient}/preprocessing/recalibrated/")),
        md_cram = temp(output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.cram"),
        md_bam = output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.bam",
        vcf = temp(output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.mutect2.filtered_snpEff_VEP.ann.vcf.gz"),
        vcf_annotated = output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz",
        depth = output_dir + "sarek/{patient}/reports/mosdepth/{patient}_tumor1/{patient}_tumor1.md.mosdepth.summary.txt"
    threads: 8
    resources:
        mem_mb=100000,
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
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        HMF_PON = config['sarek']['HMF_PON'],
        COSMIC = config['sarek']['COSMIC'],
        Mutect2_params = config['sarek']['Mutect2_params'],
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/",
        workdir=lambda wildcards: f"{config['sarek']['workdir']}/{wildcards.patient}",
        outdir=lambda wildcards: f"{output_dir}/sarek/{wildcards.patient}",
    shell:
        """
        export NXF_WORK={params.workdir}
        export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}
	export NXF_CACHE_DIR={params.workdir}

        nextflow -log {log} run {input.sarek}/{params.version}/ \
           -profile {params.profile} \
           -work-dir {params.workdir} \
           -resume \
           -c {params.Mutect2_params} \
              --input {input.samplesheet} \
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

        # Add HMF PON and COSMIC annotation
        bcftools annotate {output.vcf} -a {params.HMF_PON} -c INFO -Oz -o {params.workdir}/tmp.vcf.gz ; bcftools index -t {params.workdir}/tmp.vcf.gz
        bcftools annotate {params.workdir}/tmp.vcf.gz -a {params.COSMIC} -c INFO -Oz -o {output.vcf_annotated}; bcftools index -t {output.vcf_annotated}
        
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
        mem_mb=100000,
        gpu=0,
        runtime='30h'
        
    params:
        genome = 'hg38',
        cores = config['CopyWriteR']['cores'],
        cytobands = config['CopyWriteR']['cytobands'],
        ACE_purity_penalty = config['ACE']['penalty'],
        ACE_ploidy_penalty = config['ACE']['penploidy'],
        CNH_path = config['CNH']['path'],
        PON = lambda wildcards: config['CNA_PON'][wildcards.binsize]
    conda:
        "envs/CNA.yaml"
    script:
        'scripts/CNA_analysis.R'


#+++++++++++++++++++++++++++++++++++++++++ 3 DISTINGUISH GERMLINE-SOMATIC +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Run PureCN to call tumor purity/ploidy, classify variants and calculate CCF
rule PureCN:
    input:
        vcf =  output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz",
        Segments = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segments.txt'
    output:
        vcf = temp(output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/PureCN_{binsize}.vcf"),
        intervals = temp(output_dir + 'PureCN/{binsize}/{patient}/baits_hg19_intervals.txt'),
        PureCN_rds = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1.rds',
        PureCN_purity = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1.csv',
        variants = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_variants.csv',
        TMB = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_mutation_burden.csv',
        signatures = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_signatures.csv'
    params:
        genome = 'hg38',
        ref = config['all']['ref'],
        snake_dir = config['all']['snake_dir'],
        targets = config['sarek']['targetregions'],
        min_af = config['PureCN']['min_af'],
        min_alt = config['PureCN']['min_alt'],
        min_bq = config['PureCN']['min_bq'],
        min_cosmic_cnt = config['PureCN']['min_cosmic_cnt'],
        CNA_sdev =  config['PureCN']['CNA_sdev'],
        CNA_max_nonclonal =  config['PureCN']['CNA_max_nonclonal'],
        outdir=lambda wildcards: f"{output_dir}/PureCN/{wildcards.binsize}/{wildcards.patient}",
    conda:
        "envs/PureCN.yaml"
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

        # Modify input vcf
        python3 {params.snake_dir}/scripts/FilterVCF.py -i {input.vcf} -o {output.vcf}
        
        # Run PureCN
        Rscript {params.snake_dir}/scripts/PureCN.R \
        --out {params.outdir} \
        --sampleid {wildcards.patient}_tumor1 \
        --segfile {input.Segments} \
        --vcf {output.vcf} \
        --intervals {output.intervals} \
        --genome {params.genome} \
        --segsdev {params.CNA_sdev} \
        --min-af {params.min_af} \
        --min-base-quality {params.min_bq} \
        --min-supporting-reads {params.min_alt} \
        --min-cosmic-cnt {params.min_cosmic_cnt} \
        --max-non-clonal {params.CNA_max_nonclonal} \
        --cosmic-cnt-info-field GENOME_SCREEN_SAMPLE_COUNT
        
        # Calculate signatures/statistics
         Rscript {params.snake_dir}/scripts/Dx.R \
        --rds {output.PureCN_rds} \
        --callable {params.targets} \
        --signatures \
        --force
        """


#++++++++++++++++++++++++++++++++++++++++++++++++ 4 MERGE SAMPLE DATA +++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create samplesheet with stats
rule Create_SampleData:
    input:
        depth = expand(output_dir + "sarek/{patient}/reports/mosdepth/{patient}_tumor1/{patient}_tumor1.md.mosdepth.summary.txt",patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        CNA_stats = expand(output_dir + 'QDNAseq/{binsize}/{patient}/data/CNA_stats.txt',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        ACE_results = expand(output_dir + 'ACE/{binsize}/{patient}/ACE_fits.txt',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        CNH_results = expand(output_dir + 'CNH/{binsize}/{patient}/CNH_results.txt',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        PureCN_purity = expand(output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1.csv',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        TMB = expand(output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_mutation_burden.csv',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        variants = expand(output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_variants.csv',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
        signatures = expand(output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_signatures.csv',patient = Patients, binsize = config['CopyWriteR']['binsizes'])
    output:
        SampleData = output_dir + 'sampledata/SampleData_WES.txt'
    conda:
        "envs/R.yaml"
    script:
        'scripts/Create_SampleData.R'

#++++++++++++++++++++++++++++++++++++++++++++++++ 5 PLOT SAMPLE DATA +++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot SampleData
rule Plot_SampleData:
    input:
        SampleData = output_dir + 'sampledata/SampleData_WES.txt'
    output:
        Barplot_depth = output_dir + 'plots/Barplot_target_depth.pdf'
    conda:
        "envs/R.yaml"
    script:
        'scripts/Plot_SampleData.R'
        
        
