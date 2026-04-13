configfile: "config.yaml"
from datetime import datetime
import pandas as pd
import os
import glob
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
output_dir = config["all"]["output_dir"]
data_dir = config["all"]["data_dir"]
# Fetch Patient wildcards
Patients = pd.read_csv('samplesheet.csv')['patient'].unique() if os.path.isfile('samplesheet.csv') else []
Tumors = list(set(['_'.join(pat.split('_')[:2]) for pat in Patients]))

#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand(output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz", patient = Patients),
        expand(output_dir + 'QDNAseq/{binsize}/sWGS/{tumor}/data/QDNAseq_Segments_sWGS.txt',binsize = config['CopyWriteR']['binsizes'], tumor = Tumors),
        expand(output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segmented.Rds',binsize = config['CopyWriteR']['binsizes'], patient = Patients),

#++++++++++++++++++++++++++++++++++++++++++++ 0 CREATE SAMPLESHEET ++++++++++++++++++++++++++++++++++++++++++++++++
rule Create_Samplesheet:
    output:
        'samplesheet.csv'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_Samplesheet_GLASS.R'

        
#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Download Sarek
rule Download_Sarek:
    output:
        directory(".nf-core-sarek/")
    params:
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/"
    conda:
        "envs/nextflow.yaml"
    shell:
        """
	export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}        
	nf-core pipelines download --outdir {output} --container-system singularity --compress none -r dev sarek
        """

# 1.2 Run Sarek
rule Sarek:
    input:
        samplesheet = os.path.abspath('samplesheet.csv'),
        sarek = os.path.abspath(".nf-core-sarek/")
    output:
        samplesheet = output_dir + 'sarek/{patient}/csv/samplesheet.csv',
        mapped = temp(directory(output_dir + "sarek/{patient}/preprocessing/mapped/")),
        recal = temp(directory(output_dir + "sarek/{patient}/preprocessing/recalibrated/")),
        md_cram = temp(output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.cram"),
        md_bam = temp(output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.bam"),
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
        os.path.abspath("logs/sarek/{patient}/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".log")
    params:
        version = 'dev',
        genome = 'GATK.GRCh38',
        profile = "singularity",
        tools = "mutect2,merge",
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        HMF_PON = config['sarek']['HMF_PON'],
        Mutect2_params = os.path.abspath('params/mutect2_params.json'),
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/",
        workdir=lambda wildcards: f"{config['sarek']['workdir']}/{wildcards.patient}",
        outdir=lambda wildcards: f"{output_dir}/sarek/{wildcards.patient}",
        Mutect2_mode=lambda wildcards: "tumor1_vs_normal1" if "Paired" in wildcards.patient else "tumor1",
    shell:
        """
        export NXF_WORK={params.workdir}
        export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}
	export NXF_CACHE_DIR={params.workdir}
        
        # Subset samplesheet
        awk -F',' '$1=="patient" || $1=="{wildcards.patient}"' {input.samplesheet} > {output.samplesheet}

        nextflow -log {log} run {input.sarek}/{params.version}/ \
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
              --save_mapped \
              --wes
        
        # Save alignment as .bam (to be fixed with --save-output-as-bam in new sarek release)
        samtools view -b -o {output.md_bam} {output.md_cram}
        
         # Add HMF PON annotation
        VCF={params.outdir}/annotation/mutect2/{wildcards.patient}_{params.Mutect2_mode}/{wildcards.patient}_{params.Mutect2_mode}.mutect2.filtered_snpEff_VEP.ann.vcf.gz

        bcftools annotate $VCF -a {params.HMF_PON} -c INFO -O z -o {output.vcf_annotated}
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
# 2.1 Run QDNAseq, ACE, CNH, calculate stats and export results
rule CNA_analysis:
    input:
        bam=  output_dir + "sarek/{patient}/preprocessing/markduplicates/{patient}_tumor1/{patient}_tumor1.md.bam"
    output:
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
        mem_mb=50000,
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


rule CNA_analysis_sWGS:
    input:
        bam=  lambda wildcards: glob.glob(data_dir + 'sWGS/' + wildcards.tumor + '_*.bam')
    output:
        QDNAseq = output_dir + 'QDNAseq/{binsize}/sWGS/{tumor}/data/QDNAseq_Segmented_sWGS.Rds',
        Profile = output_dir + 'QDNAseq/{binsize}/sWGS/{tumor}/plots/QDNAseq_segmented_profile_sWGS.pdf',
        Segments = output_dir + 'QDNAseq/{binsize}/sWGS/{tumor}/data/QDNAseq_Segments_sWGS.txt',
    resources:
        mem_mb=50000,
        gpu=0,
        runtime='30h'
    conda:
        "envs/CNA.yaml"
    script:
        'scripts/CNA_analysis_sWGS.R'

#+++++++++++++++++++++++++++++++++++++++++ 3 DISTINGUISH GERMLINE-SOMATIC +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Run PureCN to call tumor purity/ploidy, classify variants and calculate CCF
rule PureCN:
    input:
        vcf =  output_dir + "sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz",
        Segments = output_dir + 'QDNAseq/{binsize}/{patient}/data/QDNAseq_Segments.txt'
    output:
        intervals = temp(output_dir + 'PureCN/{binsize}/{patient}/baits_hg19_intervals.txt'),
        PureCN_rds = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1.rds',
        variants = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_variants.csv',
        TMB = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_mutation_burden.csv',
        signatures = output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_signatures.csv'
    params:
        genome = 'hg38',
        ref = config['all']['ref'],
        targets = config['sarek']['targetregions'],
        min_af = config['PureCN']['min_af'],
        min_bq = config['PureCN']['min_bq'],        
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
        
        # Run PureCN
        Rscript $PureCN_lib/PureCN.R \
        --out {params.outdir} \
        --sampleid {wildcards.patient}_tumor1 \
        --segfile {input.Segments} \
        --vcf {input.vcf} \
        --intervals {output.intervals} \
        --genome {params.genome} \
        --min-af {params.min_af} \
        --min-base-quality {params.min_bq}
        
        # Calculate signatures/statistics
         Rscript $PureCN_lib/Dx.R \
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
        TMB = expand(output_dir + 'PureCN/{binsize}/{patient}/{patient}_tumor1_mutation_burden.csv',patient = Patients, binsize = config['CopyWriteR']['binsizes']),
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
        
        
