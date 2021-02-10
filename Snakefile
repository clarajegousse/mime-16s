# Make sure that you set the `truncLen=` option in the rule `dada2_filter_and_trim_pe` according
# to the results of the quality profile checks (after rule `dada2_quality_profile_pe` has finished on all samples).
# If in doubt, check https://benjjneb.github.io/dada2/tutorial.html#inspect-read-quality-profiles

rule all:
    input:
        # In a first run of this meta-wrapper, comment out all other inputs and only keep this one.
        # Looking at the resulting plot, adjust the `truncLen` in rule `dada2_filter_trim_pe` and then
        # rerun with all inputs uncommented.
        expand(
            "results/dada2/reports/quality-profile/{sample}-quality-profile.png",
            sample=["FX003-016-16S-V4_S58","FX003-017-16S-V4_S59"]
        ),
        "results/dada2/taxa.RDS"

rule dada2_quality_profile_pe:
    input:
        # FASTQ file without primer sequences
        expand("results/cutadapt/20190508_0074/{{sample}}_L001_{orientation}_001.fastq.gz",orientation=["R1","R2"])
    output:
        "results/dada2/reports/quality-profile/20190508_0074/{sample}-quality-profile.png"
    log:
        "logs/dada2/quality-profile/20190508_0074/{sample}-quality-profile-pe.log"
    wrapper:
        "0.70.0/bio/dada2/quality-profile"

rule dada2_filter_trim_pe:
    input:
        # Paired-end files without primer sequences
        fwd = "results/cutadapt/20190508_0074/{sample}_L001_R1_001.fastq.gz",
        rev = "results/cutadapt/20190508_0074/{sample}_L001_R2_001.fastq.gz"
    output:
        filt = "results/dada2/20190508_0074/filtered-pe/{sample}_L001_R1_001.fastq.gz",
        filt_rev = "results/dada2/20190508_0074/filtered-pe/{sample}_L001_R2_001.fastq.gz",
        stats="results/dada2/20190508_0074/filter-trim-pe/{sample}.tsv"
    params:
        # Set the maximum expected errors tolerated in filtered reads
        maxEE=1,
        # Set the number of kept bases in forward and reverse reads
        truncLen=[240,200]
    log:
        "logs/dada2/filter-trim-pe/{sample}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/filter-trim"
