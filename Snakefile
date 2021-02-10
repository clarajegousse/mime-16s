# cp results/cutadapt/*/*.fastq.gz trimmed/

# Make sure that you set the `truncLen=` option in the rule `dada2_filter_and_trim_pe` according
# to the results of the quality profile checks (after rule `dada2_quality_profile_pe` has finished on all samples).
# If in doubt, check https://benjjneb.github.io/dada2/tutorial.html#inspect-read-quality-profiles

SAMPLES = ["FX003-016-16S-V4_S58","FX003-017-16S-V4_S59", "FX003-018-16S-V4_S60"]
RUNS = ["20190508_0074"]

rule all:
    input:
        # In a first run of this meta-wrapper, comment out all other inputs and only keep this one.
        # Looking at the resulting plot, adjust the `truncLen` in rule `dada2_filter_trim_pe` and then
        # rerun with all inputs uncommented.
        expand("results/reports/dada2/filter-trim-pe/{run}/{sample}.tsv",
        sample = SAMPLES,
        run = RUNS),
        #"results/dada2/taxa.RDS"

rule cutadapt:
    input:
        fwd = "data/miseq/{run}/{sample}_L001_R1_001.fastq.gz",
        rev = "data/miseq/{run}/{sample}_L001_R2_001.fastq.gz",
    output:
        fwd = "results/trimmed/{run}/{sample}.1.fastq.gz",
        rev = "results/trimmed/{run}/{sample}.2.fastq.gz",
        report = "results/reports/cutadapt/{run}/{sample}-qc-report.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        # https://earthmicrobiome.org/protocols-and-standards/16s/
        # Updated sequences: 515F (Parada)–806R (Apprill), forward-barcoded:
        # FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
        # echo 'GGACTACNVGGGTWTCTAAT' | rev
        adapter_a = "^GTGYCAGCMGCCGCGGTAA...AATGGCGCCGMCGACYGTG",
        adapter_A = "^GGACTACNVGGGTWTCTAAT...TAATCTWTGGGVNCATCAGG",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        minimum_length = 150,
        maximum_length = 280
        #quality_cutoff = 20
    log:
        "logs/cutadapt/{run}/{sample}.log"
    shell:
        "cutadapt -a {params.adapter_a} -A {params.adapter_A} \
         -m {params.minimum_length} -M {params.maximum_length} \
         --discard-untrimmed \
         -o {output.fwd} -p {output.rev} \
          {input.fwd} {input.rev} \
          2> {output.report}"

rule dada2_quality_profile_pe:
    input:
        # FASTQ file without primer sequences
        expand("results/trimmed/{{run}}/{{sample}}.{orientation}.fastq.gz", orientation = [1,2])
    output:
        "results/reports/dada2/quality-profile/{run}/{sample}-quality-profile.png"
    log:
        "logs/dada2/quality-profile/{run}/{sample}-quality-profile-pe.log"
    wrapper:
        "0.70.0/bio/dada2/quality-profile"

rule dada2_filter_trim_pe:
    input:
        # Paired-end files without primer sequences
        fwd="results/trimmed/{run}/{sample}.1.fastq.gz",
        rev="results/trimmed/{run}/{sample}.2.fastq.gz"
    output:
        filt="results/filtered-pe/{run}/{sample}.1.fastq.gz",
        filt_rev="results/filtered-pe/{run}/{sample}.2.fastq.gz",
        stats="results/reports/dada2/filter-trim-pe/{run}/{sample}.tsv"
    params:
        # Set the maximum expected errors tolerated in filtered reads
        maxEE=[2,2],
        # Set the number of kept bases in forward and reverse reads
        truncLen=[240,200]
    log:
        "logs/dada2/filter-trim-pe/{run}/{sample}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/filter-trim"
