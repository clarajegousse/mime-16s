
# snakemake --cluster qsub -j 32 -R dada2_learn_errors --latency-wait 60

# Make sure that you set the `truncLen=` option in the rule `dada2_filter_and_trim_pe` according
# to the results of the quality profile checks (after rule `dada2_quality_profile_pe` has finished on all samples).
# If in doubt, check https://benjjneb.github.io/dada2/tutorial.html#inspect-read-quality-profiles

import pandas as pd

configfile: "config.yaml"

SampleTable = pd.read_csv(config['SAMPLETABLE'], sep = ",", header = 0, index_col = 0)
SAMPLES = list(SampleTable.index)

ORIENTATION = config["ORIENTATION"]
RUN = config["RUN"]

rule all:
    input:
        # In a first run of this meta-wrapper, comment out all other inputs and only keep this one.
        # Looking at the resulting plot, adjust the `truncLen` in rule `dada2_filter_trim_pe` and then
        # rerun with all inputs uncommented.
        expand("results/reports/cutadapt/{run}/{sample}-qc-report.txt", run = RUN, sample = SAMPLES),
        expand("results/dada2/quality-profile/{run}/{sample}-quality-profile.png", run = RUN, sample = SAMPLES),
        expand("results/dada2/filtered_trim_pe/{run}/{sample}.tsv", run = RUN, sample = SAMPLES),
        # expand("reports/dada2/learn-errors/{run}/errors_{orientation}.png", run = RUN, orientation = ORIENTATION),
        # expand("results/dada2/denoised/{run}/{sample}_{orientation}.RDS", run = RUN, sample = SAMPLES, orientation = ORIENTATION),
        # expand("results/dada2/merged/{run}/{sample}.RDS", run = RUN, sample = SAMPLES),
        # expand("results/dada2/taxa/{run}/taxa.RDS", run = RUN)

rule cutadapt:
    input:
        fwd = "data/miseq/{run}/{sample}_R1.fastq.gz",
        rev = "data/miseq/{run}/{sample}_R2.fastq.gz",
    output:
        fwd = "results/trimmed/{run}/{sample}_R1.fastq.gz",
        rev = "results/trimmed/{run}/{sample}_R2.fastq.gz",
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
        expand("results/trimmed/{{run}}/{{sample}}_{orientation}.fastq.gz", orientation = ORIENTATION)
    output:
        "results/dada2/quality-profile/{run}/{sample}-quality-profile.png"
    log:
        "logs/dada2/quality-profile/{run}/{sample}-quality-profile-pe.log"
    wrapper:
        "0.70.0/bio/dada2/quality-profile"

rule dada2_filter_trim_pe:
    input:
        # Paired-end files without primer sequences
        fwd="results/trimmed/{run}/{sample}_R1.fastq.gz",
        rev="results/trimmed/{run}/{sample}_R2.fastq.gz"
    output:
        filt="results/dada2/filtered_trim_pe/{run}/{sample}_R1.fastq.gz",
        filt_rev="results/dada2/filtered_trim_pe/{run}/{sample}_R2.fastq.gz",
        stats="results/dada2/filtered_trim_pe/{run}/{sample}.tsv"
    params:
        # Set the maximum expected errors tolerated in filtered reads
        maxEE=2,
        # Set the number of kept bases in forward and reverse reads
        truncLen=[240,200]
    log:
        "logs/dada2/filter-trim-pe/{run}/{sample}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/filter-trim"

rule dada2_learn_errors:
    input:
    # Quality filtered and trimmed forward FASTQ files (potentially compressed)
        expand("results/dada2/filtered_trim_pe/{{run}}/{sample}_{{orientation}}.fastq.gz", sample = SAMPLES)
    output:
        err="results/dada2/learn-errors/{run}/model_{orientation}.RDS",# save the error model
        plot="reports/dada2/learn-errors/{run}/errors_{orientation}.png",# plot observed and estimated rates
    params:
        randomize=True
    log:
        "logs/dada2/learn-errors/{run}/learn-errors_{orientation}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/learn-errors"

rule dada2_dereplicate_fastq:
    input:
    # Quality filtered FASTQ file
        "results/dada2/filtered_trim_pe/{run}/{fastq}.fastq.gz"
    output:
    # Dereplicated sequences stored as `derep-class` object in a RDS file
        "results/dada2/uniques/{run}/{fastq}.RDS"
    log:
        "logs/dada2/dereplicate-fastq/{run}/{fastq}.log"
    wrapper:
        "0.70.0/bio/dada2/dereplicate-fastq"

rule dada2_sample_inference:
    input:
    # Dereplicated (aka unique) sequences of the sample
        derep="results/dada2/uniques/{run}/{sample}_{orientation}.RDS",
        err="results/dada2/learn-errors/{run}/model_{orientation}.RDS" # Error model
    output:
        "results/dada2/denoised/{run}/{sample}_{orientation}.RDS" # Inferred sample composition
    log:
        "logs/dada2/sample-inference/{run}/{sample}_{orientation}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/sample-inference"

rule dada2_merge_pairs:
    input:
      dadaF="results/dada2/denoised/{run}/{sample}_R1.RDS",# Inferred composition
      dadaR="results/dada2/denoised/{run}/{sample}_R2.RDS",
      derepF="results/dada2/uniques/{run}/{sample}_R1.RDS",# Dereplicated sequences
      derepR="results/dada2/uniques/{run}/{sample}_R2.RDS"
    output:
        "results/dada2/merged/{run}/{sample}.RDS"
    log:
        "logs/dada2/merge-pairs/{run}/{sample}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/merge-pairs"

rule dada2_make_table_pe:
    input:
        expand("results/dada2/merged/{{run}}/{sample}.RDS", sample = SAMPLES)
    output:
        "results/dada2/seqtab/{run}/seqtab-pe.RDS"
    params:
        names = SAMPLES, # Sample names instead of paths
        orderBy = "nsamples" # Change the ordering of samples
    log:
        "logs/dada2/make-table/{run}/make-table-pe.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/make-table"

rule export_seqtab_to_fasta:
    input:
        "results/dada2/merged/20190508_0074/{sample}.RDS"
    output:
        "results/dada2/merged/20190508_0074/{sample}-seqtab-merged.fa"
    shell:
        "./scripts/export_seqtab_to_fasta.R {input} {output}"

# # rule kraken2:
# #     input:
# #         fasta ="results/dada2/merged/20190508_0074/{sample}-seqtab-merged.fa",
# #         db = "/users/work/cat3/db/kraken2/silva"
# #     output:
# #         report = "results/kraken2/20190508_0074/{sample}-report.txt",
# #         stdout = "results/kraken2/20190508_0074/{sample}-kraken2-stdout.txt",
# #         stderr = "results/kraken2/20190508_0074/{sample}-kraken2-stderr.txt"
# #     shell:
# #         "kraken2 --db {input.db} --threads 4 --report {output.report} {input.fasta} 1> {output.stdout} 2> {output.stderr}"
# #
rule dada2_remove_chimeras:
    input:
        "results/dada2/seqtab/{run}/seqtab-pe.RDS" # Sequence table
    output:
        "results/dada2/seqtab/{run}/seqtab.nochimeras.RDS" # Chimera-free sequence table
    log:
        "logs/dada2/remove-chimeras/{run}/remove-chimeras.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/remove-chimeras"

rule dada2_collapse_nomismatch:
    input:
        "results/dada2/seqtab/{run}/seqtab.nochimeras.RDS" # Chimera-free sequence table
    output:
        "results/dada2/seqtab/{run}/seqtab.collapsed.RDS"
    log:
        "logs/dada2/collapse-nomismatch/{run}/collapse-nomismatch.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/collapse-nomismatch"

rule dada2_assign_taxonomy:
    input:
        seqs="results/dada2/seqtab/{run}/seqtab.collapsed.RDS", # Chimera-free sequence table
        refFasta="/users/work/cat3/db/dada2/silva_nr99_v138_wSpecies_train_set.fa.gz" # Reference FASTA for taxonomy
    output:
        "results/dada2/taxa/{run}/taxa.RDS" # Taxonomic assignments
    log:
        "logs/dada2/assign-taxonomy/{run}/assign-taxonomy.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/dada2/assign-taxonomy"
