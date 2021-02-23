
# snakemake --cluster qsub -j 12 --latency-wait 60 --rerun-incomplete --printshellcmds

import pandas as pd

configfile: "config.yaml"

SampleTable = pd.read_csv(config['SAMPLETABLE'], sep = ",", header = 0, index_col = 0)
SAMPLES = list(SampleTable.index)

ORIENTATION = config["ORIENTATION"]
RUN = config["RUN"]

rule all:
    input:
        expand("results/{run}/cutadapt/{sample}-qc-report.txt", run = RUN, sample = SAMPLES),
        expand("logs/dada2/{run}/dada2-filter.log", run = RUN),
        expand("logs/dada2/{run}/dada2-inference.log", run = RUN),
        expand("logs/dada2/{run}/dada2-merge-chimera-taxo.log", run = RUN)

        #expand("results/dada2/quality-profile/{run}/{sample}-quality-profile.png", run = RUN, sample = SAMPLES),
        #expand("results/dada2/filtered_trim_pe/{run}/{sample}.tsv", run = RUN, sample = SAMPLES),
        #expand("reports/dada2/learn-errors/{run}/errors_{orientation}.png", run = RUN, orientation = ORIENTATION),
        #expand("results/dada2/denoised/{run}/{sample}_{orientation}.RDS", run = RUN, sample = SAMPLES, orientation = ORIENTATION),
        #expand("results/dada2/merged/{run}/{sample}.RDS", run = RUN, sample = SAMPLES),
        #expand("results/dada2/taxa/{run}/taxa.RDS", run = RUN),
        #expand("results/dada2/final/{run}-ASVs.fa", run = RUN),
        #expand("results/kraken2/{run}/{sample}-report.txt", run = RUN, sample = SAMPLES)

rule cutadapt:
    input:
        fwd = "data/miseq/{run}/{sample}_R1.fastq.gz",
        rev = "data/miseq/{run}/{sample}_R2.fastq.gz",
    output:
        fwd = "results/{run}/cutadapt/{sample}_R1.fastq.gz",
        rev = "results/{run}/cutadapt/{sample}_R2.fastq.gz",
        report = "results/{run}/cutadapt/{sample}-qc-report.txt"
    params:
        adapter_a = "^ACGGGGYGCAGCAGGCGCGA...AGCGCGGACGACGYGGGGCA",
        adapter_A = "^AGCGCGGACGACGYGGGGCA...ACGGGGYGCAGCAGGCGCGA",
        minimum_length = 90,
        maximum_length = 280
    log:
        "logs/cutadapt/{run}/{sample}.log"
    shell:
        "cutadapt -a {params.adapter_a} -A {params.adapter_A} \
         -m {params.minimum_length} -M {params.maximum_length} \
         -o {output.fwd} -p {output.rev} \
          {input.fwd} {input.rev} \
          2> {output.report}"

          # quick and dirty way to check the average number of reads written
          # cat snakejob.cutadapt.*.sh.o* | grep "Pairs written" | cut -f 3 -d "(" | sed 's/%)//' | awk '{ total += $1 } END { print total/NR }'

rule dada2_filter:
    input:
        path = expand("results/{run}/cutadapt/", run = RUN)
    output:
        path = directory(expand("results/{run}/dada2-filter", run = RUN))
    params:
        trunc_len_fwd = 135,
        trunc_len_rev = 135,
        maxEE = 'Inf',
        truncQ = 11,
        trimLeft = 35
    log:
        directory(expand("logs/dada2/{run}/dada2-filter.log", run = RUN))
    shell:
        #r"""./scripts/dada2-filter.r --input_path {input.path} \
        #--output_path {output.path} \
        #--trunc_len_fwd {params.trunc_len_fwd} \
        #--trunc_len_rev {params.trunc_len_rev} \
        #--maxee {params.maxEE} \
        #--truncq {params.truncQ}"""
        "./scripts/dada2-filter.r \
        --input_path {input.path} \
        --output_path {output.path} \
        --fwd_trunc_len {params.trunc_len_fwd} \
        --rev_trunc_len {params.trunc_len_rev} \
        --maxee {params.maxEE} --truncq {params.truncQ}  --trimleft {params.trimLeft}"

rule dada2_inference:
    input:
        path = expand("results/{run}/dada2-filter/", run = RUN)
    output:
        path = directory(expand("results/{run}/dada2-inference/", run = RUN))
    log:
        directory(expand("logs/dada2/{run}/dada2-inference.log", run = RUN))
    shell:
        "./scripts/dada2-inference.r \
        --input_path {input.path} \
        --output_path {output.path}"

rule dada2_merge_chimera_taxo:
    input:
        path = expand("results/{run}/dada2-inference/", run = RUN)
    output:
        path = directory(expand("results/{run}/dada2-merge-chimera-taxo/", run = RUN))
    log:
        directory(expand("logs/dada2/{run}/dada2-merge-chimera-taxo.log", run = RUN))
    shell:
        "./scripts/dada2-merge-chimera-taxo.r \
        --input_path {input.path} \
        --output_path {output.path}"

# rule dada2_quality_profile_pe:
#     input:
#         expand("results/trimmed/{{run}}/{{sample}}_{orientation}.fastq.gz", orientation = ORIENTATION)
#     output:
#         "results/dada2/quality-profile/{run}/{sample}-quality-profile.png"
#     log:
#         "logs/dada2/quality-profile/{run}/{sample}-quality-profile-pe.log"
#     wrapper:
#         "0.70.0/bio/dada2/quality-profile"
#
# rule dada2_filter_trim_pe:
#     input:
#         # Paired-end files without primer sequences
#         fwd="results/trimmed/{run}/{sample}_R1.fastq.gz",
#         rev="results/trimmed/{run}/{sample}_R2.fastq.gz"
#     output:
#         filt="results/dada2/filtered_trim_pe/{run}/{sample}_R1.fastq.gz",
#         filt_rev="results/dada2/filtered_trim_pe/{run}/{sample}_R2.fastq.gz",
#         stats="results/dada2/filtered_trim_pe/{run}/{sample}.tsv"
#     params:
#         # Set the maximum expected errors tolerated in filtered reads
#         maxEE=2,
#         # Set the number of kept bases in forward and reverse reads
#         truncLen=[240,200]
#     log:
#         "logs/dada2/filter-trim-pe/{run}/{sample}.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/filter-trim"
#
# rule dada2_learn_errors:
#     input:
#     # Quality filtered and trimmed forward FASTQ files (potentially compressed)
#         expand("results/dada2/filtered_trim_pe/{{run}}/{sample}_{{orientation}}.fastq.gz", sample = SAMPLES)
#     output:
#         err="results/dada2/learn-errors/{run}/model_{orientation}.RDS",# save the error model
#         plot="reports/dada2/learn-errors/{run}/errors_{orientation}.png",# plot observed and estimated rates
#     params:
#         randomize=True
#     log:
#         "logs/dada2/learn-errors/{run}/learn-errors_{orientation}.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/learn-errors"
#
# rule dada2_dereplicate_fastq:
#     input:
#     # Quality filtered FASTQ file
#         "results/dada2/filtered_trim_pe/{run}/{fastq}.fastq.gz"
#     output:
#     # Dereplicated sequences stored as `derep-class` object in a RDS file
#         "results/dada2/uniques/{run}/{fastq}.RDS"
#     log:
#         "logs/dada2/dereplicate-fastq/{run}/{fastq}.log"
#     wrapper:
#         "0.70.0/bio/dada2/dereplicate-fastq"
#
# rule dada2_sample_inference:
#     input:
#     # Dereplicated (aka unique) sequences of the sample
#         derep="results/dada2/uniques/{run}/{sample}_{orientation}.RDS",
#         err="results/dada2/learn-errors/{run}/model_{orientation}.RDS" # Error model
#     output:
#         "results/dada2/denoised/{run}/{sample}_{orientation}.RDS" # Inferred sample composition
#     log:
#         "logs/dada2/sample-inference/{run}/{sample}_{orientation}.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/sample-inference"
#
# rule dada2_merge_pairs:
#     input:
#       dadaF="results/dada2/denoised/{run}/{sample}_R1.RDS",# Inferred composition
#       dadaR="results/dada2/denoised/{run}/{sample}_R2.RDS",
#       derepF="results/dada2/uniques/{run}/{sample}_R1.RDS",# Dereplicated sequences
#       derepR="results/dada2/uniques/{run}/{sample}_R2.RDS"
#     output:
#         "results/dada2/merged/{run}/{sample}.RDS"
#     log:
#         "logs/dada2/merge-pairs/{run}/{sample}.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/merge-pairs"
#
# rule dada2_make_table_pe:
#     input:
#         expand("results/dada2/merged/{{run}}/{sample}.RDS", sample = SAMPLES)
#     output:
#         "results/dada2/seqtab/{run}/seqtab-pe.RDS"
#     params:
#         names = SAMPLES, # Sample names instead of paths
#         orderBy = "nsamples" # Change the ordering of samples
#     log:
#         "logs/dada2/make-table/{run}/make-table-pe.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/make-table"
#
# rule export_seqtab_to_fasta:
#     input:
#         "results/dada2/merged/{run}/{sample}.RDS"
#     output:
#         "results/dada2/merged/{run}/{sample}-seqtab-merged.fa"
#     shell:
#         "./scripts/export_seqtab_to_fasta.R {input} {output}"
#
# rule kraken2:
#     input:
#         fasta ="results/dada2/merged/{run}/{sample}-seqtab-merged.fa",
#         db = "/users/work/cat3/db/kraken2/silva"
#     output:
#         report = "results/kraken2/{run}/{sample}-report.txt",
#         stdout = "results/kraken2/{run}/{sample}-kraken2-stdout.txt",
#         stderr = "results/kraken2/{run}/{sample}-kraken2-stderr.txt"
#     shell:
#         "kraken2 --db {input.db} --threads 4 --report {output.report} {input.fasta} 1> {output.stdout} 2> {output.stderr}"
#
# rule dada2_remove_chimeras:
#     input:
#         "results/dada2/seqtab/{run}/seqtab-pe.RDS" # Sequence table
#     output:
#         "results/dada2/seqtab/{run}/seqtab.nochimeras.RDS" # Chimera-free sequence table
#     log:
#         "logs/dada2/remove-chimeras/{run}/remove-chimeras.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/remove-chimeras"
#
# rule dada2_collapse_nomismatch:
#     input:
#         "results/dada2/seqtab/{run}/seqtab.nochimeras.RDS" # Chimera-free sequence table
#     output:
#         "results/dada2/seqtab/{run}/seqtab.collapsed.RDS"
#     log:
#         "logs/dada2/collapse-nomismatch/{run}/collapse-nomismatch.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/collapse-nomismatch"
#
# rule dada2_assign_taxonomy:
#     input:
#         seqs="results/dada2/seqtab/{run}/seqtab.collapsed.RDS", # Chimera-free sequence table
#         refFasta="/users/work/cat3/db/dada2/silva_nr99_v138_wSpecies_train_set.fa.gz" # Reference FASTA for taxonomy
#     output:
#         "results/dada2/taxa/{run}/taxa.RDS" # Taxonomic assignments
#     log:
#         "logs/dada2/assign-taxonomy/{run}/assign-taxonomy.log"
#     threads: 1 # set desired number of threads here
#     wrapper:
#         "0.70.0/bio/dada2/assign-taxonomy"
#
# rule extract_dada2_results:
#     input:
#         seqtab = "results/dada2/seqtab/{run}/seqtab.nochimeras.RDS",
#         #taxo = "results/dada2/taxa/{run}/taxa.RDS"
#     output:
#         asv_seq = "results/dada2/final/{run}-ASVs.fa",
#         asv_counts = "results/dada2/final/{run}-ASVs_counts.tsv"
#     shell:
#         "./scripts/extract_dada2_results.R {input.seqtab} {output.asv_seq} {output.asv_counts}"
