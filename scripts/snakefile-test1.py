#configfile: "config.yaml"

# SAMPLES = ["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# RUNS, SAMPLES, = glob_wildcards("data/miseq/{run}/{sample}_L001_R1_001.fastq.gz")
SAMPLES, = glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")
RUNS = ["20190508_0074"]

rule all:
    input: expand("results/trimmomatic/{run}/{sample}_unpaired_L001_R2_001.fastq.gz", sample = SAMPLES, run = RUNS)

rule cutadapt:
    input:
        fwd = "data/miseq/{run}/{sample}_L001_R1_001.fastq.gz",
        rev = "data/miseq/{run}/{sample}_L001_R2_001.fastq.gz",
    output:
        fwd = "results/cutadapt/{run}/{sample}_L001_R1_001.fastq.gz",
        rev = "results/cutadapt/{run}/{sample}_L001_R2_001.fastq.gz",
        report = "results/cutadapt/{run}/{sample}-qc-report.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapter_a = "AGAGCACACGTCTGAACTCCAGTCAC",
        adapter_g = "AGATCGGAAGAGCACACGT",
        adapter_A = "AGAGCACACGTCTGAACTCCAGTCAC",
        adapter_G = "AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        minimum_length = 100
        #quality_cutoff = 20
    log:
        "logs/cutadapt/{run}/{sample}.log"
    shell:
        "cutadapt -a {params.adapter_a} -A {params.adapter_A} \
         -m {params.minimum_length} \
         -o {output.fwd} -p {output.rev} \
          {input.fwd} {input.rev} \
          2> {output.report}"

rule trimmomatic:
    input:
        fwd = "results/cutadapt/{run}/{sample}_L001_R1_001.fastq.gz",
        rev = "results/cutadapt/{run}/{sample}_L001_R2_001.fastq.gz"
    output:
        fwd_paired = "results/trimmomatic/{run}/{sample}_paired_L001_R1_001.fastq.gz",
        rev_paired = "results/trimmomatic/{run}/{sample}_paired_L001_R2_001.fastq.gz",
        fwd_unpaired = "results/trimmomatic/{run}/{sample}_unpaired_L001_R1_001.fastq.gz",
        rev_unpaired = "results/trimmomatic/{run}/{sample}_unpaired_L001_R2_001.fastq.gz"
    params:
        leading = 5,
        trailing = 5,
        slidingwindow = "4:15",
        minlen = 100
        #illuminaclip = "TruSeq3-PE.fa:2:30:10"
    log:
        "logs/trimmomatic/{run}/{sample}.log"
    shell:
        "trimmomatic PE \
        -phred33 {input.fwd} {input.rev} \
        {output.fwd_paired} {output.fwd_unpaired} \
        {output.rev_paired} {output.rev_unpaired} \
        LEADING:{params.leading} TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen}"
