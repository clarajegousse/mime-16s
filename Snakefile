#configfile: "config.yaml"
SAMPLES = ["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# SAMPLES, = glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")

rule cutadapt:
    input:
        expand(["data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz",
        "data/miseq/20190508_0074/{sample}_L001_R2_001.fastq.gz"], sample = SAMPLES)
    output:
        fastq1 = expand("results/cutadapt/20190508_0074/{sample}_L001_R1_001.fastq.gz", sample = SAMPLES),

        fastq2 = expand("results/cutadapt/20190508_0074/{sample}_L001_R2_001.fastq.gz", sample = SAMPLES),

        qc = expand("results/cutadapt/20190508_0074/{sample}.qc.txt", sample = SAMPLES)

    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 1 -q 20"
    log:
        expand("logs/cutadapt/{sample}.log", sample = SAMPLES)
    threads: 4 # set desired number of threads here
    wrapper:
        "0.70.0/bio/cutadapt/pe"
