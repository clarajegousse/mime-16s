#configfile: "config.yaml"
# SAMPLES = ["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# SAMPLES, = glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")

rule cutadapt:
    input:
        ["data/miseq/20190508_0074/FX003-016-16S-V4_S58_L001_R1_001.fastq.gz",
        "data/miseq/20190508_0074/FX003-016-16S-V4_S58_L001_R2_001.fastq.gz"])
    output:
        fastq1 ="results/cutadapt/20190508_0074/FX003-016-16S-V4_S58_L001_R1_001.fastq.gz",

        fastq2 = "results/cutadapt/20190508_0074/FX003-016-16S-V4_S58_L001_R2_001.fastq.gz",

        qc = "results/cutadapt/20190508_0074/FX003-016-16S-V4_S58.qc.txt",

    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 1 -q 20"
    log:
        "logs/cutadapt/FX003-016-16S-V4_S58.log"
    threads: 4 # set desired number of threads here
    wrapper:
        "0.70.0/bio/cutadapt/pe"
