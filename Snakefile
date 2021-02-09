#configfile: "config.yaml"
RUN = ["20190508_0074"]
SAMPLES = ["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# SAMPLES, = glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")

rule all:
    input: expand("results/cutadapt/{run}/{sample}_L001_R1_001.fastq.gz", sample = SAMPLES, run = RUN)

rule cutadapt:
    input:
        input_r1 = "data/miseq/{run}/{sample}_L001_R1_001.fastq.gz",
        input_r2 = "data/miseq/{run}/{sample}_L001_R2_001.fastq.gz",
    output:
        r1 = "results/cutadapt/{run}/{sample}_L001_R1_001.fastq.gz",
        r2 = "results/cutadapt/{run}/{sample}_L001_R2_001.fastq.gz",
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
        "cutadapt -a {params.adapter_a} -A {params.adapter_A} -m {params.minimum_length} -o {output.r1} -p {output.r2} {input.input_r1} {input.input_r2} 2> {output.report}"
