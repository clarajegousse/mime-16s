#configfile: "config.yaml"
SAMPLES = ["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# SAMPLES, = glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")

rule cutadapt:
    input:
        input_r1 = expand("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz", sample = SAMPLES),
        input_r2 = expand("data/miseq/20190508_0074/{sample}_L001_R2_001.fastq.gz", sample = SAMPLES)
    output:
        r1 = expand("results/cutadapt/20190508_0074/{sample}_L001_R1_001.fastq.gz", sample = SAMPLES),
        r2 = expand("results/cutadapt/20190508_0074/{sample}_L001_R2_001.fastq.gz", sample = SAMPLES),
        report = expand("results/cutadapt/20190508_0074/{sample}-qc-report.txt", sample = SAMPLES)
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapter_a = "AGAGCACACGTCTGAACTCCAGTCAC",
        adapter_g = "AGATCGGAAGAGCACACGT",
        adapter_A = "AGAGCACACGTCTGAACTCCAGTCAC",
        adapter_G = "AGATCGGAAGAGCACACGT"
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        #minimum_length = 1,
        #quality-cutoff = 20
    log:
        expand("logs/cutadapt/{sample}.log", sample = SAMPLES)
    shell:
        "cutadapt -a {params.adapter_a} -A {params.adapter_A} -o {output.r1} -p {output.r2} {input.input_r1} {input.input_r2} 2> {output.report}"
