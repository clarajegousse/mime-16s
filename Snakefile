#configfile: "config.yaml"
# SAMPLES = ["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# SAMPLES, = glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")

rule cutadapt:
    input:
        input_r1 = "data/miseq/20190508_0074/FX003-016-16S-V4_S58_L001_R1_001.fastq.gz",
        input_r2 = "data/miseq/20190508_0074/FX003-016-16S-V4_S58_L001_R2_001.fastq.gz"
    output:
        output_r1 = "results/cutadapt/20190508_0074/FX003-016-16S-V4_S58_L001_R1_001.fastq.gz",
        output_r2 = "results/cutadapt/20190508_0074/FX003-016-16S-V4_S58_L001_R2_001.fastq.gz",
        report = "results/cutadapt/20190508_0074/FX003-016-16S-V4-qc-report.txt"
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
        "logs/cutadapt/FX003-016-16S-V4_S58.log"
    shell:
        "cutadapt -a {adapter_a} -A {adapter_A} -o {output_r1} -p {output_r2} {input_r1} {input_r2} 2> {report}.txt"
