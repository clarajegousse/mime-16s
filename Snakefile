#configfile: "config.yaml"
# SAMPLES=["FX003-016-16S-V4_S58", "FX003-017-16S-V4_S59"]
# SAMPLES,=glob_wildcards("data/miseq/20190508_0074/{sample}_L001_R1_001.fastq.gz")

rule cutadapt:
    input:
        input_r1="data/miseq/20190508_0074/FX003-016-16S-V4_S58_L001_R1_001.fastq.gz",
        input_r2="data/miseq/20190508_0074/FX003-016-16S-V4_S58_L001_R2_001.fastq.gz"
    output:
        output_r1="results/cutadapt/20190508_0074/FX003-016-16S-V4_S58_L001_R1_001.fastq.gz",
        output_r2="results/cutadapt/20190508_0074/FX003-016-16S-V4_S58_L001_R2_001.fastq.gz",
        report="results/cutadapt/20190508_0074/FX003-016-16S-V4-qc-report.txt"
    params:
        adaptera="AGAGCACACGTCTGAACTCCAGTCAC",
        adapterg="AGATCGGAAGAGCACACGT",
        adapterA="AGAGCACACGTCTGAACTCCAGTCAC",
        adapterG="AGATCGGAAGAGCACACGT"
    log:
        "logs/cutadapt/FX003-016-16S-V4_S58.log"
    shell:
        "cutadapt -a {adaptera} -A {adapterA} -o {output_r1} -p {output_r2} {input_r1} {input_r2} 2> {report}.txt"
