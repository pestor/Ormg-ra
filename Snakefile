configfile: "files.yaml"
include: "config.py"

rule all:
    input:
        "report.html",


rule bwa_map:
    input:
        IDX,
	fwd="/mnt/lustre/Active_Projects/MODY_targeted_NGS/Public/140512_M02318_0007_000000000-A7M6L/Data/Intensities/BaseCalls/{sample}_L001_R1_001.fastq.gz",
	rev="/mnt/lustre/Active_Projects/MODY_targeted_NGS/Public/140512_M02318_0007_000000000-A7M6L/Data/Intensities/BaseCalls/{sample}_L001_R2_001.fastq.gz"
    output:
        temp("mapped_reads/{sample}.bam")
    threads: 8
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}"
    log:
        "logs/bwa_map/{sample}.log"
    shell:
        "(BWA_PATH/bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Shb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    threads: 8
    shell:
        "{SAMTOOLS_PATH}/samtools sort -@ {threads} {input} -f {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule generate_dag:
    input:
        "Snakefile"
    output:
        "dag.svg"
    shell:
        "snakemake --dag | dot -Tsvg > dag.svg"

rule haploC_call:
    input:
        fa=IDX,
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai"
    output:
        "calls/{sample}.raw.vcf"
    log:
        "logs/{sample}.gatk"
    shell:	"java -jar -Xmx32g {GATK_PATH}/GenomeAnalysisTK.jar "
    	"-T HaplotypeCaller "
    	"-R {input.fa} "
    	"-I {input.bam} "
    	"--genotyping_mode DISCOVERY "
    	"-stand_emit_conf 10 "
    	"-stand_call_conf 30 "
    	"-o {output} 2> {log}"


def _gatk_multi_arg(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)

rule combineGVCFs:
    input:
        fa=IDX,
        vcfs=expand("calls/{sample}.raw.vcf", sample=config["sample"])
    output:
        "Combined.vcf"
    run:
        vcfs = _gatk_multi_arg("--variant ",input.vcfs)
        shell(
            "java -jar {GATK_PATH}/GenomeAnalysisTK.jar "
            "-T CombineGVCFs "
            "-R {IDX} "
            "{vcfs} "
            "-o {output} ")

rule report:
    input:
        "Combined.vcf",
        "dag.svg"
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        T0

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0], T0=input[1])
