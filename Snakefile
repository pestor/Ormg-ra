configfile: "config.json"

rule all:
    input:
        "report.html",


def _get_ref(wildcards):
    return config["references"]["genome"]

def _get_units(pattern):
    def apply(wildcards):
        return expand(
            pattern, reference=wildcards.reference,
            unit=[
                unit for sample in config["samples"].values()
                for unit in sample
            ])
    return apply

rule bwa_map:
    input:
        ref=_get_ref,
	    fwd=_get_units("mapping/units/{unit}.sorted.bam")
    output:
        temp("mapped_reads/{sample}.bam")
    threads: 8
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}"
    log:
        "logs/bwa_map/{sample}.log"
    shell:
        "({config[Software][BWA_PATH]}/bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Shb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    threads: 8
    shell:
        "{config[Software][SAMTOOLS_PATH]}/samtools sort -@ {threads} {input} -f {output}"

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
        ref=_get_ref,
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai"
    output:
        "calls/{sample}.raw.vcf"
    log:
        "logs/{sample}.gatk"
    shell:	"java -jar -Xmx32g {config[Software][GATK_PATH]}/GenomeAnalysisTK.jar "
    	"-T HaplotypeCaller "
    	"-R {input.ref} "
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
        ref=_get_ref,
        vcfs=expand("calls/{sample}.raw.vcf", sample=config["samples"])
    output:
        "Combined.vcf"
    run:
        vcfs = _gatk_multi_arg("--variant ",input.vcfs)
        shell(
            "java -jar {config[Software][GATK_PATH]}/GenomeAnalysisTK.jar "
            "-T CombineGVCFs "
            "-R {input.ref} "
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
