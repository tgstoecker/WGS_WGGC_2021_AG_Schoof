rule samtools_merge_final_bams:
    input:
        expand(["results/recal/{sample}.{aligner}.sorted.bam"], sample=SAMPLES, aligner=ALIGNERS)
    output:
        "results/recal/merged/merged_recal.bam"
    params:
        "-c -p" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        24     # This value - 1 will be sent to -@
    wrapper:
        "0.74.0/bio/samtools/merge"


rule samtools_index:
    input:
        "results/recal/merged/merged_recal.bam"
    output:
        "results/recal/merged/merged_recal.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.74.0/bio/samtools/index"


scattergather:
    calling=128,

rule scatter_merged_candidates:
    input:
        "results/candidates_norm/all_norm.bcf",
    output:
        scatter.calling(
            "results/candidates_norm/{scatteritem}.bcf"
        ),
    log:
        "logs/scatter-candidates/scatter_merged_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


def get_candidate_calls(ext="bcf"):
    def inner(wildcards):
        return expand(
            "results/candidates_norm/{{scatteritem}}.{ext}",
            ext=ext,
        )

    return inner


rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        candidates=get_candidate_calls(),
        bam="results/recal/merged/merged_recal.bam",
        bai="results/recal/merged/merged_recal.bam.bai",
    output:
        "results/observations/observations.{scatteritem}.bcf",
    log:
        "logs/varlociraptor/preprocess/preprocess.{scatteritem}.log",
    benchmark:
        "benchmarks/varlociraptor_preprocess/varlociraptor_preprocess.{scatteritem}.benchmark"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --omit-insert-size --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"


rule sort_scattered_preprocessed:
    input:
        "results/observations/observations.{scatteritem}.bcf"
    output:
        "results/observations/sorted_observations.{scatteritem}.bcf"
    log:
        "logs/observations/sort/sorted_observations.{scatteritem}.log"
    params:
        tmp_dir = "`mktemp -d`",
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf = False,
        # Extra arguments
        extras = ""
    resources:
        mem_mb = 8000
    wrapper:
        "0.77.0/bio/bcftools/sort"


def get_gather_calls_input(ext="bcf"):
    pattern = "results/observations/sorted_observations.{{scatteritem}}.bcf"
    return gather.calling(pattern.format(ext=ext))


rule merge_scattered_preprocessed:
    input:
        calls=get_gather_calls_input(),
    output:
        "results/observations/merged_sorted_observations.bcf",
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    threads: 4
    resources:
        mem_mb=10,
    wrapper:
        "0.77.0/bio/bcftools/concat"


rule varlociraptor_call:
    input:
        obs="results/observations/merged_sorted_observations.bcf",
        scenario="config/scenario.yaml",
    output:
        "results/varlociraptor_calls/calls.bcf"
    log:
        "logs/varlociraptor/call/calls.log",
    params:
        extra=config["params"]["varlociraptor_call"],
    conda:
        "../envs/varlociraptor.yaml"
    benchmark:
        "benchmarks/varlociraptor_call.benchmark"
    shell:
        "varlociraptor "
        "call variants generic --obs NA12878={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

