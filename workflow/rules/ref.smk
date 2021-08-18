rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.59.2/bio/samtools/faidx"

rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "

### Novoalign index
rule register_novoalign_license:
    input:
        "config/licenses/novoalign.lic"
    output:
        touch("config/licenses/novoalign.lic.check")
    conda:
        "../envs/novoalign.yaml"
    shell:
        "novoalign-license-register {input}"

rule novoalign_index:
    input:
        license="config/licenses/novoalign.lic.check",
        genome="resources/genome.fasta",
    output:
        "resources/genome.novoalign.idx"
    conda:
        "../envs/novoalign.yaml"
    shell:
        "novoindex {output} {input.genome}"

### Download pre-process known variants

rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        type="all",
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-variation"

rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    benchmark:
        "benchmarks/remove_iupac_codes.benchmark"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz",
    output:
        "resources/{prefix}.{format}.gz.tbi",
    log:
        "logs/tabix/{prefix}.{format}.log",
#    params:
#        get_tabix_params,
    cache: True
    wrapper:
        "0.59.2/bio/tabix"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.59.2/bio/bwa/index"
