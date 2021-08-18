rule gatk_calling:
    input:
        # single or list of bam files
        bam="results/recal/{sample}.{aligner}.sorted.bam",
        ref="resources/genome.fasta",
#        known="resources/variation.noiupac.vcf.gz"
    output:
        vcf=temp("results/candidate_calls/haplotypecaller/{sample}.{aligner}.vcf"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{aligner}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx4G -XX:ParallelGCThreads=10", # optional
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/gatk_calling-{sample}-{aligner}.benchmark"
    shell:
        "gatk --java-options '{params.java_opts}' HaplotypeCaller {params.extra} "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.vcf} > {log} 2>&1"


rule gatk_vcf_to_bcf:
    input:
        "results/candidate_calls/haplotypecaller/{sample}.{aligner}.vcf"
    output:
        temp("results/candidate_calls/haplotypecaller/{sample}.{aligner}.bcf")
    conda: "../envs/bcftools.yaml"
    shell:
        "bcftools view {input} -Ob > {output}"


rule gatk_bcftools_reheader:
    input:
        vcf="results/candidate_calls/haplotypecaller/{sample}.{aligner}.bcf",
        ## new header, can be omitted if "samples" is set
        #header="header.txt",
        ## file containing new sample names, can be omitted if "header" is set
        samples="samples.tsv"
    output:
        temp("results/candidate_calls/haplotypecaller/{sample}.{aligner}.reheader.bcf")
    params:
        extra="",  # optional parameters for bcftools reheader
        view_extra="-O b"  # add output format for internal bcftools view call
    wrapper:
        "0.74.0/bio/bcftools/reheader"



rule index_gatk_haplotypecaller:
    input:
        "results/candidate_calls/haplotypecaller/{sample}.{aligner}.reheader.bcf"
    output:
        temp("results/candidate_calls/haplotypecaller/{sample}.{aligner}.reheader.bcf.csi")
    params:
        extra=""  # optional parameters for bcftools index
    wrapper:
        "0.74.0/bio/bcftools/index"


rule deep_variant_calling:
    input:
        bam="results/recal/{sample}.{aligner}.sorted.bam",
        ref="resources/genome.fasta"
    output:
        vcf=temp("results/candidate_calls/deep_variant/{sample}.{aligner}.vcf.gz")
    params:
        model="wes",   # {wgs, wes, pacbio, hybrid}
        #sample_name=lambda w: w.sample, # optional
        extra=""
    threads: 24
    log:
        "logs/deepvariant/{sample}.{aligner}/stdout.log"
    benchmark:
        "benchmarks/deep_variant_calling-{sample}-{aligner}.benchmark"
    wrapper:
        "0.74.0/bio/deepvariant"


rule deepvariant_bcftools_reheader:
    input:
        vcf="results/candidate_calls/deep_variant/{sample}.{aligner}.vcf.gz",
        ## new header, can be omitted if "samples" is set
        #header="header.txt",
        ## file containing new sample names, can be omitted if "header" is set
        samples="samples.tsv"
    output:
        temp("results/candidate_calls/deep_variant/{sample}.{aligner}.reheader.bcf")
    params:
        extra="",  # optional parameters for bcftools reheader
        view_extra="-O b"  # add output format for internal bcftools view call
    wrapper:
        "0.74.0/bio/bcftools/reheader"


rule index_deep_variant:
    input:
        "results/candidate_calls/deep_variant/{sample}.{aligner}.reheader.bcf"
    output:
        temp("results/candidate_calls/deep_variant/{sample}.{aligner}.reheader.bcf.csi")
    params:
        extra=""  # optional parameters for bcftools index
    wrapper:
        "0.74.0/bio/bcftools/index"


rule freebayes:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        #regions=get_regions(),
        # you can have a list of samples here
        samples="results/recal/{sample}.{aligner}.sorted.bam",
        index="results/recal/{sample}.{aligner}.sorted.bai",
    output:
        temp("results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.bcf"),
    log:
        "logs/freebayes/{sample}.{aligner}.log",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra="--pooled-continuous --min-alternate-count 1 --min-alternate-fraction {}".format(
            config["params"]["freebayes"].get("min_alternate_fraction", "0.05")
        ),
    threads: 24
    benchmark:
        "benchmarks/freebayes-{sample}-{aligner}.benchmark"
    wrapper:
        "0.68.0/bio/freebayes"


rule freebayes_bcftools_reheader:
    input:
        vcf="results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.bcf",
        ## new header, can be omitted if "samples" is set
        #header="header.txt",
        ## file containing new sample names, can be omitted if "header" is set
        samples="samples.tsv"
    output:
        temp("results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.reheader.bcf")
    params:
        extra="",  # optional parameters for bcftools reheader
        view_extra="-O b"  # add output format for internal bcftools view call
    wrapper:
        "0.74.0/bio/bcftools/reheader"


rule index_freebayes:
    input:
        "results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.reheader.bcf"
    output:
        temp("results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.reheader.bcf.csi")
    params:
        extra=""  # optional parameters for bcftools index
    wrapper:
        "0.74.0/bio/bcftools/index"


rule delly:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        samples="results/recal/{sample}.{aligner}.sorted.bam",
        index="results/recal/{sample}.{aligner}.sorted.bai",
        #exclude=get_excluded_regions(),
    output:
        temp("results/candidate_calls/delly/{sample}.{aligner}.delly.bcf"),
    log:
        "logs/delly/{sample}.{aligner}.log",
    params:
        extra=config["params"].get("delly", ""),
    threads: 24
    benchmark:
        "benchmarks/delly-{sample}-{aligner}.benchmark"
    wrapper:
        "0.68.0/bio/delly"


rule delly_bcftools_reheader:
    input:
        vcf="results/candidate_calls/delly/{sample}.{aligner}.delly.bcf",
        ## new header, can be omitted if "samples" is set
        #header="header.txt",
        ## file containing new sample names, can be omitted if "header" is set
        samples="samples.tsv"
    output:
        temp("results/candidate_calls/delly/{sample}.{aligner}.delly.reheader.bcf")
    params:
        extra="",  # optional parameters for bcftools reheader
        view_extra="-O b"  # add output format for internal bcftools view call
    wrapper:
        "0.74.0/bio/bcftools/reheader"


rule index_delly:
    input:
        "results/candidate_calls/delly/{sample}.{aligner}.delly.reheader.bcf"
    output:
        temp("results/candidate_calls/delly/{sample}.{aligner}.delly.reheader.bcf.csi")
    params:
        extra=""  # optional parameters for bcftools index
    wrapper:
        "0.74.0/bio/bcftools/index"
