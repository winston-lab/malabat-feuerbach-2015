#!/usr/bin/env python

configfile: "config.yaml"

CMSAMPLES = {k:v for (k,v) in config["samples"].items() if v["group"]=="malabat"}
SDSAMPLES = {k:v for (k,v) in config["samples"].items() if v["group"] in ["doris","viktorovskaya"]}

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")


localrules: wig_split_strands, wig_to_bigwig, bigwig_to_bedgraph,
    make_stranded_bedgraph, stranded_bedgraph_to_bigwig, make_stranded_annotations, cat_matrices

rule all:
    input:
        # expand("reformatted/{sample}-{strand}.{fmt}", sample=CMSAMPLES, strand=["plus","minus"], fmt=["wig", "bw", "bedgraph"]),
        # expand("reformatted/{sample}-SENSE.{fmt}", sample=CMSAMPLES, fmt=["bedgraph", "bw"]),
        expand("figures/correlation/tss-seq-v-malabat-window-{windowsize}-tss-correlations.svg", windowsize=config["corr-windowsizes"]),
        # expand("figures/{annotation}/{annotation}-{sample}-SENSE-melted.tsv.gz", annotation=config["annotations"], sample=config["samples"]),
        # expand("figures/{annotation}/allsamples-{annotation}-SENSE.tsv.gz", annotation=config["annotations"]),
        expand("figures/{annotation}/{annotation}-tss-seq-v-malabat.svg", annotation=config["annotations"])

rule wig_split_strands:
    input:
        wig = lambda wildcards: CMSAMPLES[wildcards.sample]["path"]
    output:
        plus = temp("reformatted/{sample}-plus.wig"),
        minus = temp("reformatted/{sample}-minus.wig")
    conda: "envs/malabat15_default.yaml"
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} $1=="variableStep chrom=pombe_AB325691 span=1" {{exit}} $1~/^track/ {{print $0 > "{output.plus}"; print $0 > "{output.minus}"; next}} $1~/^variableStep/ {{print $0 > "{output.plus}"; print $0 > "{output.minus}"; next}} $2>=0 {{print $0 > "{output.plus}"}} $2<0 {{print $1, -$2 > "{output.minus}"}}' {input}
        """

rule wig_to_bigwig:
    input:
        wig = "reformatted/{sample}-{strand}.wig",
        fasta = config["fasta"]
    output:
        bw = "reformatted/{sample}-{strand,plus|minus}.bw"
    conda: "envs/malabat15_default.yaml"
    shell: """
        wigToBigWig <(sed 's/chrmt/chrM/g' {input.wig}) <(faidx {input.fasta} -i chromsizes) {output.bw}
        """

rule bigwig_to_bedgraph:
    input:
        bw = "reformatted/{sample}-{strand}.bw"
    output:
        bg = temp("reformatted/{sample}-{strand,plus|minus}.bedgraph")
    conda: "envs/malabat15_default.yaml"
    shell: """
        bigWigToBedGraph {input.bw} {output.bg}
        """

rule make_stranded_bedgraph:
    input:
        plus = "reformatted/{sample}-plus.bedgraph",
        minus = "reformatted/{sample}-minus.bedgraph",
    output:
        sense = "reformatted/{sample}-SENSE.bedgraph"
    conda: "envs/malabat15_default.yaml"
    shell: """
        bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}
        """

rule map_to_windows:
    input:
        bg = lambda wildcards: SDSAMPLES[wildcards.sample]["path"] if wildcards.sample in SDSAMPLES else "reformatted/" + wildcards.sample + "-SENSE.bedgraph",
        fasta = config["fasta"]
    output:
        exp = temp("reformatted/{sample}-window-{windowsize}-coverage.bedgraph"),
    conda: "envs/malabat15_default.yaml"
    shell: """
        bedtools makewindows -g <(awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' <(faidx {input.fasta} -i chromsizes)) -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}$4=="."{{$4=0}}{{print $0}}' > {output.exp}
        """

rule join_window_counts:
    input:
        coverage = expand("reformatted/{sample}-window-{{windowsize}}-coverage.bedgraph", sample=list(SDSAMPLES.keys())+list(CMSAMPLES.keys())),
    params:
        names = list(SDSAMPLES.keys()) + list(CMSAMPLES.keys())
    output:
        "figures/correlation/union-bedgraph-window-{windowsize}-allsamples.tsv.gz"
    conda: "envs/malabat15_default.yaml"
    shell: """
        bedtools unionbedg -i {input.coverage} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output}
        """

rule plotcorrelations:
    input:
        "figures/correlation/union-bedgraph-window-{windowsize}-allsamples.tsv.gz"
    output:
        "figures/correlation/tss-seq-v-malabat-window-{windowsize}-tss-correlations.svg"
    params:
        pcount = 0.1,
        samplelist = list(SDSAMPLES.keys()) + list(CMSAMPLES.keys())
    conda: "envs/malabat15_tidyverse.yaml"
    script:
        "scripts/plotcorr.R"

rule stranded_bedgraph_to_bigwig:
    input:
        bg = "reformatted/{sample}-SENSE.bedgraph",
        fasta = config["fasta"]
    output:
        "reformatted/{sample}-SENSE.bw"
    conda: "envs/malabat15_default.yaml"
    shell: """
        bedGraphToBigWig {input.bg} <(awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' <(faidx {input.fasta} -i chromsizes)) {output}
        """

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "annotations/{annotation}.bed"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    conda: "envs/malabat15_default.yaml"
    shell: """
        (cut -f1-6 | bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule deeptools_matrix:
    input:
        annotation = "annotations/{annotation}.bed",
        bw = lambda wildcards: "reformatted/" + wildcards.sample + "-SENSE.bw" if config["samples"][wildcards.sample]["group"]=="malabat" else config["samples"][wildcards.sample]["bw"]
    output:
        dtfile = temp("figures/{annotation}/{annotation}_{sample}-SENSE.mat"),
        matrix = temp("figures/{annotation}/{annotation}_{sample}-SENSE.tsv"),
        matrix_gz = "figures/{annotation}/{annotation}_{sample}-SENSE.tsv.gz",
    params:
        scaled_length = lambda wildcards: config["annotations"][wildcards.annotation]["scaled-length"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"] + config["annotations"][wildcards.annotation]["binsize"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"] + config["annotations"][wildcards.annotation]["binsize"],
        # unscaled_5p = lambda wildcards: config["annotations"][wildcards.annotation]["unscaled-5p"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}_{sample}.log"
    conda: "envs/malabat15_default.yaml"
    shell: """
        (computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}
        """

rule melt_matrix:
    input:
        matrix = "figures/{annotation}/{annotation}_{sample}-SENSE.tsv.gz",
    output:
        temp("figures/{annotation}/{annotation}_{sample}-SENSE-melted.tsv.gz")
    params:
        group = lambda wildcards : config["samples"][wildcards.sample]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
    conda: "envs/malabat15_tidyverse.yaml"
    script:
        "scripts/melt_matrix.R"

rule cat_matrices:
    input:
        expand("figures/{{annotation}}/{{annotation}}_{sample}-SENSE-melted.tsv.gz", sample=config["samples"])
    output:
        "figures/{annotation}/allsamples-{annotation}-SENSE.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}.log"
    conda: "envs/malabat15_default.yaml"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_meta:
    input:
        "figures/{annotation}/allsamples-{annotation}-SENSE.tsv.gz"
    params:
        trim_pct = 0.01,
        scaled_length = lambda wildcards: config["annotations"][wildcards.annotation]["scaled-length"],
        ylabel = lambda wildcards: config["annotations"][wildcards.annotation]["ylabel"],
    output:
        "figures/{annotation}/{annotation}-tss-seq-v-malabat.svg"
    conda: "envs/malabat15_tidyverse.yaml"
    script: "scripts/tss_metagene.R"

