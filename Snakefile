#!/usr/bin/env python

configfile: "config.yaml"

CMSAMPLES = {k:v for (k,v) in config["samples"].items() if v["group"]=="CM"}
SDSAMPLES = {k:v for (k,v) in config["samples"].items() if v["group"]=="SD"}

rule all:
    input:
        expand("reformatted/{sample}-{strand}.{fmt}", sample=CMSAMPLES, strand=["plus","minus"], fmt=["wig", "bw", "bedgraph"]),
        expand("reformatted/{sample}-SENSE.bedgraph", sample=CMSAMPLES),
        "union-bedgraph-allsamples.tsv.gz",
        "doris-v-malabat-tss-correlations.svg"

rule wig_split_strands:
    input:
        wig = lambda wildcards: CMSAMPLES[wildcards.sample]["path"]
    output:
        plus = "reformatted/{sample}-plus.wig",
        minus = "reformatted/{sample}-minus.wig"
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} $1=="variableStep chrom=pombe_AB325691 span=1" {{exit}} $1~/^track/ {{print $0 > "{output.plus}"; print $0 > "{output.minus}"; next}} $1~/^variableStep/ {{print $0 > "{output.plus}"; print $0 > "{output.minus}"; next}} $2>=0 {{print $0 > "{output.plus}"}} $2<0 {{print $1, -$2 > "{output.minus}"}}' {input}
        """

rule wig_to_bigwig:
    input:
        wig = "reformatted/{sample}-{strand}.wig",
        chrsizes = config["chrsizes"]
    output:
        bw = "reformatted/{sample}-{strand}.bw"
    shell: """
        wigToBigWig <(sed 's/chrmt/chrM/g' {input.wig}) {input.chrsizes} {output.bw}
        """

rule bigwig_to_bedgraph:
    input:
        bw = "reformatted/{sample}-{strand}.bw"
    output:
        bg = "reformatted/{sample}-{strand}.bedgraph"
    shell: """
        bigWigToBedGraph {input.bw} {output.bg}
        """

rule make_stranded_bedgraph:
    input:
        plus = "reformatted/{sample}-plus.bedgraph",
        minus = "reformatted/{sample}-minus.bedgraph",
    output:
        sense = "reformatted/{sample}-SENSE.bedgraph"
    shell: """
        bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}
        """

rule union_bedgraph:
    input:
        sd = [v["path"] for k,v in SDSAMPLES.items()],
        cm = expand("reformatted/{sample}-SENSE.bedgraph", sample=CMSAMPLES)
    params:
        names = list(SDSAMPLES.keys()) + list(CMSAMPLES.keys())
    output:
        "union-bedgraph-allsamples.tsv.gz"
    shell: """
        bedtools unionbedg -i {input.sd} {input.cm} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output}
        """

rule plotcorrelations:
    input:
        "union-bedgraph-allsamples.tsv.gz"
    output:
        "doris-v-malabat-tss-correlations.svg"
    params:
        pcount = 0.1,
        samplelist = list(SDSAMPLES.keys()) + list(CMSAMPLES.keys())
    script:
        "scripts/plotcorr.R"

rule deeptools_matrix:
    input:
        annotation = "../genome/annotations/stranded/{annotation}-STRANDED.bed",
        bw = "coverage/{norm}/bw/{sample}-tss-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-{norm}-{strand}.log"
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
