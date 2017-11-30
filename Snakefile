#!/usr/bin/env python

configfile: "config.yaml"

CMSAMPLES = {k:v for (k,v) in config["samples"].items() if v["group"]=="CM"}

rule all:
    input:
        expand("reformatted/{sample}-{strand}.{fmt}", sample=CMSAMPLES, strand=["plus","minus"], fmt=["wig", "bw"])

rule wig_split_strands:
    input:
        wig = lambda wildcards: CMSAMPLES[wildcards.sample]["path"]
    output:
        plus = "reformatted/{sample}-plus.wig",
        minus = "reformatted/{sample}-minus.wig"
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} $1~/^track/ {{print $0 > "{output.plus}"; print $0 > "{output.minus}"; next}} $1~/^variableStep/ {{print $0 > "{output.plus}"; print $0 > "{output.minus}"; next}} $2>=0 {{print $0 > "{output.plus}"}} $2<0 {{print $1, -$2 > "{output.minus}"}}' {input}
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
