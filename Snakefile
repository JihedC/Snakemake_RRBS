################## Import libraries ##################


import pandas as pd
import os
import sys
from subprocess import call
import itertools
from snakemake.utils import R


################## Configuration file and PATHS ##################

configfile: "config.yaml"

WORKING_DIR             = config["working_dir"]
SAMPLE_DIR              = config["sample_dir"]
GENOME_DIR              = config["genome_dir"]
RESULT_DIR              = config["result_dir"]
DATA_DIR                = config["data_dir"]

################## Samples ##################

SAMPLES, = glob_wildcards(SAMPLE_DIR + "{sample}_1.fastq")

rule all:
	message:
		"All done!"
	input:
			expand(DIR + "methyl/{sample}_1_val_1_bismark_bt2_pe.bedGraph.gz", sample=SAMPLES)
			c
rule trim:
	input:
        DIR + "sample_bs/{sample}_1.fastq",
		DIR + "sample_bs/{sample}_2.fastq"
	output:
        sample1 = "trimmed/{sample}_1_val_1.fq",
		sample2 = "trimmed/{sample}_2_val_2.fq"
    conda:
		"envs.yaml"
	shell:
        "trim_galore --rrbs --paired -o trimmed/ {input} 2> trim.log"


rule CpGisland_finder:
	input:
        GENOME_DIR + "mm10.reference.fa"
	output:
        DIR + "CpGislands.gtf"
	shell:
        "python CpGislandsFind.py {input} > {output}"


rule gtf_to_bed:
	input:
        DIR + "CpGislands.gtf"
	output:
        DIR + "CpGislands.bed"
	shell:
        """
		 awk -F "\\t" '{{OFS="\\t"; print $2, $4, $5, "CpGisland", $6, $7}}' {input} > {output}
	    """

rule genome_prep:
	input:
        GENOME_DIR
	output:
        GENOME_DIR + "Bisulfite_Genome"
    conda:
		"envs.yaml"
	shell:
        "bismark_genome_preparation --verbose {input}"

rule align:
	input:
        sample1 = "trimmed/{sample}_1_val_1.fq",
		sample2 = "trimmed/{sample}_2_val_2.fq",
		genome = GENOME_DIR
	output:
        DIR + "{sample}_1_val_1_bismark_bt2_pe.bam",
		DIR + "{sample}_1_val_1_bismark_bt2_PE_report.txt"
    conda:
		"envs.yaml"
	shell:
        "bismark --genome {input.genome} -1 {input.sample1} -2 {input.sample2} 2> alig.log"

rule clean:
	input:
        DIR + "{sample}_1_val_1_bismark_bt2_PE_report.txt", DIR + "{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        DIR + "BAMs/{sample}_1_val_1_bismark_bt2_SE_report.txt", DIR + "BAMs/{sample}_1_val_1_bismark_bt2.bam"
	shell:
        "mv {input} BAMs/"

rule methyl_ex:
	input:
        sample = DIR + "{sample}_1_val_1_bismark_bt2_pe.bam", genome = GENOME_DIR
	output:
        DIR + "methyl/{sample}_1_val_1_bismark_bt2_pe.bedGraph.gz"
    conda:
		"envs.yaml"
	shell:
        "bismark_methylation_extractor --paired-end --zero_based --remove_spaces --ignore_r2 2 --bedGraph --cytosine_report --buffer_size 10G --genome_folder {input.genome} -o methyl/ {input.sample} 2> met_ex.log"
