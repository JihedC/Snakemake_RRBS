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

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.samples), ["fq1", "fq2"]].dropna()

##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "[A-Za-z0-9]+"


rule all:
	message:
		"All done!"
	input:
			expand(RESULT_DIR + "methyl/{sample}_1_val_1_bismark_bt2_pe.bedGraph.gz", sample=SAMPLES),


rule trim:
	input:
        	get_fastq
	output:
		sample1 = RESULT_DIR + "trimmed/{sample}_1_val_1.fq",
		sample2 = RESULT_DIR + "trimmed/{sample}_2_val_2.fq"
	conda:
		"envs.yaml"
	shell:
        	"trim_galore --rrbs --paired -o trimmed/ {input.read} 2> trim.log"


rule CpGisland_finder:
	input:
        	GENOME_DIR + "mm10.reference.fa"
	output:
        	RESULT_DIR + "CpGislands.gtf"
	shell:
        	"python CpGislandsFind.py {input} > {output}"


rule gtf_to_bed:
	input:
        	RESULT_DIR + "CpGislands.gtf"
	output:
		RESULT_DIR + "CpGislands.bed"
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
		sample1 = RESULT_DIR + "trimmed/{sample}_1_val_1.fq",
		sample2 = RESULT_DIR + "trimmed/{sample}_2_val_2.fq",
		genome = GENOME_DIR
	output:
		RESULT_DIR + DIR + "bismark/{sample}_1_val_1_bismark_bt2_pe.bam",
		RESULT_DIR + DIR + "bismark/{sample}_1_val_1_bismark_bt2_PE_report.txt"
	conda:
		"envs.yaml"
	shell:
		"bismark --genome {input.genome} -1 {input.sample1} -2 {input.sample2} 2> alig.log"


rule methyl_ex:
	input:
        	RESULT_DIR + DIR + "bismark/{sample}_1_val_1_bismark_bt2_pe.bam",
		genome = GENOME_DIR
	output:
        	RESULT_DIR + "methyl/{sample}_1_val_1_bismark_bt2_pe.bedGraph.gz"
	conda:
		"envs.yaml"
	shell:
		"bismark_methylation_extractor --paired-end --zero_based --remove_spaces --ignore_r2 2 --bedGraph --cytosine_report --buffer_size 10G --genome_folder {input.genome} -o methyl/ {input.sample} 2> met_ex.log"
