import pandas as pd
import os

# dataset and samples
metadata = pd.read_table(config["metadata"], dtype={"sample": str}, sep = " ").set_index(["sample"], drop=False)
sample_set = metadata["sample"].tolist()

# target rule
rule all:
	input:
		dedup_summary = config['result_dir'] + "all_dedup_summary.out" if config["remove_duplicates"] else [],
		summary = config['result_dir'] + "kaiju_summary.out"
	params: 
		qsub=config["qsub_default"]
	
### 
# Duplicate Removal with cd-hit-dup (optional) ####
###

rule cd_hit_dup:
	input:
               	read1 = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read1"],
               	read2 = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read2"]
	params:
		qsub=config["qsub_default"],
		result_dir = config['result_dir'] 
	output:
		read1_nodup="fastq_files/{sample}_1_nodup.fastq",
		read2_nodup="fastq_files/{sample}_2_nodup.fastq",
		dedup_summary = config['result_dir'] + "{sample}_dedup_summary.out"
	shell:
		'''
		#!/bin/bash
	
		# Dedup search
		mkdir -p fastq_files
		zcat {input.read1} > fastq_files/{wildcards.sample}_1.fastq
		zcat {input.read2} > fastq_files/{wildcards.sample}_2.fastq
		cd-hit-dup -i fastq_files/{wildcards.sample}_1.fastq -i2 fastq_files/{wildcards.sample}_2.fastq -o {output.read1_nodup} -o2 {output.read2_nodup} 
	
		# Dedup summary
		mkdir -p {params.result_dir}
		total_reads=$(wc -l < fastq_files/{wildcards.sample}_1.fastq)
		dedup_reads=$(wc -l < {output.read1_nodup})
		echo $((total_reads/4)) $((dedup_reads/4)) >> {output.dedup_summary}

		# remove temporary fastq files
		rm fastq_files/{wildcards.sample}_1.fastq fastq_files/{wildcards.sample}_2.fastq 
		rm fastq_files/{wildcards.sample}_*.clstr
		'''
  
rule cd_hit_dup_summary:
	input:
		dedup_summary = expand(config['result_dir'] + "{sample}_dedup_summary.out", sample = sample_set) if config["remove_duplicates"] else [],
	params:
		qsub=config["qsub_default"],
		result_dir = config['result_dir']
	output:
	  config['result_dir'] + "all_dedup_summary.out"
	script: 
		"scripts/dedup_summary.py"
		
### 
# Kaiju ####
###

# select fastq files, select files that exist
def get_fastq_files(sample, readn):
        # cur_file = os.path.isfile("fastq_files/" + sample + "_" + readn + "_nodup.fastq")

        if config["remove_duplicates"]:
                return "fastq_files/" + sample + "_" + readn + "_nodup.fastq"
        else:
                return config["fastq_dir"] + metadata.loc[sample, "read" + readn]

rule kaiju:
	input:
		read1 = lambda wildcards: get_fastq_files(wildcards.sample, "1"),
		read2 = lambda wildcards: get_fastq_files(wildcards.sample, "2")
	params:
		qsub=config["qsub_kaiju"],
		ref=config["reference"],
		rname=config["rname"], 
		kaiju_params=config["kaiju_params"]
	output:
		kaiju_out = config['result_dir'] + "{sample}_kaiju.out"
	shell:
		'''
		kaiju {params.kaiju_params} -t {params.ref}/nodes.dmp -f {params.ref}/kaiju_db_{params.rname}.fmi -i {input.read1} -j {input.read2} -o {output[0]}
		'''

rule kaiju_table:
        input:
                kaiju_out = config['result_dir'] + "{sample}_kaiju.out",
	params:
		qsub=config["qsub_default"],
		ref=config["reference"],
		result_dir = config['result_dir'],
		kaiju_addnames_params=config["kaiju_addnames_params"],
		kaiju_table_params=config["kaiju_table_params"]
	output:
		table = config['result_dir'] + "{sample}_kaiju_table.csv"
	shell:
		'''
		~/.conda/envs/r-base-4.0/bin/Rscript scripts/kaiju_table.R {input.kaiju_out} {output.table}
		'''

rule kaiju_summary:
        input:
                table = expand(config['result_dir'] + "{sample}_kaiju_table.csv", sample = sample_set)
	params:
		qsub=config["qsub_kaiju_summary"]
	output:
		summary = config['result_dir'] + "kaiju_summary.out"
	shell:
		"""
		~/.conda/envs/r-base-4.0/bin/Rscript scripts/kaiju_summary.R "{input.table}" {output.summary}
		"""
