import pandas as pd
import os

# dataset and samples
metadata = pd.read_table(config["metadata"], dtype={"sample": str}, sep = " ").set_index(["sample"], drop=False)
sample_set = metadata["sample"].tolist()

# parameters
number_of_folds_pyfasta = int(config['number_of_folds_pyfasta'])
if number_of_folds_pyfasta > 9:
  pyfasta_fold_indices = ["%02d" % i for i in range(0,int(config['number_of_folds_pyfasta']))]
else:
  pyfasta_fold_indices = [i for i in range(0,int(config['number_of_folds_pyfasta']))]
  
  
# target rule
rule all:
  input:
    cctyper_summary = expand(config["result_dir"] + "{sample}_hmmer.tab", sample = sample_set)
  params:
    qsub = config['qsub.default']

rule divide_fasta:
        input:
                read = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read"]
        params:
                qsub = config['qsub.default'],
                nfolds = config['number_of_folds_pyfasta']
        output:
                temp(expand("fasta_files/{{sample}}.{fold}.fasta", fold = pyfasta_fold_indices))
        shell:
                '''
                mkdir -p fasta_files
                ln -s {input.read} fasta_files/{wildcards.sample}.fasta
                ~/.conda/envs/crisprcastyper/bin/pyfasta split -n {params.nfolds} fasta_files/{wildcards.sample}.fasta
        
                # rm intermediate files
                rm fasta_files/{wildcards.sample}.fasta.flat
                rm fasta_files/{wildcards.sample}.fasta.gdx
                '''

rule run_cctyper:
        input:
                read = "fasta_files/{sample}.{fold}.fasta"
        params:
                qsub = config['qsub.cctyper'],
                result_dir = config["result_dir"],
                reference_hmmer = config["reference_hmmer"],
                cctyper_params = config["cctyper_params"]
        output:
                sample_fold_hmmer = temp(config["result_dir"] + "{sample}_fold{fold}_hmmer.tab")
        shell:
                '''
                export PATH="$HOME/.conda/envs/crisprcastyper/bin/:$PATH"
                mkdir -p {params.result_dir}
                export CCTYPER_DB="{params.reference_hmmer}"
                cctyper {input.read} {params.result_dir}/{wildcards.sample}_{wildcards.fold} {params.cctyper_params}   
                cp {params.result_dir}/{wildcards.sample}_{wildcards.fold}/hmmer.tab {output.sample_fold_hmmer}
                '''

rule aggregate_run_cctyper:
        input:
                sample_fold_hmmer = expand(config["result_dir"] + "{{sample}}_fold{fold}_hmmer.tab", fold = pyfasta_fold_indices)
        params:
                qsub = config['qsub.default'],
        output:
                sample_hmmer = config["result_dir"] + "{sample}_hmmer.tab"
        shell:
                '''
                IFS=', ' read -r -a array <<< "{input.sample_fold_hmmer}"
                head -1 ${{array[0]}} > {output.sample_hmmer} 
                for fold in "${{array[@]}}"
                do
                '''
