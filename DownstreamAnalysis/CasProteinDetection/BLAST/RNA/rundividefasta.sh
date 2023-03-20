#!/bin/sh

#$ -N rundividefasta 
#$ -o rundividefasta.out 
#$ -e rundividefasta.err
#$ -cwd			    
#$ -V 			    

# Guix
export GUIX_PROFILE=$HOME/.guix-profile
source $GUIX_PROFILE/etc/profile

# Main Job
~/.conda/envs/crisprcastyper/bin/pyfasta split -n 100 fasta_files/Trinity_RNA_sense_proteins_filtered_cas.ffa
