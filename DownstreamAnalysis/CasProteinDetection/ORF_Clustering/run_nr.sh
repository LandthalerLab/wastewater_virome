#!/bin/sh

#$ -N run_nr_cdhit 
#$ -l m_mem_free=10G
#$ -cwd			    
#$ -V 

# Guix
export GUIX_PROFILE=$HOME/.guix-profile
source $GUIX_PROFILE/etc/profile

# Main Job
cd-hit -M 9000 -c 0.8 -d 0 -bak 1 -sc 1 -sf 1 -p 1 -i Trinity_sense_hmmer_filtered.fasta -o results
