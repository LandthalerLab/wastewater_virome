#!/bin/sh

#$ -N runblast 
#$ -l m_mem_free=200G
#$ -cwd			    
#$ -V 
#$ -t 1-40			    
#$ -j y

# Guix
export GUIX_PROFILE=$HOME/.guix-profile
source $GUIX_PROFILE/etc/profile

# Main Job
file=$(head -$SGE_TASK_ID main_list.txt | tail -1)
~/.conda/envs/rBLAST/bin/Rscript runblast.R $file
