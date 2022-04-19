#!/bin/bash

#SBATCH --ntasks=150
#SBATCH --time=48:00:00
#SBATCH --job-name=nochange_WF
#SBATCH --output=nochange_WF_out_%j
#SBATCH --error=nochange_WF_error_%j
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=email@school.edu

# By: Ethan Gyllenhaal (egyllenhaal@unm.edu)
# Last updated: 12April2022
# 
# Slurm script for running static Samoa SLiM simulations in parallel
# First initializes output files for nucleotide diversity, population size (not reported in paper), and Fst
# Then runs 100 replicates of SLiM simulation over 8 dispersal parameters
# Slim uses .png rasters per interval generated from tiffs produced by PleistoDist


# Load miniconda and parallel modules
module load miniconda3/4.10.3-an4v
module load parallel/20210922-cfec

# Activate the environment and load env_parallel
source activate slim-env
source $(which env_parallel.bash)

# Change to working directory, make a node list for parallelizing across nodes
cd $SLURM_SUBMIT_DIR
scontrol show hostname > ./node_list_${SLURM_JOB_ID}

# initialize nucleotide diversity output files
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp5.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp10.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp25.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp50.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp100.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp200.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp400.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_PI_disp800.txt

# initialize population size output files
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp5.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp10.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp25.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp50.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp100.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp200.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp400.txt
echo -e 'savaii\tupolu\ttutuila\tofuolo\ttau' > output/samoa_nochange_WF_size_disp800.txt

# initialize FST output files
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\tutuila/ofuolo\ttutuila/tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp5.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila/tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp10.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila/tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp25.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila/tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp50.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila/tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp100.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila/tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp200.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila\tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp400.txt
echo -e 'savaii/upolu\tsavaii/tutuila\tsavaii/ofuolo\tsavaii/tau\tupolu/savaii\tupolu/tutuila\tupolu/ofuolo\tupolu/tau\ttutuila/savaii\ttutuila/upolu\ttutuila/ofuolo\ttutuila\tau\tofuolo/savaii\tofuolo/upolu\tofuolo/tutuila\tofuolo/tau\ttau/savaii\ttau/upolu\ttau/tutuila\ttau/ofuolo' > output/samoa_nochange_WF_FST_disp800.txt

# Parallel call this script is for!
# 1 second delay to lower chance of errors writing to file
# 8 jobs per 8 core node, memory usage is low so isn't a concern, and each job uses a single core
# Runs over the 8 dispersal values 100 times each (1..100 isn't ever used
env_parallel -j 8 --delay 1 --sshloginfile ./node_list_${SLURM_JOB_ID} \
    'slim -d mean_dispersal={1} $SLURM_SUBMIT_DIR/samoa_static.slim' \
    ::: 5 10 25 50 100 200 400 800 ::: {1..100}
