#!/bin/bash

#SBATCH --ntasks=100
#SBATCH --time=48:00:00
#SBATCH --job-name=Fiji_msp
#SBATCH --output=Fiji_out_%j
#SBATCH --error=Fiji_error_%j
#SBATCH --partition=debug
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=email@school.edu

# By: Ethan Gyllenhaal (egyllenhaal@unm.edu)
# Last updated: 12April2022
# 
# Slurm script for running static and dynamic Fiji msprime simulations in parallel
# First runs dynamic, then static, for 100 replicates with each parameter value
# param_list is comma delimited file with one line per combo, then format alpha,mean_dispersal
# alpha is number of propagules per m^2, and mean_dispersal is mean dispersal distance in meters
# Setting alpha to zero means no gene flow, we used a value of 0.0000000028023 for the paper
# Also takes in input csv with distance, target sizes, and island sizes per interval from PleistoDist
# Header is: Interval,Distance,Kadavu_target,Viti_target,Kadavu_size,Viti_size

# Load miniconda and parallel modules
module load miniconda3/4.10.3-an4v
module load parallel/20210922-cfec

# Activate the environment and load env_parallel
source activate pleisto_msp-env
source $(which env_parallel.bash)

# Set a directory path to working directory, then make node list for parallelizing across nodes
dir=$SLURM_SUBMIT_DIR
scontrol show hostname > $dir/node_list_${SLURM_JOB_ID}

# While loop reads in parameter list to get alpha and dispersal values, fed into unix variables
while IFS=, read alpha disp; do
	# Make the output file for that parameter pair
    echo -e 'Fst\tViti div\t Kadavu div' > $dir/output/fiji_dynamic_meandisp${disp}.out
    # Run parallel for parameter pair for dynamic, 8 jobes per 8 core node, 1 second delay to avoid output errors
	# Note this maxes out at running as many tasks at a time as you have replicates
	# If you want to run more tasks at a time than you have replicates, you'll have to modify this
	env_parallel -j 8 --delay 1 --sshloginfile $dir/node_list_${SLURM_JOB_ID} \
	'echo {}
         python $dir/msp_fiji_fluctuate.py -i $dir/msprime_input.csv \
                -o $dir/output/fiji_dynamic_meandisp${disp}.out -a $alpha -d $disp' ::: {1..100}
done < $dir/param_list

# repeat for static, otherwise the same as above
while IFS=, read alpha disp; do
    echo -e 'Fst\tViti div\t Kadavu div' > $dir/output/fiji_static_meandisp${disp}.out
    env_parallel -j 8 --delay 1 --sshloginfile $dir/node_list_${SLURM_JOB_ID} \
        'echo {}
         python $dir/msp_fiji_static.py -i $dir/msprime_input.csv \
                -o $dir/output/fiji_static_meandisp${disp}.out -a $alpha -d $disp' ::: {1..100}
done < $dir/param_list
