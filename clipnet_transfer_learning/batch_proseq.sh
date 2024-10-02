#!/bin/bash -l
#SBATCH --job-name=k562_proseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tP100:1
#SBATCH --mail-type=all
#SBATCH --mail-user=ayh8@cornell.edu
#SBATCH --array=5-9

# SLURM cript to run the training of the model on the CBSU GPU node

#mamba init
conda activate clipnet
cd /home2/ayh8/clipnet_k562/clipnet_transfer_learning
time python transfer_learn_k562_proseq.py $SLURM_ARRAY_TASK_ID 0