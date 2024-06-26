#!/bin/bash
#SBATCH --job-name=sbs_1-hdf
#SBATCH --partition=20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --time=8:00:00
#SBATCH --output sbs_1-hdf-%j.out

# Load any necessary modules or activate virtual environments
source /lab/barcheese01/mdiberna/OpticalPooledScreens/venv_ops_new/bin/activate

# Change to the directory where your Snakefiles are located
cd /lab/barcheese01/screens/baker

# Execute Snakemake
python3 sbs_1-hdf.py