#!/bin/bash
#BATCH -J af2_job
#SBATCH -o myprefix.%j.output
#SBATCH -e myprefix.%j.err
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH --account=was136
#SBATCH --mem=8G
______________________________________________________________________
conda activate <conda_enviroment>
cd PPIScreenML
python get_classification.py <path_to_models> <number_of_chains_protein_B> <name_for_output_csv_file>
