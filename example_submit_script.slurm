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
source ~/.bashrc
conda activate <conda_enviroment>
cd PPIscreenML
python get_classification.py --working_directory <path_to_folder_with_AF_models> --protein1_chains_input A, B --protein2_chains_input C --csv_name interaction_anlaysis
