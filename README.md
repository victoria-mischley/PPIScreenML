# PPIScreenML

Run the get_classification.py script to use the model and get the results. 

Requirements:

- The ColabFold input folder must contain AMBER relaxed models and the corresponding .json file for each model.
- This identifies an interacting protein with respect to the other proteins in the model. You must decide what chains in the pdb file are considered protein 1 and what chains are considered protein 2.
  - For example, if I have a four protein model (chains: A, B, C, D) and want to look at chain B's interaction with the rest of the model, protein 1 wouuld be considered chains A, C, D and protein 2 would be considered chain B.

Software needed to generate models:
- ColabFold
  - link to colabfold github: https://github.com/sokrypton/ColabFold
    
Software needed to generate features:
- pyRosetta
  - link to get pyrosetta academic license: https://els2.comotion.uw.edu/product/pyrosetta
  - link to install pyrosetta: https://www.pyrosetta.org/downloads
  - Instructions for installation of pyrosetta on expanse:
    - create a new conda enviroment: conda create --name newenv
    - activate conda enviroment: conda activate newenv
    - edit ~/.condarc: vi ~/.condarc
    - add in info listed here (use WEST cost): https://www.pyrosetta.org/downloads#h.c0px19b8kvuw
    - run command: conda install pyrosetta

Packages Needed to generate features:
- biopandas
- numpy
- pandas
- argparse
- xgboost
- json
- os
   
- download other python packages into the same enviroment into which pyrosetta was downloaded.

Arguments needed to run the script:
1. working_directory: location of the ColabFold folder on your computer
2. protein1_chains_input: list of the names of the chains that will be considered protein 1
3. protein2_chains_input: list of the names of the chains that will be considered protein 2
4. csv_name: desired name for the output csvfile with the results

Example command:
**you must change into the PPIScreenML directory first***

Command 1: cd <path_to_PPIScreenML_folder>

Command 2: python get_classification.py --working_directory <path_to_folder_with_AF_models> --protein1_chains_input A B --protein2_chains_input C  --csv_name interaction_anlaysis
Outputs:
The feature CSV file. This contains the features extracted for each model. It can be used to looked at individual metrics.
The results CSV file. This contains the scores, the predicted label, and the final prediction of interacting vs. non-interacting.

Optional functionality:
If you are running a high throughput screen with different numbers of proteins, but are consistently modeling the protein for which you are looking for the interaction as the last protein, you can change the flags --protein1_chains_input and --protein2_chains_input to 1. 
- This will all chains other than the last as protein 1 and the last chain as protein 2.
  - For example:
    - Two prtoeins: A, B. Each protein has one chain.
     - This setting will look at the interaction of A to B
  - Two proteins: A, B. Protein A has 2 chains (denoted by ":") and protein B has 1 chain. So in the pdb file, chains A and B are protein 1 and chain C is protein 2
     - This setting will look at A and B against C.
  - Three proteins: A, B, C. Each protein has one chain. You are interested in looking at the interaction of protein A & B with protein C.
     - This setting will look at A and B against C.

Example command:
**you must change into the PPIScreenML directory first***

Command 1: cd <path_to_PPIScreenML_folder>

Command 2: python get_classification.py --working_directory <path_to_folder_with_AF_models> --protein1_chains_input 1 --protein2_chains_input 1  --csv_name interaction_anlaysis

Example command 1: cd PPIScreenML

Example command 2: python get_classification.py ~/Downloads/AFV3_test_models 1 AFV3_test_models_resutls

Outputs:
The feature CSV file. This contains the features extracted for each model. It can be used to looked at individual metrics.
The results CSV file. This contains the scores, the predicted label, and the final prediction of interacting vs. non-interacting.

Other functionality:
- you can use the fix_formatting script alone and this will output the pdb file with the modified chains. For example, if in ColabFold you modeled protein 1 as chains A, B, and C and protein 2 as chain D, this will output a pdb file with chains A, B, and C as A and chain D as B.
- you can run get_interface_dataframe seperatley to get a csv file with the interacting residues of the two proteins
