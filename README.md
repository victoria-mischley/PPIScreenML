# PPIScreenML

Run the get_classification.py script to use the model and get the results. 

Requirements:

- The ColabFold input folder must contain AMBER relaxed models and the corresponding .json file for each model.
- The protein that you are wishing to check for interactions must be the last protein in the sequence. If it is a 2 protein model, it must be protein 2, if it is a three protein model, it must be protein 3. 

Packages Needed:
- biopandas
- numpy
- pandas
- argparse
- xgboost
- json
- os


Software needed:
- ColabFold
- pyRosetta
  - link to get pyrosetta academic license: https://els2.comotion.uw.edu/product/pyrosetta
  - link to install pyrosetta: https://www.pyrosetta.org/downloads
  - Instructions for installation of pyrosetta on expanse:
    - create a new conda enviroment: conda create --name newenv
    - activate conda enviroment: conda activate newenv
    - edit ~/.condarc: vi ~/.condarc
    - add in info listed here: https://www.pyrosetta.org/downloads#h.c0px19b8kvuw
    - run command: conda install pyrosetta
   
- download other python packages into the same enviroment into which pyrosetta was downloaded.

Arguments needed to run the script:

1. Working directory: location of the ColabFold folder on your computer
2. Number of chains of protein B: this is the number of chains of your last protein. The script will change everything everything but the number of chains of protein B to "A" and the number of chains of protein B to "B". It will then calculate the interface between A and B. Here are examples:
  - Two proteins: A, B. Each protein has one chain.
     - Number of chain protiens = 1
  - Two proteins: A, B. Protein A has 2 chains (denoted by ":") and protein B has 1 chain.
     - Number of interacting protiens = 1
  - Two proteins: A, B. Protein A has 1 chain (denoted by ":") and protein B has 2 chains.
     - Number of interacting protiens = 2
  - Three proteins: A, B, C. Each protein has one chain. You are interested in looking at the interaction of protein A & B with protein C.
     - Number of interacting protiens = 1
  - Three proteins: A, B, C. Each protein has one chain. You are interested in looking at the interaction of protein B & C with protein A.
     - Number of interacting protiens = 2
  - Three proteins: A, B, C. Each protein has one chain. You are interested in looking at the interaction of protein A & C with protein B.
     - Not possible, you would have to reorder the proteins in your fasta file. 3. The name that you want for the output csv files.

**you must change into the PPIScreenML directory first***

Example command: cd PPIScreenML
Example command: python get_classification.py ~/Downloads/AFV3_test_models 1 AFV3_test_models_resutls

Outputs:

The feature CSV file. This contains the features extracted for each model. It can be used to looked at individual metrics.
The results CSV file. This contains the scores, the predicted label, and the final prediction of interacting vs. non-interacting.
