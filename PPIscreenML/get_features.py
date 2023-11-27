import biopandas_features
import rosetta_features
import alphafold_features
import pandas as pd
import argparse
from pathlib import Path
import Fix_formatting
import os
import get_interface_dataframe
########################################################################################################################################
# Create a parser object
# Create a parser object



def args():
	parser = argparse.ArgumentParser(description='Combine dataframes')

	# Add command line arguments
	parser.add_argument('--working_directory', type=str, help='path to location of AF folder output')
	parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
	parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
	parser.add_argument('--csv_name', type=str, help='the name that you would like your csv file to be saved as')

	# Parse the command line arguments
	args = parser.parse_args()

	return args


###########################################################################################################################################

def main(working_directory, protein1_chains_input, protein2_chains_input, csv_name):
	directory = Path(working_directory)
	pdb_files  = list(directory.glob('*.pdb'))
	all_json_files = list(directory.glob('*.json'))
	json_files = [file for file in all_json_files if 'error' not in file.name]
	dict_list = []
	final_df = pd.DataFrame()
	for i,file in enumerate(pdb_files):
		results = []
		if file != ".DS_Store":
			##Match pdb name with appropriate .json file#####
			file_name = file.stem
			print(file_name)
			pdb_prefix = str(file_name).split('_unrelaxed_')[0] if '_unrelaxed_' in str(file_name) else str(file_name).split('_relaxed_')[0]
			pdb_suffix = str(file_name).split('_unrelaxed_')[1] if '_unrelaxed_' in str(file_name) else str(file_name).split('_relaxed_')[1]
			json_file_string = next((str(file) for file in json_files if pdb_prefix in str(file) and pdb_suffix in str(file)), None)
			file_path = file
			json_file_path = str(json_file_string)

			###Fix formating based on chain IDs for protein 1 and protein 2
			AF_fixed, AF_chainA, AF_chainB, fixed_json_numbering_dictionary = Fix_formatting.main(file_path, protein1_chains_input, protein2_chains_input)

			AF_chainA_CB_df = (AF_chainA.df["ATOM"])
			AF_chainB_CB_df = (AF_chainB.df["ATOM"])

			AF_chainA_CB = AF_chainA_CB_df.loc[(AF_chainA_CB_df['atom_name'] == 'CB') | ((AF_chainA_CB_df['residue_name'] == 'GLY') & (AF_chainA_CB_df['atom_name'] == 'CA'))]
			AF_chainB_CB = AF_chainB_CB_df.loc[(AF_chainB_CB_df['atom_name'] == 'CB') | ((AF_chainB_CB_df['residue_name'] == 'GLY') & (AF_chainB_CB_df['atom_name'] == 'CA'))]

			######get interface df###
			results_df = get_interface_dataframe.main(AF_chainA_CB, AF_chainB_CB)

			####Get features
			alphafold_features_dict = alphafold_features.main(results_df,AF_chainA_CB, AF_chainB_CB, json_file_path, fixed_json_numbering_dictionary)
			rosetta_features_dict = rosetta_features.main(working_directory, file_name, results_df, AF_fixed, AF_chainA, AF_chainB)
			biopandas_features_dict = biopandas_features.main(results_df)

			###combine features for a single file
			data = {'File_name': file_name}
			data.update(alphafold_features_dict)
			data.update(rosetta_features_dict)
			data.update(biopandas_features_dict)

####Append to dict_list to get features from all files in the folder
			dict_list.append(data)

######Export the final features df with features from all files in the folder
	final_df = pd.DataFrame.from_dict(dict_list)
	final_df.reset_index(drop = True, inplace = True)
	
	csv_file_path = f"{working_directory}/{csv_name}_features.csv"
	final_df.to_csv(csv_file_path)

	return final_df

			



			




if __name__ == '__main__':
	args = args()
	main(args.working_directory, args.protein1_chains_input, args.protein2_chains_input, args.csv_name)

	# Assign variables based on command line arguments

