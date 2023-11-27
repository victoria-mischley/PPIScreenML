import pandas as pd
import argparse
from pathlib import Path
import Fix_formatting
import os
import numpy as np

########################################################################################################################################
# Create a parser object
# Create a parser object


def args():
	parser = argparse.ArgumentParser(description='Combine dataframes')

	# Add command line arguments
	parser.add_argument('--file_path', type=str, help='path to location of AF folder output')
	parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
	parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
	
    # Parse the command line arguments
	args = parser.parse_args()
	return args
def main(AF_chainA_CB, AF_chainB_CB):
	interface_cutoff = 12
	results = []
	for index, row in AF_chainA_CB.iterrows():
		X_coord_1 = row['x_coord']
		Y_coord_1 = row['y_coord']
		Z_coord_1 = row['z_coord']
		plDDt_1 = row['b_factor']
		resname_1 = row['residue_name']
		res_number_1 = row['residue_number']
		atom_number_1 = row['atom_number']
		for index2, row2 in AF_chainB_CB.iterrows():
			X_coord_2 = row2['x_coord']
			Y_coord_2 = row2['y_coord']
			Z_coord_2 = row2['z_coord']
			plDDt_2 = row2['b_factor']
			resname_2 = row2['residue_name']
			res_number_2 = row2['residue_number']
			atom_number_2 = row2['atom_number']
			distance = np.sqrt((X_coord_1 - X_coord_2) **2 + (Y_coord_1 - Y_coord_2) **2 + (Z_coord_1 - Z_coord_2) **2)
			if distance <= int(interface_cutoff):
				results.append({'resname_A': resname_1, 'res_numberA': res_number_1, 'atom_numberA':atom_number_1, 'resname_B': resname_2, 'res_numberB': res_number_2, 'atom_numberB':atom_number_2, 'plDDt_A': plDDt_1, 'plDDt_B': plDDt_2, 'distance': distance})
	results_df = pd.DataFrame(results)
	return results_df
			

if __name__ == '__main__':
	args = args()
	AF_fixed, AF_chainA, AF_chainB, new_atom_num_json_num = Fix_formatting.main(args.file_path, args.protein1_chains_input, args.protein2_chains_input)
	AF_fixed_df = AF_fixed.df["ATOM"]
	AF_chainA_CB_df = (AF_chainA.df["ATOM"])
	AF_chainB_CB_df = (AF_chainB.df["ATOM"])
	AF_chainA_CB = AF_chainA_CB_df.loc[(AF_chainA_CB_df['atom_name'] == 'CB') | ((AF_chainA_CB_df['residue_name'] == 'GLY') & (AF_chainA_CB_df['atom_name'] == 'CA'))]
	AF_chainB_CB = AF_chainB_CB_df.loc[(AF_chainB_CB_df['atom_name'] == 'CB') | ((AF_chainB_CB_df['residue_name'] == 'GLY') & (AF_chainB_CB_df['atom_name'] == 'CA'))]
	AF_fixed_CB = AF_fixed_df.loc[(AF_fixed_df['atom_name'] == 'CB') | ((AF_fixed_df['residue_name'] == 'GLY') & (AF_fixed_df['atom_name'] == 'CA'))]
	main(AF_chainA_CB, AF_chainB_CB)
	
	file_path = Path(args.file_path)
	working_directory = file_path.parent
	csv_file_path = f"{working_directory}/interface_dataframe.csv"
	results_df.to_csv(csv_file_path, index=False)
