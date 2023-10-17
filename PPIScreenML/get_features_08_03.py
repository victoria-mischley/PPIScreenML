import biopandas_features_07_13 
import rosetta_features_07_27 
import alphafold_features_07_13 
import pandas as pd
import argparse
########################################################################################################################################
# Create a parser object
# Create a parser object



def args():
	parser = argparse.ArgumentParser(description='Combine dataframes')

	# Add command line arguments
	parser.add_argument('working_directory', type=str, help='path to location of AF folder output')
	parser.add_argument('num_chains_protein_B', type=str, help='number of chains of the last protein')
	parser.add_argument('csv_name', type=str, help='the name that you would like your csv file to be saved as')

	# Parse the command line arguments
	args = parser.parse_args()

	return args


###########################################################################################################################################

def get_feature_csv(working_directory, num_chains_protein_B, interface_cutoff):

	combined_df = get_features(working_directory, num_chains_protein_B, interface_cutoff)


	return combined_df

	
##########################################################################################################################################

def get_features(working_directory, num_chains_protein_B, interface_cutoff):
	
	df_af = alphafold_features_07_13.main(working_directory, num_chains_protein_B, interface_cutoff)
	df_rosetta = rosetta_features_07_27.main(working_directory, num_chains_protein_B, interface_cutoff)
	df_biopandas = biopandas_features_07_13.main(working_directory, num_chains_protein_B, interface_cutoff)

	df_af_ros = pd.merge(df_af, df_rosetta, on='File_name', right_index=False)
	combined_df = pd.merge(df_af_ros, df_biopandas, on='File_name', right_index=False)

	return combined_df


if __name__ == '__main__':
	args = args()
	interface_cutoff = int(12)
	working_directory = args.working_directory
	csv_name = args.csv_name
	num_chains_protein_B = args.num_chains_protein_B
	combined_csv_name = f"{working_directory}/{csv_name}"
	df_all_info = get_feature_csv(args.working_directory, args.num_chains_protein_B, interface_cutoff)
	df_all_info.to_csv(f"{combined_csv_name}.csv", mode='w', header=True, index=False)
	# Assign variables based on command line arguments

