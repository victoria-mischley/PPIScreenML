from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
import json
from pathlib import Path
import argparse
import get_interface_dataframe
import Fix_formatting


###For two proteins, must indicate how many chains the second protein contains. This allows for all proteins other than the amount of 
### chains in the second protein to be changed to 'A' and all chains in the second protein to be changed to 'B'
def args():
	parser = argparse.ArgumentParser(description='Combine dataframes')

	# Add command line arguments
	parser.add_argument('--file_path', type=str, help='path to location of AF folder output')
	parser.add_argument('--json_file_path', type=str, help='path to location of AF folder output')
	parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
	parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
	
    # Parse the command line arguments
	args = parser.parse_args()
	return args


def main(results_df,AF_chainA_CB, AF_chainB_CB, json_file_path, fixed_json_numbering_dictionary):
	dict_list =[]
	###Caclulate plddt
	all_plDDt = []
	for a in AF_chainA_CB.index:
		b_factor_A = AF_chainA_CB['b_factor'][a]
		all_plDDt.append(b_factor_A)
	for b in AF_chainB_CB.index:
		b_factor_B = AF_chainB_CB['b_factor'][b]
		all_plDDt.append(b_factor_B)
	unique_plDDT = np.unique(all_plDDt)
	final_plddt = np.round(np.average(unique_plDDT), 4)

	####Load the .json file######
	f = open(json_file_path, "r")
	json_data = json.loads(f.read())
	PAE = json_data['pae']

			###Calculate ptm
	ptm = json_data['ptm']

	if len(results_df) != 0:
		###Caclulate iplddt:
		all_iplddt = []
		for x in results_df.index:
			IF_A_plddt = results_df['plDDt_A'][x]
			IF_B_plddt = results_df['plDDt_B'][x]
			all_iplddt.append(IF_A_plddt)
			all_iplddt.append(IF_B_plddt)
			iplddt_final = np.round(np.average(all_iplddt), 2)
		

		#####get old numbering and new numbering for residue numbers###
		tPAEs = []
	
		for index, row in AF_chainA_CB.iterrows():
			protein1_res = row['residue_number']
			protein1_atom_number = row['atom_number']
			protein1_json_number = fixed_json_numbering_dictionary[protein1_atom_number]
			for index2, row2 in AF_chainB_CB.iterrows():
						protein2_res = row2['residue_number']
						protein2_atom_number = row2['atom_number']
						protein2_json_number = fixed_json_numbering_dictionary[protein2_atom_number]
						protein1_tPAE = PAE[protein1_json_number]
						protein2_tPAE = PAE[protein2_json_number]
						tPAE1 = protein1_tPAE[protein2_json_number]
						tPAE2 = protein2_tPAE[protein1_json_number]
						tPAEs.append(tPAE1)
						tPAEs.append(tPAE2)
		tPAE_final = np.round((np.average(tPAEs)), 4)
		tPAE_1_4 = int(np.round(((len(tPAEs)) / 4), 0))
		tPAEs_sorted = np.sort(tPAEs)
		tPAEs_top_1_4_sorted = tPAEs_sorted[0:tPAE_1_4]
		tPAEs_top_1_4_final = np.round(np.average(tPAEs_top_1_4_sorted), 4)

				###Calculate iPAE and top_1_4_iPAE:
				### First calculate length of chain A:

		iPAEs = []
		json_number_check = []
		PAE_check = []
				#### Read the json file containing PAE value

		for index,row3 in results_df.iterrows():
			protein1_atom_number = row3['atom_numberA']
			protein2_atom_number =  row3['atom_numberB']
			protein1_json_number = fixed_json_numbering_dictionary[protein1_atom_number]
			json_number_check.append(protein1_json_number)
			protein2_json_number = fixed_json_numbering_dictionary[protein2_atom_number]
			json_number_check.append(protein2_json_number)
			protein1_PAE = PAE[protein1_json_number]
			PAE_check.append(protein1_PAE[0])
			protein2_PAE = PAE[protein2_json_number]
			PAE_check.append(protein2_PAE[0])
			iPAE1 = protein1_PAE[protein2_json_number]
			iPAE2 = protein2_PAE[protein1_json_number]
			iPAEs.append(iPAE1)
			iPAEs.append(iPAE2)
		iPAE_final = np.round((np.average(iPAEs)), 4)
		iPAEs_sorted = np.sort(iPAEs)
		if len(iPAEs) >= 8 :
			if_1_4 = int(np.round(((len(iPAEs)) / 4), 0))
			top_1_4_iPAE_sorted = iPAEs_sorted[0:(if_1_4 -1)] 
			top_1_4_iPAE_final = np.round(np.average(top_1_4_iPAE_sorted), 4)
		if len(iPAEs) < 8:
			top_1_4_iPAE_sorted = np.array([iPAEs_sorted])
			top_1_4_iPAE_final = np.round(np.average(top_1_4_iPAE_sorted), 4)


				####Export file name and pDockQ score
		data = {'ptm': ptm,
				'plddt': final_plddt,
				'iplddt': iplddt_final,
				'iPAE': iPAE_final,
				'iPAE_1_4': top_1_4_iPAE_final,
				'tPAE': tPAE_final,
				'tPAE_1_4': tPAEs_top_1_4_final}

		dict_list.append(data)

	else:
		data = {'ptm': ptm,
				'plddt': final_plddt,
				'iplddt': None ,
				'iPAE': None,
				'iPAE_1_4':None ,
				'tPAE': tPAE_final,
				'tPAE_1_4': tPAEs_top_1_4_final}

		dict_list.append(data)


	alphafold_df = pd.DataFrame.from_dict(dict_list)
	
	return data



if __name__ == '__main__':
	args = args()
	file_path = args.file_path
	AF_fixed, AF_chainA, AF_chainB, fixed_json_numbering_dictionary = Fix_formatting.main(file_path, args.protein1_chains_input, args.protein2_chains_input)
	AF_chainA_CB_df = (AF_chainA.df["ATOM"])
	AF_chainB_CB_df = (AF_chainB.df["ATOM"])
	AF_chainA_CB = AF_chainA_CB_df.loc[(AF_chainA_CB_df['atom_name'] == 'CB') | ((AF_chainA_CB_df['residue_name'] == 'GLY') & (AF_chainA_CB_df['atom_name'] == 'CA'))]
	AF_chainB_CB = AF_chainB_CB_df.loc[(AF_chainB_CB_df['atom_name'] == 'CB') | ((AF_chainB_CB_df['residue_name'] == 'GLY') & (AF_chainB_CB_df['atom_name'] == 'CA'))]
	results_df = get_interface_dataframe.main(AF_chainA_CB, AF_chainB_CB)
	main(results_df,AF_chainA_CB, AF_chainB_CB, args.json_file_path, fixed_json_numbering_dictionary)

