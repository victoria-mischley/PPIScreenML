from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
import json
from pathlib import Path
import argparse
import utils


###For two proteins, must indicate how many chains the second protein contains. This allows for all proteins other than the amount of 
### chains in the second protein to be changed to 'A' and all chains in the second protein to be changed to 'B'
def args():
	parser = argparse.ArgumentParser()
	# Add command line arguments
	parser.add_argument('working_directory', type=str, help='path to location of AF folder output')
	parser.add_argument('num_chains_protein_B', type=str, help='number of chains of the last protein')
	parser.add_argument('interface_cutoff', type=str, help='distance that defines the interface')
	parser.add_argument('csv_name', type=str, help='the name that you would like your csv file to be saved as')


	# Parse the command line arguments
	args = parser.parse_args()

	return args


def main(working_directory, num_chains_protein_B, interface_cutoff):
	directory = Path(working_directory)
	pdb_files  = list(directory.glob('*.pdb'))
	all_json_files = list(directory.glob('*.json'))
	json_files = [file for file in all_json_files if 'error' not in file.name]
	dict_list = []
	for i,file in enumerate(pdb_files):
		results = []
		if file != ".DS_Store":
			file_name = file.stem
			print(file_name)
			pdb_prefix = str(file_name).split('_unrelaxed_')[0] if '_unrelaxed_' in str(file_name) else str(file_name).split('_relaxed_')[0]
			pdb_suffix = str(file_name).split('_unrelaxed_')[1] if '_unrelaxed_' in str(file_name) else str(file_name).split('_relaxed_')[1]
			pdb_suffix_mod = str(pdb_suffix).split("_minimized")[0] 
			json_file_string = next((str(file) for file in json_files if pdb_prefix in str(file) and pdb_suffix in str(file)), None)
			file_path = file
			json_file_path = str(json_file_string)
			AF_chainA_pandas, AF_chainB_pandas = utils.split_chains(file_path, 1)
			AF_chainA = (AF_chainA_pandas.df["ATOM"])
			AF_chainB = (AF_chainB_pandas.df["ATOM"])
	
			### Seperate out the CB atoms or the CA for glycine atoms
			AF_chainA_CB = AF_chainA.loc[(AF_chainA['atom_name'] == 'CB') | ((AF_chainA['residue_name'] == 'GLY') & (AF_chainA['atom_name'] == 'CA'))]
			AF_chainB_CB = AF_chainB.loc[(AF_chainB['atom_name'] == 'CB') | ((AF_chainB['residue_name'] == 'GLY') & (AF_chainB['atom_name'] == 'CA'))]
			### for every atom in A, caluclate the distance from the atom in B and calculate the distance. If distance is less than 8, store in data frame
		
			results_df = utils.get_interface(file_path, num_chains_protein_B, interface_cutoff)
		
			chainA_residue_num = []
			for residue in AF_chainA_CB['residue_name'].index:
				chainA_residue_num.append(AF_chainA['residue_number'][residue])

					###Caclulate plddt:
			all_plDDt = []
			for a in AF_chainA_CB.index:
				b_factor_A = AF_chainA_CB['b_factor'][a]
				all_plDDt.append(b_factor_A)
			for b in AF_chainB_CB.index:
				b_factor_B = AF_chainB_CB['b_factor'][b]
				all_plDDt.append(b_factor_B)
			unique_plDDT = np.unique(all_plDDt)
			final_plddt = np.round(np.average(unique_plDDT), 4)

			f = open(json_file_path, "r")
			json_data = json.loads(f.read())
			PAE = json_data['pae']
			tPAEs = []
			length = (len(chainA_residue_num))
			for resA in AF_chainA_CB['residue_number'].index:
				chainA_res = AF_chainA_CB['residue_number'][resA] -1
				for resB in AF_chainB_CB['residue_number'].index:
					chainB_res = (AF_chainB_CB['residue_number'][resB] -1) + (len(chainA_residue_num)) 
					A_tPAE = PAE[chainA_res]
					B_tPAE = PAE[chainB_res]
					tPAE1 = A_tPAE[chainB_res]
					tPAE2 = B_tPAE[chainA_res]
					tPAEs.append(tPAE1)
					tPAEs.append(tPAE2)
			tPAE_final = np.round((np.average(tPAEs)), 4)
			tPAE_1_4 = int(np.round(((len(tPAEs)) / 4), 0))
			tPAEs_sorted = np.sort(tPAEs)
			tPAEs_top_1_4_sorted = tPAEs_sorted[0:tPAE_1_4]
			tPAEs_top_1_4_final = np.round(np.average(tPAEs_top_1_4_sorted), 4)

			###Calculate ptm
			ptm = json_data['ptm']

		
			if len(results_df) != 0:
				plDDT_A_unique = []
				for val in results_df['plDDt_A'].unique():
					plDDT_A_unique.append(val)
				plDDT_B_unique = []
				for val_B in results_df['plDDt_B'].unique():
					plDDT_B_unique.append(val_B) 

				plDDT_total_values = plDDT_A_unique + plDDT_B_unique
				plDDT_avg = np.average(plDDT_total_values)
				num_of_IF_contacts = len(results)

				###Caclulate iplddt:
				all_iplddt = []
				for x in results_df.index:
					IF_A_plddt = results_df['plDDt_A'][x]
					IF_B_plddt = results_df['plDDt_B'][x]
					all_iplddt.append(IF_A_plddt)
					all_iplddt.append(IF_B_plddt)
				iplddt_unique = np.unique(all_iplddt)
				iplddt_final = np.round(np.average(iplddt_unique), 4)


				###Caclulate iplddt:
				all_iplddt = []
				for x in results_df.index:
					IF_A_plddt = results_df['plDDt_A'][x]
					IF_B_plddt = results_df['plDDt_B'][x]
					all_iplddt.append(IF_A_plddt)
					all_iplddt.append(IF_B_plddt)
				iplddt_final = np.round(np.average(all_iplddt), 2)


				###Calculate iPAE and top_1_4_iPAE:
				### First calculate length of chain A:

				iPAEs = []
				#### Read the json file containing PAE values
				f = open(json_file_path, "r")

				for contact in results_df.index:
					x_fixed = (results_df['res_numberA'][contact]) -1
					y_fixed= (results_df['res_numberB'][contact] -1) + (len(chainA_residue_num)) 
					x_PAE = PAE[x_fixed]
					y_PAE = PAE[y_fixed]
					iPAE1 = x_PAE[y_fixed]
					iPAE2 = y_PAE[x_fixed]
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
				
				
				interface_cutoff_str = str(interface_cutoff)
				ptm_name = f"ptm_{interface_cutoff_str}A"
				plddt_name = f"pLDDT_{interface_cutoff_str}A"
				iplddt_name = f"ipLDDT_{interface_cutoff_str}A"
				iPAE_name = f"iPAE_{interface_cutoff_str}A"
				iPAE_1_4_name = f"iPAE_top_1_4_{interface_cutoff_str}A"
				tPAE_name = f"tPAE_{interface_cutoff_str}A"
				tPAE_1_4_name = f"tPAE_top_1_4_{interface_cutoff_str}A"

				####Export file name and pDockQ score
				data = {'File_name': file_name,
						ptm_name: ptm,
						plddt_name: final_plddt,
						iplddt_name: iplddt_final,
						iPAE_name: iPAE_final,
						iPAE_1_4_name: top_1_4_iPAE_final,
						tPAE_name: tPAE_final,
						tPAE_1_4_name: tPAEs_top_1_4_final}

				dict_list.append(data)

			else:

				interface_cutoff_str = str(interface_cutoff)
				ptm_name = f"ptm_{interface_cutoff_str}A"
				plddt_name = f"pLDDT_{interface_cutoff_str}A"
				iplddt_name = f"ipLDDT_{interface_cutoff_str}A"
				iPAE_name = f"iPAE_{interface_cutoff_str}A"
				iPAE_1_4_name = f"iPAE_top_1_4_{interface_cutoff_str}A"
				tPAE_name = f"tPAE_{interface_cutoff_str}A"
				tPAE_1_4_name = f"tPAE_top_1_4_{interface_cutoff_str}A"

				data = {'File_name': file_name,
						ptm_name: ptm,
						plddt_name: final_plddt,
						iplddt_name: None ,
						iPAE_name: None,
						iPAE_1_4_name:None ,
						tPAE_name: tPAE_final,
						tPAE_1_4_name: tPAEs_top_1_4_final}

				dict_list.append(data)


	alphafold_df = pd.DataFrame.from_dict(dict_list)
	
	return alphafold_df



if __name__ == '__main__':
	args = args()
	features = main(args.working_directory, args.num_chains_protein_B, args.interface_cutoff)
	
	# Assign variables based on command line arguments
	working_directory = args.working_directory
	csv_name = args.csv_name
	# working_directory = "Desktop/PhD/Pipeline/myPPI_Dataset_03_19"
	combined_csv_name = f"{working_directory}/{csv_name}"

	features.to_csv(f"{combined_csv_name}.csv", mode='w', header=True, index=False)
