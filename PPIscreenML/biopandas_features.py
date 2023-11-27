
from biopandas.pdb import PandasPdb
import numpy as np
from pathlib import Path
import pandas as pd
from collections import Counter
import argparse
import Fix_formatting

def args():
	parser = argparse.ArgumentParser()
	# Add command line arguments
	parser.add_argument('--file_path', type=str, help='path to location of AF folder output')
	parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
	parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
	

	# Parse the command line arguments
	args = parser.parse_args()

	return args

def main(results_df):
	dict_list = []
	AA_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
	polar_uncharged_AAs = ['SER', 'THR', 'TYR', 'ASN', 'GLN']
	nonpolar_AAs = ['GLY', 'ALA', 'VAL', 'CYS', 'PRO', 'LEU', 'ILE', 'MET', 'TRP', 'PHE']
	charged_AAs = ['LYS', 'ARG', 'HIS', 'ASP', 'GLU']
	neg_AAs = ['APS', 'GLU']
	pos_AAs = ['LYS', 'ARG', 'HIS']
	all_polar = ['SER', 'THR', 'TYR', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']


            ###Get number of residues that are involved in the interface
	unique_residue_numbers_A = set()
	unique_residue_numbers_B = set()				
	var_dict = {}
	if len(results_df) != 0:
		for x in results_df.index:
				res_num_A = results_df['res_numberA'][x]
				res_num_B = results_df['res_numberB'][x]
				unique_residue_numbers_A.add(res_num_A)
				unique_residue_numbers_B.add(res_num_B)

		num_unique_residues_at_interface = len(unique_residue_numbers_A) + len(unique_residue_numbers_B)

		results_df_noduplicates_A = results_df.drop_duplicates(subset =['res_numberA'])
		unique_residues_names_A = results_df_noduplicates_A.query('res_numberA in @unique_residue_numbers_A')['resname_A'].tolist()

		results_df_noduplicates_B = results_df.drop_duplicates(subset =['res_numberB'])
		unique_residues_names_B = results_df_noduplicates_B.query('res_numberB in @unique_residue_numbers_B')['resname_B'].tolist()
		unique_residues_names = unique_residues_names_A + unique_residues_names_B
				
		count_dict = {}
				
		for AA in AA_list:
			count = unique_residues_names.count(AA)
			count_dict[AA] = count
			var_name = 'num_{}_IF'.format(AA)
			var_name_edited = f"{var_name}"
			var_dict[var_name_edited] = count

				

		num_of_IF_contacts = len(results_df)
		num_unique_residues_at_interface = num_unique_residues_at_interface
		polar_uncharged_list = [item for item in unique_residues_names if item in polar_uncharged_AAs]
		num_polar_uncharged_AAs_IF = len(polar_uncharged_list)
		nonpolar_list = [item for item in unique_residues_names if item in nonpolar_AAs]
		num_nonpolar_AAs_IF = len(nonpolar_list)
		charged_AAs_list = [item for item in unique_residues_names if item in charged_AAs]
		num_charged_AAs_IF = len(charged_AAs_list)
		neg_AAs_list = [item for item in unique_residues_names if item in neg_AAs]
		num_neg_AAs_IF = len(neg_AAs_list)
		pos_AAs_list = [item for item in unique_residues_names if item in pos_AAs]
		num_pos_AAs_IF = len(pos_AAs_list)
		all_charged_list = [item for item in unique_residues_names if item in all_polar]
		num_all_charged_IF = len(all_charged_list)
		num_all_polar_all_polar_aa_contacts_atIF = len(results_df[results_df['resname_A'].isin(all_polar) & results_df['resname_B'].isin(all_polar)])
		num_polar_polar_aa_contacts_atIF = len(results_df[results_df['resname_A'].isin(polar_uncharged_AAs) & results_df['resname_B'].isin(polar_uncharged_AAs)])
		num_charged_charged_aa_contacts_atIF = len(results_df[results_df['resname_A'].isin(charged_AAs) & results_df['resname_B'].isin(charged_AAs)])
		num_nonpolar_nonpolar_aa_conatacts_atIF = len(results_df[results_df['resname_A'].isin(nonpolar_AAs) & results_df['resname_B'].isin(nonpolar_AAs)])
		num_pos_neg_aa_contact_atIF = len(results_df[results_df['resname_A'].isin(pos_AAs) & results_df['resname_B'].isin(neg_AAs)])
		num_neg_pos_aa_contact_atIF = len(results_df[results_df['resname_A'].isin(neg_AAs) & results_df['resname_B'].isin(pos_AAs)])


		IF_analysis_data = {'num_of_IF_contacts': num_of_IF_contacts,
						'num_unique_residues_at_interface': num_unique_residues_at_interface, 
				        'num_polar_uncharged_AAs_IF': num_polar_uncharged_AAs_IF,
				        'num_nonpolar_AAs_IF': num_nonpolar_AAs_IF,
				        'num_charged_AAs_IF': num_charged_AAs_IF,
				        'num_neg_AAs_IF': num_neg_AAs_IF,
				        'num_pos_AAs_IF': num_pos_AAs_IF,
				        'num_all_charged_IF': num_all_charged_IF,
				        'num_all_polar_all_polar_aa_contacts_atIF': num_all_polar_all_polar_aa_contacts_atIF,
				        'num_polar_polar_aa_contacts_atIF': num_polar_polar_aa_contacts_atIF,
				        'num_charged_charged_aa_contacts_atIF': num_charged_charged_aa_contacts_atIF,
				        'num_nonpolar_nonpolar_aa_conatacts_atIF': num_nonpolar_nonpolar_aa_conatacts_atIF,
				        'num_pos_neg_aa_contact_atIF': num_pos_neg_aa_contact_atIF,
				        'num_neg_pos_aa_contact_atIF': num_neg_pos_aa_contact_atIF
				       }

		IF_analysis_data.update(var_dict)

		dict_list.append(IF_analysis_data)
				
				
				
	else:
		IF_analysis_data = {'num_of_IF_contacts': 0,
					'num_unique_residues_at_interface': 0, 
				    'num_polar_uncharged_AAs_IF': 0,
				    'num_nonpolar_AAs_IF': 0,
				    'num_charged_AAs_IF': 0,
				    'num_neg_AAs_IF': 0,
				    'num_pos_AAs_IF': 0,
				    'num_all_charged_IF': 0,
				    'num_all_polar_all_polar_aa_contacts_atIF': 0,
				    'num_polar_polar_aa_contacts_atIF': 0,
				    'num_charged_charged_aa_contacts_atIF': 0,
				    'num_nonpolar_nonpolar_aa_conatacts_atIF': 0,
				    'num_pos_neg_aa_contact_atIF': 0,
				    'num_neg_pos_aa_contact_atIF': 0
				    }
			

		IF_analysis_data.update(var_dict)
				
		dict_list.append(IF_analysis_data)
		
	biopandas_df = pd.DataFrame.from_dict(dict_list)
	


	return IF_analysis_data		


if __name__ == '__main__':
    args = args()
    file_path = args.file_path
    AF_fixed, AF_chainA, AF_chainB, fixed_json_numbering_dictionary = Fix_formatting.main(file_path, args.protein1_chains_input, args.protein2_chains_input)
    AF_chainA_CB_df = (AF_chainA.df["ATOM"])
    AF_chainB_CB_df = (AF_chainB.df["ATOM"])
    AF_chainA_CB = AF_chainA_CB_df.loc[(AF_chainA_CB_df['atom_name'] == 'CB') | ((AF_chainA_CB_df['residue_name'] == 'GLY') & (AF_chainA_CB_df['atom_name'] == 'CA'))]
    AF_chainB_CB = AF_chainB_CB_df.loc[(AF_chainB_CB_df['atom_name'] == 'CB') | ((AF_chainB_CB_df['residue_name'] == 'GLY') & (AF_chainB_CB_df['atom_name'] == 'CA'))]
    results_df = get_interface_dataframe.main(AF_chainA_CB, AF_chainB_CB)
    file_path = Path(args.file_path)
    working_directory = file_path.parent
    file_name = file_path.name

            
    
    main(results_df)

