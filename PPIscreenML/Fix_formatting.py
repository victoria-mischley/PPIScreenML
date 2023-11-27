import pandas as pd
import argparse
from pathlib import Path
import numpy as np
from biopandas.pdb import PandasPdb
import copy

def args():
    parser = argparse.ArgumentParser()
    # Add command line arguments
    parser.add_argument('--file_path', type=str, help='path to location of AF folder output')
    parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
    parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
    # parser.add_argument('num_chains_protein_B', type=str, help='number of chains of the last protein')
    # Parse the command line arguments
    args = parser.parse_args()
    return args

def main(file_path, protein1_chains_input, protein2_chains_input):
    #####Read in file
    AF = PandasPdb().read_pdb(str(file_path))
    print(file_path)

    ###this allows one to use 1 to take everything but the last chain as protein 1 and the last chain as protein2
    ### Allows to run in one run proteins with different number of chains, as long as the last protein is only one chain
    meet_criteria = ['1']
    chains = []
    if protein1_chains_input and protein2_chains_input == meet_criteria:
        for chain_ID in AF.df['ATOM']['chain_id'].unique():
            chains.append(chain_ID)
        protein1_chains = chains[:-1]
        protein2_chains = chains[-1:]

    else:
        protein1_chains = protein1_chains_input
        protein2_chains = protein2_chains_input
    


    chains_old_order = []
    for chain_ID in AF.df['ATOM']['chain_id'].unique():
        chains_old_order.append(chain_ID)

    chains_new_order = protein1_chains + protein2_chains
    
    counter7 = 0
    old_chain_index_dict = {}
    for val in chains_old_order:
        old_chain_index_dict[counter7] = val
        counter7 += 1

    counter8 = 0
    new_chain_index_dict = {}
    for val in chains_new_order:
        new_chain_index_dict[counter8] = val
        counter8 += 1

    
    ###get old atom_number_residue_pairing
    old_atomnum_resnumb_dict = {}
    for ind in AF.df['ATOM']['atom_number'].index:
        atom_number = AF.df['ATOM']['atom_number'][ind]
        residue_number = AF.df['ATOM']['residue_number'][ind]
        old_atomnum_resnumb_dict[atom_number] = residue_number
    
    ###get dict list for all chains in protein old atom_num and res_num
    old_resnum_atom_num_dict_list = []
    old_unique_res_set_list = []
    chain_length_dict_list = []
    chain_length_dict = {}
    for val in chains_old_order:
        dictionary = {}
        residue_numbers = set()
        chain_length_temp = {}
        AF_chain = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == val]
        for ind in AF_chain.index:
            atom_number = AF_chain['atom_number'][ind]
            residue_number = AF_chain['residue_number'][ind]
            dictionary[atom_number] = residue_number
            residue_numbers.add(residue_number)
        old_resnum_atom_num_dict_list.append(dictionary)
        old_unique_res_set_list.append(residue_numbers)

        chain_length = len(residue_numbers)
        chain_length_temp[val] = chain_length
        chain_length_dict[val] = chain_length
        chain_length_dict_list.append(chain_length_temp)


    ##Get old atom_number_residue_pairing for protein 1
    protein1_old_resnum_atom_num_dict_list = []
    protein1_old_unique_res_set_list = []
    protein_1_old_atom_numbers = set()
    
    for val in protein1_chains:
        dictionary = {}
        residue_numbers = set()
        AF_chain = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == val]
        for ind in AF_chain.index:
            atom_number = AF_chain['atom_number'][ind]
            residue_number = AF_chain['residue_number'][ind]
            dictionary[atom_number] = residue_number
            residue_numbers.add(residue_number)
            protein_1_old_atom_numbers.add(atom_number)
        protein1_old_resnum_atom_num_dict_list.append(dictionary)
        protein1_old_unique_res_set_list.append(residue_numbers)
    
    ##Get old atom_number_residue_pairing for protein 2
    protein2_old_resnum_atom_num_dict_list = []
    protein2_old_unique_res_set_list = []
    protein_2_old_atom_numbers = set()
    for val in protein2_chains:
        dictionary = {}
        residue_numbers = set()
        AF_chain = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == val]
        for ind in AF_chain.index:
            atom_number = AF_chain['atom_number'][ind]
            residue_number = AF_chain['residue_number'][ind]
            dictionary[atom_number] = residue_number
            residue_numbers.add(residue_number)
            protein_2_old_atom_numbers.add(atom_number)
        protein2_old_resnum_atom_num_dict_list.append(dictionary)
        protein2_old_unique_res_set_list.append(residue_numbers)
    
    ###get new residue numbers protein 1
    offset = 0
    counter5 = 0
    protein1_new_atomnum_res_num_dict = {}
    for val in protein1_old_resnum_atom_num_dict_list:
        chain_dict = val
        counter10 = 1
        current_res_number = None  # Initialize to None to handle the first iteration
        for key, value in chain_dict.items():
            atom_number = key
            old_residue_number = value

        # Check if the residue number has changed
            if current_res_number is not None and current_res_number != old_residue_number:
                counter10 += 1
            new_residue_number = counter10 + offset
            protein1_new_atomnum_res_num_dict[atom_number] = new_residue_number
        # Update the old_residue_number for the next iteration
            current_res_number = old_residue_number
        offset += counter10
        counter5 +=1


    ###get new residue numbers protein 2
    offset = 0
    counter6 = 0
    protein2_new_atomnum_res_num_dict = {}
    for val in protein2_old_resnum_atom_num_dict_list:
        chain_dict = val
        counter11 = 1
        current_res_number = None  # Initialize to None to handle the first iteration
        for key, value in chain_dict.items():
            atom_number = key
            old_residue_number = value
        # Check if the residue number has changed
            if current_res_number is not None and current_res_number != old_residue_number:
                counter11 += 1
            new_residue_number = counter11 + offset
            protein2_new_atomnum_res_num_dict[atom_number] = new_residue_number
        # Update the old_residue_number for the next iteration
            current_res_number = old_residue_number
        offset += counter11
        counter6 +=1

    combined_new_atomnum_resnum = {**protein1_new_atomnum_res_num_dict, **protein2_new_atomnum_res_num_dict}

    ###get numbering for json file.
    old_atom_num_json_num = {}
    counter15 = 0
    current_res_number = None
    old_atom_number_list = []
    for ind in AF.df['ATOM']['atom_number'].index:
        old_atom_number = AF.df['ATOM']['atom_number'][ind]
        old_atom_number_list.append(old_atom_number)
        old_residue_number = AF.df['ATOM']['residue_number'][ind]
        json_number = counter15
        old_atom_num_json_num[old_atom_number] = json_number
        if current_res_number is not None and current_res_number != old_residue_number:
            counter15 +=1
        current_res_number = old_residue_number

    
    ###change the numbering###

    for ind in AF.df['ATOM']['atom_number'].index:
        atom_number = AF.df['ATOM']['atom_number'][ind]
        old_residue_number = AF.df['ATOM']['residue_number']
        new_residue_number = combined_new_atomnum_resnum[atom_number]
        AF.df['ATOM'].at[ind, 'residue_number'] = new_residue_number
    
    ####Change the  chain IDS###
    for ind in AF.df['ATOM']['chain_id'].index:
        current_chain_ID = AF.df['ATOM']['chain_id'][ind]
        if current_chain_ID in protein1_chains:
            AF.df['ATOM'].at[ind, 'chain_id'] = 'A'
        else:
            AF.df['ATOM'].at[ind, 'chain_id'] = 'B'

    ###change the atom_number
    counter12 = 1
    counter13 = 1
    for ind in AF.df['ATOM']['chain_id'].index:
        current_chain_ID = AF.df['ATOM']['chain_id'][ind]
        if current_chain_ID == 'A':
            AF.df['ATOM'].at[ind, 'atom_number'] = counter12
            counter12 += 1

    for ind in AF.df['ATOM']['chain_id'].index:
        current_chain_ID = AF.df['ATOM']['chain_id'][ind]
        if current_chain_ID == 'B':
            new_atom = counter12 + counter13 -1
            AF.df['ATOM'].at[ind, 'atom_number'] = new_atom
            counter13 += 1

  
    old_atom_df = AF.df['ATOM']
    sorted_atom_df = AF.df['ATOM'].sort_values(by='atom_number').reset_index(drop=True)
    

    # Update the PandasPdb object with the sorted DataFrame
    AF.df['ATOM'] = sorted_atom_df

    ##Seperate into chainA dataframe and chainB dataframe
    AF_chainA = copy.deepcopy(AF)
    AF_chainB = copy.deepcopy(AF)

    AF_chainA.df['ATOM'] = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'A']
    AF_chainB.df['ATOM'] = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'B']

    new_atom_new_res = {}
    new_atom_number_list = []
    for ind in sorted_atom_df['atom_number'].index:
        new_atom_number = sorted_atom_df['atom_number'][ind]
        new_residue_number = sorted_atom_df['residue_number'][ind]
        new_atom_new_res[new_atom_number] = new_residue_number
        new_atom_number_list.append(new_atom_number)

    old_atom_num_new_atom_num = dict(zip(old_atom_number_list, new_atom_number_list))

    new_atom_num_json_num = {}
    for key, value in combined_new_atomnum_resnum.items():
        old_atom_num = key
        json_number = old_atom_num_json_num[old_atom_num]
        new_atom_num = old_atom_num_new_atom_num[old_atom_num]
        new_atom_num_json_num[new_atom_num] = json_number
    

    return AF, AF_chainA, AF_chainB, new_atom_num_json_num

if __name__ == '__main__':
    args = args()
    AF, AF_chainA, AF_chainB, new_atom_num_json_num = main(args.file_path, args.protein1_chains_input, args.protein2_chains_input)
    file_path = Path(args.file_path)
    file_name = file_path.name[:-4]
    working_directory = file_path.parent
    fixed_PDB_output_file_path = f"{working_directory}/{file_name}_fixed.pdb"
    AF.to_pdb(path=fixed_PDB_output_file_path, records=None, gz=False, append_newline=True)
