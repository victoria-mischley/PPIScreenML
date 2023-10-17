import pandas as pd
import argparse
from pathlib import Path
import numpy as np
from biopandas.pdb import PandasPdb
import copy

def args():
    parser = argparse.ArgumentParser()
    # Add command line arguments
    parser.add_argument('file_path', type=str, help='path to location of AF folder output')
    parser.add_argument('num_chains_protein_B', type=str, help='number of chains of the last protein')
    parser.add_argument('interface_cutoff', type=str, help='distance that defines the interface')
    # Parse the command line arguments
    args = parser.parse_args()
    return args

def get_interface(file_path, num_chains_protein_B, interface_cutoff):
    results = []
    AF = PandasPdb().read_pdb(str(file_path))
    chains_AF = []
    for chain_ID in AF.df['ATOM']['chain_id'].unique():
        chains_AF.append(chain_ID)
    
    total_number_chains = len(chains_AF)
    num_chains_protein_A = total_number_chains - int(num_chains_protein_B)
    chains_AF_A = chains_AF[:num_chains_protein_A]
    chains_AF_B = chains_AF[num_chains_protein_A:]
    
    for ind in AF.df['ATOM']['chain_id'].index:
        current_chain_ID = AF.df['ATOM']['chain_id'][ind]
        if current_chain_ID in chains_AF_A:
            AF.df['ATOM'].at[ind, 'chain_id'] = 'A'
        else:
            AF.df['ATOM'].at[ind, 'chain_id'] = 'B'
        
    AF_chainA = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'A']
    AF_chainB = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'B']
    
    AF_chainA_CB = AF_chainA.loc[(AF_chainA['atom_name'] == 'CB') | ((AF_chainA['residue_name'] == 'GLY') & (AF_chainA['atom_name'] == 'CA'))]
    AF_chainB_CB = AF_chainB.loc[(AF_chainB['atom_name'] == 'CB') | ((AF_chainB['residue_name'] == 'GLY') & (AF_chainB['atom_name'] == 'CA'))]
    
    for ind in AF_chainA_CB.index:
        X_coord_A = AF_chainA_CB['x_coord'][ind]
        Y_coord_A = AF_chainA_CB['y_coord'][ind]
        Z_coord_A = AF_chainA_CB['z_coord'][ind]
        plDDt_A = AF_chainA_CB['b_factor'][ind]
        resname_A = AF_chainA_CB['residue_name'][ind]
        res_number_A = AF_chainA_CB['residue_number'][ind]
        for B_ind in AF_chainB_CB.index:
            X_coord_B = AF_chainB_CB['x_coord'][B_ind]
            Y_coord_B = AF_chainB_CB['y_coord'][B_ind]
            Z_coord_B = AF_chainB_CB['z_coord'][B_ind]
            plDDt_B = AF_chainB_CB['b_factor'][B_ind]
            resname_B = AF_chainB_CB['residue_name'][B_ind]
            res_number_B = AF_chainB_CB['residue_number'][B_ind]
            distance = np.sqrt((X_coord_A - X_coord_B) **2 + (Y_coord_A - Y_coord_B) **2 + (Z_coord_A - Z_coord_B) **2)
            if distance <= int(interface_cutoff):
                results.append({'resname_A': resname_A, 'res_numberA': res_number_A, 'resname_B': resname_B, 'res_numberB': res_number_B, 'plDDt_A': plDDt_A, 'plDDt_B': plDDt_B, 'distance': distance})
    results_df = pd.DataFrame(results)
    
    return results_df

def split_chains(file_path, num_chains_protein_B):
    AF = PandasPdb().read_pdb(str(file_path))

    chains = AF.df["ATOM"]["chain_id"].unique()
    
    chain_A = chains[:-int(num_chains_protein_B)]
    chain_B = chains[-int(num_chains_protein_B):]
    
    AF_chainA = copy.deepcopy(AF)
    AF_chainB = copy.deepcopy(AF)

    AF_chainA.df["ATOM"] = AF_chainA.df["ATOM"].query("`chain_id` in @chain_A")
    AF_chainB.df["ATOM"] = AF_chainB.df["ATOM"].query("`chain_id` in @chain_B")

    AF_chainA.df["ATOM"].loc[:, "chain_id"] = "A"
    AF_chainB.df["ATOM"].loc[:, "chain_id"] = "B"

    return AF_chainA, AF_chainB

def get_interface_2(file_path, num_chains_protein_B, interface_cutoff):
    results = []
    AF = PandasPdb().read_pdb(str(file_path))
    chains_AF = []
    for chain_ID in AF.df['ATOM']['chain_id'].unique():
        chains_AF.append(chain_ID)
    
    total_number_chains = len(chains_AF)
    num_chains_protein_A = total_number_chains - int(num_chains_protein_B)
    chains_AF_A = chains_AF[:num_chains_protein_A]
    chains_AF_B = chains_AF[num_chains_protein_A:]
    
    for ind in AF.df['ATOM']['chain_id'].index:
        current_chain_ID = AF.df['ATOM']['chain_id'][ind]
        if current_chain_ID in chains_AF_A:
            AF.df['ATOM'].at[ind, 'chain_id'] = 'A'
        else:
            AF.df['ATOM'].at[ind, 'chain_id'] = 'B'
        
    AF_chainA = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'A']
    AF_chainB = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'B']
    
    AF_chainA_CB = AF_chainA.loc[(AF_chainA['atom_name'] == 'CB') | ((AF_chainA['residue_name'] == 'GLY') & (AF_chainA['atom_name'] == 'CA'))]
    AF_chainB_CB = AF_chainB.loc[(AF_chainB['atom_name'] == 'CB') | ((AF_chainB['residue_name'] == 'GLY') & (AF_chainB['atom_name'] == 'CA'))]
    
    for ind in AF_chainA_CB.index:
        X_coord_A = AF_chainA_CB['x_coord'][ind]
        Y_coord_A = AF_chainA_CB['y_coord'][ind]
        Z_coord_A = AF_chainA_CB['z_coord'][ind]
        plDDt_A = AF_chainA_CB['b_factor'][ind]
        resname_A = AF_chainA_CB['residue_name'][ind]
        res_number_A = AF_chainA_CB['residue_number'][ind]
        atom_number_A = AF_chainA_CB['atom_number'][ind]
        for B_ind in AF_chainB_CB.index:
            X_coord_B = AF_chainB_CB['x_coord'][B_ind]
            Y_coord_B = AF_chainB_CB['y_coord'][B_ind]
            Z_coord_B = AF_chainB_CB['z_coord'][B_ind]
            plDDt_B = AF_chainB_CB['b_factor'][B_ind]
            resname_B = AF_chainB_CB['residue_name'][B_ind]
            res_number_B = AF_chainB_CB['residue_number'][B_ind]
            atom_number_B = AF_chainB_CB['atom_number'][ind]
            distance = np.sqrt((X_coord_A - X_coord_B) **2 + (Y_coord_A - Y_coord_B) **2 + (Z_coord_A - Z_coord_B) **2)
            if distance <= int(interface_cutoff):
                results.append({'resname_A': resname_A, 'res_numberA': res_number_A, 'atom_numberA': atom_number_A, 'resname_B': resname_B, 'res_numberB': res_number_B,'atom_numberB': atom_number_B, 'plDDt_A': plDDt_A, 'plDDt_B': plDDt_B, 'distance': distance})
    results_df = pd.DataFrame(results)




def split_chains_2(file_path, num_chains_protein_B):
    AF = PandasPdb().read_pdb(str(file_path))
    AF_chainA = PandasPdb().read_pdb(file_path)
    AF_chainB = PandasPdb().read_pdb(file_path)

    chains = AF.df["ATOM"]["chain_id"].unique()
    chain_A = chains[:-int(num_chains_protein_B)]
    chain_B = chains[-int(num_chains_protein_B):]

    chains_AF = []
    for chain_ID in AF.df['ATOM']['chain_id'].unique():
        chains_AF.append(chain_ID)
        total_number_chains = len(chains)
        num_chains_protein_A = total_number_chains - int(num_chains_protein_B)
        chains_AF_A = chains_AF[:num_chains_protein_A]
        chains_AF_B = chains_AF[num_chains_protein_A:]
        ####Replace chain IDS to be 'A' and 'B'
        for ind in AF.df['ATOM']['chain_id'].index:
            current_chain_ID = AF.df['ATOM']['chain_id'][ind]
            if current_chain_ID in chains_AF_A:
                AF.df['ATOM'].at[ind, 'chain_id'] = 'A'
            if current_chain_ID in chains_AF_B:
                AF.df['ATOM'].at[ind, 'chain_id'] = 'B'
    AF_chainA.df['ATOM'] = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'A']
    AF_chainB.df['ATOM'] = AF.df['ATOM'][AF.df['ATOM']['chain_id'] == 'B']
    return AF_chainA, AF_chainB, AF

def atom_num_to_res_num(file_path, num_chains_protein_B):
     AF = PandasPdb().read_pdb(str(file_path))




if __name__ == '__main__':
    args = args()
    chainsA, chainsB, AF = split_chains_2(args.file_path, args.num_chains_protein_B)
 
    # interface = get_interface(args.file_path, args.num_chains_protein_B, args.interface_cutoff)
    # working_directory = "test_cases"
    # csv_name = "test_df"
    
    # combined_csv_name = f"{working_directory}/{csv_name}"
    # interface.to_csv(f"{combined_csv_name}.csv", mode='w', header=True, index=False)