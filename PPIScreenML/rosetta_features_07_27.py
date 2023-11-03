#!/usr/bin/env python
# coding: utf-8

import os
import argparse
from pathlib import Path
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
from pyrosetta import pose_from_sequence
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
import pyrosetta.toolbox
import pyrosetta; pyrosetta.init()
from pyrosetta import *
from pyrosetta.teaching import *
init('-add_orbitals')
import utils

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
    dict_list = []
    if not os.path.exists(f"{working_directory}/chains/"):
        os.mkdir(f"{working_directory}/chains/")
    for i,file in enumerate(pdb_files):
        if file != ".DS_Store":
            file_name = file.stem
            print(file_name)

            file_path = f"{working_directory}/{file_name}.pdb"
            AF = PandasPdb().read_pdb(file_path)
            
            AF_chainA_pandas, AF_chainB_pandas, AF = utils.split_chains_2(file_path, int(num_chains_protein_B))
            
            chainA_file = f"{working_directory}/chains/{file_name}chainA.pdb"
            chainB_file = f"{working_directory}/chains/{file_name}chainB.pdb"
            complex_file = f"{working_directory}/chains/{file_name}complex.pdb"

            fp = open(chainA_file, 'w')
            fp.close()
            fp2 = open(chainB_file, 'w')
            fp2.close()
            
            AF.to_pdb(complex_file, records = None, gz=False, append_newline = False)
            AF_chainA_pandas.to_pdb(chainA_file, records = None, gz=False, append_newline = False)
            AF_chainB_pandas.to_pdb(chainB_file, records = None, gz=False, append_newline = False)

            interface_df = utils.get_interface(file_path, num_chains_protein_B, interface_cutoff)
            ####format pandas dataframes to fit into pyrossetta
            pose = pose_from_pdb(f"{working_directory}/chains/{file_name}complex.pdb")
            chain_A_pose = pose_from_pdb(chainA_file)
            chain_B_pose = pose_from_pdb(chainB_file)

            sfxn = get_score_function(True)
            sfxn(pose)


            ####Calculate binding energy
            BE_complex = sfxn(pose)
            BE_chainA = sfxn(chain_A_pose)
            BE_chainB = sfxn(chain_B_pose)
            Binding_energy = BE_complex - BE_chainA - BE_chainB

            fa_atr_complex = pose.energies().total_energies()[ScoreType.fa_atr]
            fa_atr_chainA = chain_A_pose.energies().total_energies()[ScoreType.fa_atr]
            fa_atr_chainB = chain_B_pose.energies().total_energies()[ScoreType.fa_atr]
            fa_atr_binding = fa_atr_complex - fa_atr_chainA - fa_atr_chainB

            fa_rep_complex = pose.energies().total_energies()[ScoreType.fa_rep]
            fa_rep_chainA = chain_A_pose.energies().total_energies()[ScoreType.fa_rep]
            fa_rep_chainB = chain_B_pose.energies().total_energies()[ScoreType.fa_rep]
            fa_rep_binding = fa_rep_complex - fa_rep_chainA - fa_rep_chainB


            fa_sol_complex = pose.energies().total_energies()[ScoreType.fa_sol]
            fa_sol_chainA = chain_A_pose.energies().total_energies()[ScoreType.fa_sol]
            fa_sol_chainB = chain_B_pose.energies().total_energies()[ScoreType.fa_sol]
            fa_sol_binding = fa_sol_complex - fa_sol_chainA - fa_sol_chainB


            fa_elec_complex = pose.energies().total_energies()[ScoreType.fa_elec]
            fa_elec_chainA = chain_A_pose.energies().total_energies()[ScoreType.fa_elec]
            fa_elec_chainB = chain_B_pose.energies().total_energies()[ScoreType.fa_elec]
            fa_elec_binding = fa_elec_complex - fa_elec_chainA - fa_elec_chainB

            lk_ball_wtd_complex = pose.energies().total_energies()[ScoreType.lk_ball_wtd]
            lk_ball_wtd_chainA = chain_A_pose.energies().total_energies()[ScoreType.lk_ball_wtd]
            lk_ball_wtd_chainB = chain_B_pose.energies().total_energies()[ScoreType.lk_ball_wtd]
            lk_ball_wtd_binding = lk_ball_wtd_complex - lk_ball_wtd_chainA - lk_ball_wtd_chainB

            pro_close_complex = pose.energies().total_energies()[ScoreType.pro_close]
            p_aa_pp_complex = pose.energies().total_energies()[ScoreType.p_aa_pp]
            fa_dun_complex = pose.energies().total_energies()[ScoreType.fa_dun]
            rama_prepro_complex = pose.energies().total_energies()[ScoreType.rama_prepro]


            ###Get number of H bonds formed between interface residue in A vs B
            hbond_set = pose.get_hbonds()
            Num_IF_H_bonds = 0
            for ind in interface_df.index:
                resnum_A = interface_df['res_numberA'][ind]
                resnum_B = interface_df['res_numberB'][ind]
                resA = pose.pdb_info().pdb2pose("A", resnum_A)
                resB = pose.pdb_info().pdb2pose("B", resnum_B)
                for hbond in hbond_set.hbonds():
                    if hbond.don_res() == resA and hbond.acc_res() == resB:
                        Num_IF_H_bonds +=1
                    if hbond.don_res() == resB and hbond.acc_res() == resA:
                        Num_IF_H_bonds +=1
                    else:
                        continue

        
            ####Get energy between residues in interface
            list_atr = []
            list_rep = []
            list_sol = []
            list_elec = []
            list_fa_dun = []
            for ind in interface_df.index:
                resnum_A = interface_df['res_numberA'][ind]
                resnum_B = interface_df['res_numberB'][ind]
                resA = pose.pdb_info().pdb2pose("A", resnum_A)
                resB = pose.pdb_info().pdb2pose("B", resnum_B)
                emap = EMapVector()
                sfxn.eval_ci_2b(pose.residue(resA), pose.residue(resB), pose, emap)
                values_atr = emap[fa_atr]
                values_rep = emap[fa_rep]
                values_sol = emap[fa_sol]
                values_elec = emap[fa_elec]
                

                fa_dun_score_resA = pose.energies().residue_total_energies(resA)[fa_dun]
                fa_dun_score_resB = pose.energies().residue_total_energies(resB)[fa_dun]
            
                list_atr.append(values_atr)
                list_rep.append(values_rep)
                list_sol.append(values_sol)
                list_elec.append(values_elec)
                list_fa_dun.append(fa_dun_score_resA)
                list_fa_dun.append(fa_dun_score_resB)


            avg_fa_dun = np.average(list_fa_dun)
            avg_atr = np.average(list_atr)
            avg_rep = np.average(list_rep)
            avg_sol = np.average(list_sol)
            avg_elec = np.average(list_elec)

        
            list_hbond_bb_sc = []
            for val in interface_df.index:
                resnum_A = interface_df['res_numberA'][val]
                resnum_B = interface_df['res_numberB'][val]
                residue_energy_map_A = pose.energies().residue_total_energies(resnum_A)
                residue_energy_map_B = pose.energies().residue_total_energies(resnum_B)
           
                hbond_bb_sc_energy_A = residue_energy_map_A[ScoreType.hbond_bb_sc]
                hbond_bb_sc_energy_B = residue_energy_map_B[ScoreType.hbond_bb_sc]
                
                
             
                list_hbond_bb_sc.append(hbond_bb_sc_energy_A)
                list_hbond_bb_sc.append(hbond_bb_sc_energy_B)



            avg_hbond_bb_sc = np.average(list_hbond_bb_sc)

            unique_interface_res = set()
            for ind in interface_df.index:
                resnum_A = interface_df['res_numberA'][ind]
                resnum_B = interface_df['res_numberB'][ind]
                resA = pose.pdb_info().pdb2pose("A", resnum_A)
                resB = pose.pdb_info().pdb2pose("B", resnum_B)
                unique_interface_res.add(resA)
                unique_interface_res.add(resB)
            

            ###Get secondary structure characteristics 
            DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
            DSSP.apply(pose)
            ss = pose.secstruct()
            interface_SS = []
            for val in unique_interface_res:
                value = val -1 
                struct = ss[value]
                interface_SS.append(struct)
            IF_SS_helix = 0
            IF_SS_sheet = 0
            IF_SS_loop = 0
            for x in interface_SS:
                if x == 'H':
                    IF_SS_helix += 1
                if x == 'E':
                    IF_SS_sheet += 1
                if x == 'L':
                    IF_SS_loop += 1

            if IF_SS_helix != 0 or IF_SS_sheet != 0 or IF_SS_loop !=0:
                prop_IF_SS_helix = IF_SS_helix / len(interface_SS)
                prop_IF_SS_sheet = IF_SS_sheet / len(interface_SS)
                prop_IF_SS_loop = IF_SS_loop /len(interface_SS)
            else: 
                prop_IF_SS_helix = 0
                prop_IF_SS_sheet = 0 
                prop_IF_SS_loop = 0               


            # In[12]:


            ####Get the SASA
            sasa_calc = sasa.SasaCalc()
            SASA_complex = sasa_calc.calculate(pose)
            SASA_chainA = sasa_calc.calculate(chain_A_pose)
            SASA_chainB = sasa_calc.calculate(chain_B_pose)
            buried_SASA =  (SASA_chainA + SASA_chainB) - SASA_complex 



            # In[24]:

            interface_cutoff_str = str(interface_cutoff)

            BE_complex_name = f"BE_complex"
            BE_bining_energy_name = f"BE_bining_energy"
            fa_atr_complex_name = f"fa_atr_complex"
            fa_atr_binding_name = f"fa_atr_binding"
            fa_rep_complex_name = f"fa_rep_complex"
            fa_rep_binding_name = f"fa_rep_binding"
            fa_sol_complex_name = f"fa_sol_complex"
            fa_sol_binding_name = f"fa_sol_binding"
            fa_elec_complex_name = f"fa_elec_complex"
            fa_elec_binding_name = f"fa_elec_binding"
            lk_ball_wtd_complex_name = f"lk_ball_wtd_complex"
            lk_ball_wtd_binding_name = f"lk_ball_wtd_binding"
            pro_close_complex_name = f"pro_close_complex"
            p_aa_pp_complex_name = f"p_aa_pp_complex"
            fa_dun_complex_name = f"fa_dun_complex"
            rama_prepro_complex_name = f"rama_prepro_complex"
            Num_IF_H_bonds_name = f"Num_IF_H_bonds"
            avg_IF_atr_name = f"avg_IF_atr"
            avg_IF_rep_name = f"avg_IF_rep"
            avg_IF_sol_name = f"avg_IF_sol"
            avg_IF_elec_name = f"avg_IF_elec"
            avg_fa_dun_name = f"avg_IF_fa_dun"
            avg_hbond_bb_sc_name = f"avg_hbond_bb_sc"
            prop_IF_SS_helix_name = f"prop_IF_SS_helix"
            prop_IF_SS_sheet_name = f"prop_IF_SS_sheet"
            prop_IF_SS_loop_name = f"prop_IF_SS_loop"
            SASA_complex_name = f"SASA_complex"
            buried_SASA_name = f"buried_SASA"


            data = {'File_name': file_name,
                    BE_complex_name: BE_complex,
                    BE_bining_energy_name: Binding_energy,
                    fa_atr_complex_name: fa_atr_complex,
                    fa_atr_binding_name: fa_atr_binding,
                    fa_rep_complex_name: fa_rep_complex,
                    fa_rep_binding_name: fa_rep_binding,
                    fa_sol_complex_name: fa_sol_complex,
                    fa_sol_binding_name: fa_sol_binding,
                    fa_elec_complex_name: fa_elec_complex,
                    fa_elec_binding_name: fa_elec_binding,
                    lk_ball_wtd_complex_name: lk_ball_wtd_complex,
                    lk_ball_wtd_binding_name: lk_ball_wtd_binding,
                    pro_close_complex_name: pro_close_complex,
                    p_aa_pp_complex_name: p_aa_pp_complex,
                    fa_dun_complex_name: fa_dun_complex,
                    rama_prepro_complex_name: rama_prepro_complex,
                    Num_IF_H_bonds_name: Num_IF_H_bonds,
                    avg_IF_atr_name: avg_atr,
                    avg_IF_rep_name: avg_rep,
                    avg_IF_sol_name: avg_sol,
                    avg_IF_elec_name: avg_elec,
                    avg_fa_dun_name: avg_fa_dun,
                    avg_hbond_bb_sc_name: avg_hbond_bb_sc,
                    prop_IF_SS_helix_name: prop_IF_SS_helix,
                    prop_IF_SS_sheet_name: prop_IF_SS_sheet,
                    prop_IF_SS_loop_name: prop_IF_SS_loop,
                    SASA_complex_name: SASA_complex,
                    buried_SASA_name: buried_SASA}

            dict_list.append(data)

            os.remove(chainA_file)
            os.remove(chainB_file)
            os.remove(complex_file)
            


    rosetta_df = pd.DataFrame.from_dict(dict_list)


    return rosetta_df

if __name__ == '__main__':
    args = args()
    features = main(args.working_directory, args.num_chains_protein_B, args.interface_cutoff)
    
    # Assign variables based on command line arguments
    working_directory = args.working_directory
    csv_name = args.csv_name
    combined_csv_name = f"{working_directory}/{csv_name}"
    features.to_csv(f"{combined_csv_name}.csv", mode='w', header=True, index=False)
