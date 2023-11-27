#!/usr/bin/env python
# coding: utf-8
import Fix_formatting
import get_interface_dataframe
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

def args():
    parser = argparse.ArgumentParser(description='Combine dataframes')

    # Add command line arguments
    #parser.add_argument('--working_directory', type=str, help='path to location of AF folder output')
    parser.add_argument('--file_path', type=str, help='path to location of AF folder output')
    parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
    parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
    
    # Parse the command line arguments
    args = parser.parse_args()
    return args

def main(working_directory, file_name, results_df, AF_fixed, AF_chainA, AF_chainB):
    
    if not os.path.exists(f"{working_directory}/temp/"):
        os.mkdir(f"{working_directory}/temp/")
    chainA_file = f"{working_directory}/temp/{file_name}chainA.pdb"
    chainB_file = f"{working_directory}/temp/{file_name}chainB.pdb"
    complex_file = f"{working_directory}/temp/{file_name}complex.pdb"

    fp = open(chainA_file, 'w')
    fp.close()
    fp2 = open(chainB_file, 'w')
    fp2.close()
            
    AF_fixed.to_pdb(complex_file, records = None, gz=False, append_newline = False)
    AF_chainA.to_pdb(chainA_file, records = None, gz=False, append_newline = False)
    AF_chainB.to_pdb(chainB_file, records = None, gz=False, append_newline = False)

    dict_list = []
            ####format pandas dataframes to fit into pyrossetta
    pose = pose_from_pdb(complex_file)
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

     ####Get the SASA
    sasa_calc = sasa.SasaCalc()
    SASA_complex = sasa_calc.calculate(pose)
    SASA_chainA = sasa_calc.calculate(chain_A_pose)
    SASA_chainB = sasa_calc.calculate(chain_B_pose)
    buried_SASA =  (SASA_chainA + SASA_chainB) - SASA_complex 


            ###Get number of H bonds formed between interface residue in A vs B
    if len(results_df) != 0:
        hbond_set = pose.get_hbonds()
        Num_IF_H_bonds = 0
        for ind in results_df.index:
            resnum_A = results_df['res_numberA'][ind]
            resnum_B = results_df['res_numberB'][ind]
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
        for ind in results_df.index:
            resnum_A = results_df['res_numberA'][ind]
            resnum_B = results_df['res_numberB'][ind]
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
        for val in results_df.index:
            resnum_A = results_df['res_numberA'][val]
            resnum_B = results_df['res_numberB'][val]
            residue_energy_map_A = pose.energies().residue_total_energies(resnum_A)
            residue_energy_map_B = pose.energies().residue_total_energies(resnum_B)
            
            hbond_bb_sc_energy_A = residue_energy_map_A[ScoreType.hbond_bb_sc]
            hbond_bb_sc_energy_B = residue_energy_map_B[ScoreType.hbond_bb_sc]
                    
                    
                
            list_hbond_bb_sc.append(hbond_bb_sc_energy_A)
            list_hbond_bb_sc.append(hbond_bb_sc_energy_B)



        avg_hbond_bb_sc = np.average(list_hbond_bb_sc)

        unique_interface_res = set()
        for ind in results_df.index:
            resnum_A = results_df['res_numberA'][ind]
            resnum_B = results_df['res_numberB'][ind]
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
        
        data = { 'Binding_energy': Binding_energy,
                'fa_atr_binding': fa_atr_binding,
                'fa_rep_binding': fa_rep_binding,
                'fa_sol_binding': fa_sol_binding,
                'fa_elec_binding': fa_elec_binding,
                'lk_ball_wtd_binding': lk_ball_wtd_binding,
                'Num_IF_H_bonds': Num_IF_H_bonds,
                'avg_IF_atr': avg_atr,
                'avg_IF_rep': avg_rep,
                'avg_IF_sol': avg_sol,
                'avg_IF_elec': avg_elec,
                'avg_IF_fa_dun': avg_fa_dun,
                'avg_hbond_bb_sc': avg_hbond_bb_sc,
                'prop_IF_SS_helix': prop_IF_SS_helix,
                'prop_IF_SS_sheet': prop_IF_SS_sheet,
                'prop_IF_SS_loop': prop_IF_SS_loop,
                'buried_SASA_name': buried_SASA}               

    else:
        data = { 'Binding_energy': Binding_energy,
                'fa_atr_binding': fa_atr_binding,
                'fa_rep_binding': fa_rep_binding,
                'fa_sol_binding': fa_sol_binding,
                'fa_elec_binding': fa_elec_binding,
                'lk_ball_wtd_binding': lk_ball_wtd_binding,
                'Num_IF_H_bonds': None,
                'avg_IF_atr': None,
                'avg_IF_rep': None,
                'avg_IF_sol': None,
                'avg_IF_elec': None,
                'avg_IF_fa_dun': None,
                'avg_hbond_bb_sc': None,
                'prop_IF_SS_helix': None,
                'prop_IF_SS_sheet': None,
                'prop_IF_SS_loop': None,
                'buried_SASA_name': buried_SASA}  

            # In[12]:








    dict_list.append(data)
    
    rosetta_df = pd.DataFrame.from_dict(dict_list)

    os.remove(chainA_file)
    os.remove(chainB_file)
    os.remove(complex_file)



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
    file_path = Path(args.file_path)
    working_directory = file_path.parent
    file_name = file_path.name

            
    
    main(working_directory, file_name, results_df, AF_fixed, AF_chainA, AF_chainB)


            
