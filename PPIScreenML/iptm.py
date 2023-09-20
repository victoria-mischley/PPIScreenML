import pandas as pd
import json
from pathlib import Path
import argparse
import utils

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('working_directory', type=str, help='path to location of AF folder output')
    parser.add_argument('num_chains_protein_B', type=str, help='number of chains of the last protein')
    parser.add_argument('interface_cutoff', type=str, help='distance that defines the interface')
    parser.add_argument('csv_name', type=str, help='the name that you would like your csv file to be saved as')
    return parser.parse_args()

def main(working_directory):
    directory = Path(working_directory)
    pdb_files  = list(directory.glob('*.pdb'))
    all_json_files = list(directory.glob('*.json'))
    json_files = [file for file in all_json_files if 'error' not in file.name]

    dict_list = []

    for file in pdb_files:
        if file.name != ".DS_Store":
            file_name = file.stem
            pdb_prefix = str(file_name).split('_unrelaxed_')[0] if '_unrelaxed_' in str(file_name) else str(file_name).split('_relaxed_')[0]
            pdb_suffix = str(file_name).split('_unrelaxed_')[1] if '_unrelaxed_' in str(file_name) else str(file_name).split('_relaxed_')[1]
            json_file_string = next((str(file) for file in json_files if pdb_prefix in str(file) and pdb_suffix in str(file)), None)
            file_path = file
            json_file_path = str(json_file_string)
            print(file_name)
            f = open(json_file_path, "r")
            json_data = json.loads(f.read())
            iptm = json_data['iptm']
            iptm_name = "iptm"
            data = {'File_name': file_name, iptm_name: iptm}
            dict_list.append(data)

    return pd.DataFrame.from_dict(dict_list)

if __name__ == '__main__':
    args = parse_args()
    features = main(args.working_directory)
    combined_csv_name = f"{args.working_directory}/{args.csv_name}"
    features.to_csv(f"{combined_csv_name}.csv", mode='w', header=True, index=False)
