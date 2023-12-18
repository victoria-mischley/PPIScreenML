import get_features
import pandas as pd
import argparse
import xgboost as xgb
import numpy as np

def args():
	parser = argparse.ArgumentParser(description='Combine dataframes')

	# Add command line arguments
	parser.add_argument('--working_directory', type=str, help='path to location of AF folder output')
	parser.add_argument('--protein1_chains_input', nargs='+', help='List of chains for Protein 1')
	parser.add_argument('--protein2_chains_input', nargs='+', help='List of chains for Protein 2')
	parser.add_argument('--csv_name', type=str, help='the name that you would like your csv file to be saved as')

	# Parse the command line arguments
	args = parser.parse_args()

	return args

############################################################################################################################


def get_scores(X):
	# final_model = joblib.load(f"PPI_xgboost_afv3.pkl")
    data_dmatrix = xgb.DMatrix(X)
    final_model = xgb.Booster()
    final_model.load_model('PPIScreenML_model.model')
    raw_scores = final_model.predict(data_dmatrix, output_margin=True)
    model_scores = 1 / (1 + np.exp(-raw_scores))
    threshold = 0.5
    predicted_label = (model_scores >= threshold).astype(int)
    
    return model_scores, predicted_label


def main(working_directory, protein1_chains_input, protein2_chains_input, csv_name, features_file_path):
    features_df = pd.read_csv(features_file_path)
    ###Get columns for reduced feature model#### Df is offset by one
    columns_red_feats =  [4, 14, 15, 16, 17, 21, 28]
    df_columns_red_feats = [-1, 4, 14, 15, 16, 17, 21, 28]

    columns_red_feats_std = [val + 1 for val in columns_red_feats]
    df_columns_red_feats_std = [val + 1 for val in df_columns_red_feats]

    ###Drop the filenames that don't have an interface- these will automatically be predicted as not interacting
    nonans_red_feats = features_df.iloc[:, columns_red_feats_std].dropna()
    X = nonans_red_feats.values
    X_df = features_df.iloc[:, df_columns_red_feats_std].dropna()
    print(X_df)
    ####Get the names of the files that had no interface and were dropped and assign as 0
    file_names_all = features_df['File_name']
    file_names_X_df = X_df['File_name']
    file_names_dropped = file_names_all[~file_names_all.isin(file_names_X_df)]

    results_dropped_files = pd.DataFrame()
    if len(file_names_dropped) != 0:
        for val in file_names_dropped:
            data = {
                'File_name': val,
                'score': 0,
                'predicted_label': 0
            }
        results_dropped_files = results_dropped_files.append(data, ignore_index=True)


    ####Get predicted score
    model_scores, predicted_label = get_scores(X)
    
    ####get predicted results for model
    results_not_dropped = pd.DataFrame({
        'File_name': X_df['File_name'],
        'score': model_scores,
        'predicted_label': predicted_label
    })

    ###Add in the ones that did not have any interface
    results = pd.concat([results_dropped_files, results_not_dropped], axis=0, ignore_index=True)
    results['predicted_label'] = results['predicted_label'].round().astype(int)
    results.to_csv(f"{combined_csv_name}_classification_results.csv", mode='w', header=True, index=False)


    def label_to_prediction(label):
        return 'interacting' if label == 1 else 'not interacting'
    results['prediction'] = results['predicted_label'].apply(label_to_prediction)
    ###Group all 5 models together based on name
    results['common_name'] = results['File_name'].apply(lambda x: x.split("_relaxed_")[0])
    # Rank by score within each common_name group
    results['prob_rank'] = results.groupby('common_name')['score'].rank(method='first', ascending=False).astype(int)
    # Filter rows where prob_rank is 1
    filtered_df = results[results['prob_rank'] == 1]
    final_df = filtered_df.drop(columns=['prob_rank', 'common_name'])
    final_df.to_csv(f"{combined_csv_name}_classification_results_topmodel.csv", mode='w', header=True, index=False)
    
    print(results)

if __name__ == '__main__':
    args = args()
    working_directory = args.working_directory
    csv_name = args.csv_name
    combined_csv_name = f"{working_directory}/{csv_name}"
        ####get featurews df        
    features_file_path = f"{working_directory}/{csv_name}_features.csv"
    main(working_directory, args.protein1_chains_input, args.protein2_chains_input, csv_name, features_file_path)

    
      
