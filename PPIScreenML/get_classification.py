from get_features_08_03 import get_features
import pandas as pd
import argparse
import xgboost as xgb
import numpy as np

def args():
	parser = argparse.ArgumentParser(description='Combine dataframes')

	# Add command line arguments
	parser.add_argument('working_directory', type=str, help='path to location of ColabFold folder output')
	parser.add_argument('num_chains_protein_B', type=str, help='number of chains of the last protein')
	parser.add_argument('csv_name', type=str, help='the name that you would like your csv file to be saved as')

	# Parse the command line arguments
	args = parser.parse_args()

	return args

############################################################################################################################

def get_red_features(working_directory, num_chains_protein_B, interface_cutoff, combined_csv_name):
    df = get_features(working_directory, num_chains_protein_B, interface_cutoff)
    df.to_csv(f"{combined_csv_name}_features.csv", mode='w', header=True, index=False)

    columns_red_feats =  [2, 4, 5, 8, 9, 15]
    df_columns_red_feats = [-1, 2, 4, 5, 8, 9, 15]

    columns_red_feats_std = [val + 1 for val in columns_red_feats]
    df_columns_red_feats_std = [val + 1 for val in df_columns_red_feats]

    nonans_red_feats = df.iloc[:, columns_red_feats_std].dropna()
    X = nonans_red_feats.values
    X_df = df.iloc[:, df_columns_red_feats_std].dropna()

    return X, X_df

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

def main(working_directory, num_chains_protein_B, interface_cutoff, combined_csv_name):
    X, X_df = get_red_features(working_directory, num_chains_protein_B, interface_cutoff, combined_csv_name)
    model_scores, predicted_label = get_scores(X)
    
    results = pd.DataFrame({
        'file_name': X_df['File_name'],
        'score': model_scores,
        'predicted_label': predicted_label
    })

    def label_to_prediction(label):
        return 'interacting' if label == 1 else 'not interacting'
    results['prediction'] = results['predicted_label'].apply(label_to_prediction)
    ###Group all 5 models together based on name
    results['common_name'] = results['file_name'].apply(lambda x: x.split("_relaxed_")[0])
    # Rank by score within each common_name group
    results['prob_rank'] = results.groupby('common_name')['score'].rank(method='first', ascending=False).astype(int)
    # Filter rows where prob_rank is 1
    filtered_df = results[results['prob_rank'] == 1]
    final_df = filtered_df.drop(columns=['prob_rank', 'common_name'])
    final_df.to_csv(f"{combined_csv_name}_classification_results.csv", mode='w', header=True, index=False)


if __name__ == '__main__':
    args = args()
    interface_cutoff = int(12)
    working_directory = args.working_directory
    csv_name = args.csv_name
    num_chains_protein_B = args.num_chains_protein_B
    combined_csv_name = f"{working_directory}/{csv_name}"
    main(working_directory, num_chains_protein_B, interface_cutoff, combined_csv_name)
    
      
