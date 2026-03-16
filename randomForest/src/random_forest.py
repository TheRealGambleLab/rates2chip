import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

from pathlib import Path

from scipy.stats import spearmanr, gaussian_kde

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import quantile_transform

def clean_data(df:pd.DataFrame,target:str,log_columns:list[str],nan_columns:list[str]) -> pd.DataFrame:
	# Switch data to log scale
	for c in log_columns:
		df.loc[df[c] <= 0,c] = np.nan
		df[c] = np.log10(df[c])

	# Change NaN to mean value
	for c in nan_columns:
		nan_idx = df[c].isna()
		df.loc[nan_idx,c] = np.mean(df.loc[~nan_idx,c])

	# Drop nan and reset index
	df.dropna(subset=[target],inplace=True)
	df.reset_index(drop=True,inplace=True)

	return df

def randomForest(df:pd.DataFrame,target:str,include_spearman:bool=True,seed:int=42,quantile_target:bool=False,prediction_plot:Path|None=None) -> tuple[float,float,pd.DataFrame]:
	if quantile_target:
		n_q = int(len(df.index))
		df[target] = quantile_transform(df[target].to_frame(), n_quantiles=n_q,output_distribution='uniform',copy=True,subsample=n_q).squeeze()

	# Clean data
	x = df.drop(columns=[target],axis=1)
	y = df[target]

	# Split into train/test
	x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.1, random_state=seed)

	# Fit model
	model = RandomForestRegressor(n_estimators=200, random_state=seed, min_samples_split=10, min_samples_leaf=5, n_jobs=-1)
	model.fit(x_train, y_train)

	# Predict
	y_pred = model.predict(x_test)

	if prediction_plot:
		plt.figure(figsize=(8, 6))
		density = gaussian_kde([y_test,y_pred])([y_test,y_pred])
		sb.scatterplot(x=10**y_test, y=10**y_pred, hue=density, alpha=0.6,legend=False)
		plt.xscale('log')
		plt.yscale('log')
		plt.plot([10**y_test.min(), 10**y_test.max()], [10**y_test.min(), 10**y_test.max()], 'k--')
		plt.xlabel("Observed")
		plt.ylabel("Predicted")
		plt.title("Observed vs Predicted")
		plt.savefig(prediction_plot)
		plt.close()

	# Evaluate
	r2 = r2_score(y_test, y_pred)
	print(r2)
	mse = mean_squared_error(y_test,y_pred)
	weights = permutation_importance(model,x_train,y_train,n_repeats=100, random_state=seed,n_jobs=-1)

	# Store in dataframe
	data = pd.DataFrame({'Feature':x.columns,'Importance':model.feature_importances_,'Mean_Importance':weights.importances_mean,'Standard_Deviation':weights.importances_std})

	# Calculate correlations
	if include_spearman:
		data['Correlation'] = 0.0
		data['Pvalue'] = 0.0
		for c in x.columns:
			spearman = spearmanr(df[target].to_numpy(),df[c].to_numpy())
			idx = data['Feature'] == c
			data.loc[idx,'Correlation'] = spearman.statistic
			data.loc[idx,'Pvalue'] = spearman.pvalue

	return r2, mse, data