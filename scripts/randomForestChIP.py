import argparse
import numpy as np
import pandas as pd

from typing import cast
from pathlib import Path

from rates2chip.utilities import import_pandas, export_pandas
from rates2chip.random_forest import randomForest, clean_data

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='input rate path')
	parser.add_argument('--target',help='target column')
	parser.add_argument('--features',nargs='+',help='specific features to include')
	parser.add_argument('--exclude',nargs='+',help='if features is empty, use all columns excluding these')
	parser.add_argument('--filters',default=[],nargs='+',help='List of filters')
	parser.add_argument('--plot_path',required=False,help='path to create plot of observed vs expected')
	parser.add_argument('--seed',requried=False,type=int,help='seed for random forest regression')
	parser.add_argument('--out',help='output path')
	args = parser.parse_args()

	# Import data
	df = import_pandas(Path(args.input))

	# Subset data
	for x in args.filters: 
		if '==' in x:
			df = df[df[x.split('==')[0]] == float(x.split('==')[1])]
		elif '<=' in x:
			df = df[df[x.split('<=')[0]] == float(x.split('<=')[1])]
		elif '>=' in x:
			df = df[df[x.split('>=')[0]] == float(x.split('>=')[1])]
		elif '!=' in x:
			df = df[df[x.split('!=')[0]] != float(x.split('!=')[1])]
		elif '>' in x:
			df = df[df[x.split('>')[0]] > float(x.split('>')[1])]
		elif '<' in x:
			df = df[df[x.split('<')[0]] == float(x.split('<')[1])]
		else: raise Exception('Invalid filter')

	# Identify target and features
	target = args.target
	if len(args.features) > 0: features = [f for f in args.features if f in df.columns]
	else: features = [f for f in df.columns if f not in args.exclude]
	include_columns = [target] + features
	
	# Subset data
	df = df[include_columns]
	df = clean_data(df,target,[target],features)

	# Collect properties for random forest run
	if args.plot_path: plot_path = Path(args.plot_path)
	else: plot_path = None
	if args.seed: seed = args.seed
	else: seed = cast(int,np.random.SeedSequence().entropy)

	# Run regression
	r2, mse, out_df = randomForest(df,target,True,seed,False,plot_path)
	out_df.sort_values(by='Importance', ascending=False, inplace=True)
	print(f'R2: {r2:.3f}')
	print(f'MSE: {mse:.3f}')

	# Save outputs
	export_pandas(out_df,Path(args.out))