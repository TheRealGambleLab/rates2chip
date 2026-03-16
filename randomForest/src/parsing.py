import pyBigWig
import numpy as np
import pandas as pd

from pathlib import Path

def melt_rates(rate_df:pd.DataFrame) -> pd.DataFrame:
	columns = ['gene','start','stop','strand','chromosome','Gb','Gs','spawn','drb','spawnGs','contamination','unreleased','alpha','rate','first_rate','first_cleave','total_cleave','median_splice','mean_splice','geomean_splice']
	data = {x:[] for x in columns}
	for gene,df in rate_df.groupby('gene',observed=True):
		local_data = {}

		# Get basic properties
		idx = df['type'] == 'tauC'
		for c in ['gene','start','stop','chromosome','strand']:local_data[c] = df.loc[idx,c].to_numpy()[0]

		# Get avg elongation rate
		valid_chk = all([df.loc[idx,c].to_numpy()[0] for c in ['replicated','converged','valid_dependencies']]) and not any([df.loc[idx,c].to_numpy()[0] for c in ['upper','lower']])
		if valid_chk: local_data['rate'] = (local_data['stop'] - local_data['start']) / df.loc[idx,'value'].to_numpy()[0]
		else: local_data['rate'] = np.nan

		# Get basic rates
		for x in ['Gb','Gs','spawn','drb']:
			idx = df['type'] == x
			val = df.loc[idx,'value'].to_numpy()[0]
			valid_chk = all([df.loc[idx,c].to_numpy()[0] for c in ['replicated','converged','valid_dependencies']]) and not any([df.loc[idx,c].to_numpy()[0] for c in ['upper','lower']])
			if valid_chk: local_data[x] = val
			else: local_data[x] = np.nan

		# Get contamination
		for x in ['contamination','unreleased','alpha']:
			idx = df['type'] == x
			val = df.loc[idx,'value'].to_numpy()[0]
			valid_chk = all([df.loc[idx,c].to_numpy()[0] for c in ['replicated','converged']]) and not any([df.loc[idx,c].to_numpy()[0] for c in ['upper','lower']])
			if valid_chk: local_data[x] = val
			else: local_data[x] = np.nan

		# Get spawnGs
		idx_Gs = df['type'] == 'Gs'
		idx_spawn = df['type'] == 'spawn'
		valid_chk =  all([df.loc[idx_Gs,c].to_numpy()[0] for c in ['replicated','converged']]) and all([df.loc[idx_spawn,c].to_numpy()[0] for c in ['replicated','converged']])
		if valid_chk: local_data['spawnGs'] = df.loc[idx_Gs,'value'].to_numpy()[0] * df.loc[idx_spawn,'value'].to_numpy()[0]
		else: local_data['spawnGs'] = np.nan

		# Get first elongation rate
		idx = df['type'] == 'elongation'
		if local_data['strand'] == '+': first_idx = np.argsort(df.loc[idx,'start'].to_numpy())[0]
		else: first_idx = np.argsort(df.loc[idx,'stop'].to_numpy())[-1]
		val = df.loc[idx,'value'].to_numpy()[first_idx]
		valid_chk = all([df.loc[idx,c].to_numpy()[first_idx] for c in ['replicated','converged']]) and not any([df.loc[idx,c].to_numpy()[first_idx] for c in ['upper','lower']])
		if valid_chk: local_data['first_rate'] = val
		else: local_data['first_rate'] = np.nan

		# Get first cleavage rate
		idx = df['type'] == 'cleave'
		if local_data['strand'] == '+': first_idx = np.argsort(df.loc[idx,'start'].to_numpy())[0]
		else: first_idx = np.argsort(df.loc[idx,'stop'].to_numpy())[-1]
		val = df.loc[idx,'value'].to_numpy()[first_idx]
		valid_chk = all([df.loc[idx,c].to_numpy()[first_idx] for c in ['replicated','converged','valid_dependencies']]) and not any([df.loc[idx,c].to_numpy()[first_idx] for c in ['upper','lower']])
		if valid_chk: local_data['first_cleave'] = val
		else: local_data['first_cleave'] = np.nan

		# Get total cleavage rate
		val = np.sum(df.loc[idx,'value'].to_numpy())
		valid_chk = all([all(df.loc[idx,c].to_numpy()) for c in ['replicated','converged','valid_dependencies']]) and not any([any(df.loc[idx,c].to_numpy()) for c in ['upper','lower']])
		if valid_chk: local_data['total_cleave'] = val
		else: local_data['total_cleave'] = np.nan

		# Get splicing rates
		idx = df['type'] == 'splice'
		valid_chk = [True for _ in range(sum(idx))]
		if len(valid_chk) > 0:
			for c in ['replicated','converged','valid_dependencies']:
				temp = df.loc[idx,c].to_numpy()
				for i in range(len(valid_chk)): valid_chk[i] = valid_chk[i] and temp[i]
			for c in ['upper','lower']:
				temp = df.loc[idx,c].to_numpy()
				for i in range(len(valid_chk)): valid_chk[i] = valid_chk[i] and not temp[i]
			rates = df.loc[idx,'value'].to_numpy()[valid_chk]
		else: rates = np.array([])

		if len(rates) > 0:
			local_data['median_splice'] = np.median(rates)
			local_data['mean_splice'] = np.mean(rates)
			local_data['geomean_splice'] = 10**np.mean(np.log10(rates))
		else:
			for c in ['median_splice','mean_splice','geomean_splice']: local_data[c] = np.nan

		# Save data
		for x,y in local_data.items():
			data[x].append(y)
	return pd.DataFrame(data)

def getRegions(gene_df:pd.DataFrame,method:str,bw_path:Path) -> pd.DataFrame:
	# Get annotation
	pos = gene_df['strand'] == '+'
	gene_df['region_start'] = 0
	gene_df['region_stop'] = 0

	# Get regional positions
	if method == 'promoter1500':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start'] - 1000
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start'] + 500

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop'] - 500
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop'] + 1000
	elif method == 'promoter3000':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start'] - 2000
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start'] + 1000

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop'] - 1000
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop'] + 2000
	elif method == 'initial500':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start']
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start'] + 500

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop'] - 500
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop']
	elif method == 'initial1000':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start']
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start'] + 1000

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop'] - 1000
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop']
	elif method == 'initial1500':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start']
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start'] + 1500

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop'] - 1500
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop']
	elif method == 'initial2000':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start']
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start'] + 2000

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop'] - 2000
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop']
	elif method == 'upstream500':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start'] - 500
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start']

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop']
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop'] + 500
	elif method == 'upstream1000':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start'] - 1000
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start']

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop']
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop'] + 1000
	elif method == 'upstream1500':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start'] - 1500
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start']

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop']
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop'] + 1500
	elif method == 'upstream2000':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'start'] - 2000
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'start']

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'stop']
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'stop'] + 2000
	elif method == 'term1500':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'stop'] - 500
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'stop'] + 1000

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'start'] - 1000
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'start'] + 500
	elif method == 'term3000':
		gene_df.loc[pos,'region_start'] = gene_df.loc[pos,'stop'] - 1000
		gene_df.loc[pos,'region_stop'] = gene_df.loc[pos,'stop'] + 2000

		gene_df.loc[~pos,'region_start'] = gene_df.loc[~pos,'start'] - 2000
		gene_df.loc[~pos,'region_stop'] = gene_df.loc[~pos,'start'] + 1000
	elif method == 'base':
		gene_df['region_start'] = gene_df['start']
		gene_df['region_stop'] = gene_df['stop']

	# Set lower bound at 0
	gene_df.loc[gene_df['region_start'] < 0, 'region_start'] = 0

	# Set upper bound to chromosome length
	with pyBigWig.open(str(bw_path.absolute())) as bw:
		chrom_dict = bw.chroms()

		for c in gene_df['chromosome'].unique():
			max_length = chrom_dict[c]
			chrom_idx = gene_df['chromosome'] == c
			bound_idx = gene_df['region_stop'] > max_length
			idx = chrom_idx & bound_idx
			gene_df.loc[idx,'region_stop'] = max_length

	# Create output column
	col = bw_path.stem
	gene_df[col] = 0.0

	# Get bigwig coverage
	with pyBigWig.open(str(bw_path.absolute())) as bw:
		for i,row in gene_df.iterrows():
			gene_df.loc[i,col] = bw.stats(row['chromosome'],row['region_start'],row['region_stop'])

	gene_df = gene_df[['gene',col]]

	return gene_df