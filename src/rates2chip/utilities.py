import pandas as pd

from pathlib import Path

def import_pandas(path:Path) -> pd.DataFrame:
	if path.suffix == '.parquet': dataframe = pd.read_parquet(path)
	elif path.suffix == '.obj' or path.suffix == '.pkl': dataframe = pd.read_pickle(path)
	elif path.suffix == '.csv': dataframe = pd.read_csv(path)
	else: dataframe = pd.read_csv(path,sep='\t')
	return dataframe

def export_pandas(dataframe:pd.DataFrame,out_path:Path) -> None:
	if out_path.suffix == '.parquet': dataframe.to_parquet(out_path,index=False)
	elif out_path.suffix == '.obj' or out_path.suffix == '.pkl': dataframe.to_pickle(out_path)
	elif out_path.suffix == '.csv': dataframe.to_csv(out_path,index=False)
	else: dataframe.to_csv(out_path,index=False,sep='\t')

def subset_pandas(dataframe:pd.DataFrame,subset:list[str]) -> pd.DataFrame:
	temp = dataframe
	for x in subset: 
		if x[0] == '~': temp = temp[~temp[x[1:]]]
		else: temp = temp[temp[x]]
	return temp

def trimTssWindows(rate_df:pd.DataFrame,tss_offset:int,cps_offset:int=0) -> pd.DataFrame:
	"""Trim the SKaTER elongation windows nearest the TSS and CPS."""
	windows = rate_df[rate_df['type'] == 'elongation'].copy()
	gene_rates = rate_df[rate_df['type'] == 'rate']
	drop = set()
	for gene,gene_windows in windows.groupby('gene',observed=True):
		gene_rate = gene_rates[gene_rates['gene'] == gene].iloc[0]
		if gene_rate['strand'] == '+':
			tss_idx = gene_windows['start'].idxmin()
			mask_stop = gene_rate['start'] + tss_offset
			if windows.loc[tss_idx,'stop'] <= mask_stop: drop.add(tss_idx)
			elif windows.loc[tss_idx,'start'] < mask_stop: windows.loc[tss_idx,'start'] = mask_stop

			cps_idx = gene_windows['start'].idxmax()
			mask_start = gene_rate['stop'] - cps_offset
			if windows.loc[cps_idx,'start'] >= mask_start: drop.add(cps_idx)
			elif windows.loc[cps_idx,'stop'] > mask_start: windows.loc[cps_idx,'stop'] = mask_start
		else:
			tss_idx = gene_windows['start'].idxmax()
			mask_start = gene_rate['stop'] - tss_offset
			if windows.loc[tss_idx,'start'] >= mask_start: drop.add(tss_idx)
			elif windows.loc[tss_idx,'stop'] > mask_start: windows.loc[tss_idx,'stop'] = mask_start

			cps_idx = gene_windows['start'].idxmin()
			mask_stop = gene_rate['start'] + cps_offset
			if windows.loc[cps_idx,'stop'] <= mask_stop: drop.add(cps_idx)
			elif windows.loc[cps_idx,'start'] < mask_stop: windows.loc[cps_idx,'start'] = mask_stop

	columns = ['chromosome','start','stop','strand','gene','value']
	return windows.drop(index=list(drop)).reset_index(drop=True)[columns]

def getSkaterWindows(rate_df:pd.DataFrame,filters:list[str],tss_offset:int,cps_offset:int) -> pd.DataFrame:
	gene_rates = rate_df[rate_df['type'] == 'rate']
	windows = subset_pandas(rate_df[rate_df['type'] == 'elongation'],filters)
	return trimTssWindows(pd.concat([gene_rates,windows]),tss_offset,cps_offset)

def tauC2rate(df:pd.DataFrame) -> pd.DataFrame:
	df = df.copy()
	idx = df['type'] == 'tauC'
	df.loc[idx, 'value'] = (df.loc[idx,'stop'] - df.loc[idx,'start'])/df.loc[idx,'value']
	df['type'] = df['type'].astype('object')
	df.loc[idx, 'type'] = 'rate'
	df['type'] = df['type'].astype('category')
	return df

def filter_df(input_df:pd.DataFrame,x:str) -> pd.DataFrame:
	df = input_df.copy()
	if '==' in x:
		df = df[df[x.split('==')[0]] == float(x.split('==')[1])]
	elif '<=' in x:
		df = df[df[x.split('<=')[0]] <= float(x.split('<=')[1])]
	elif '>=' in x:
		df = df[df[x.split('>=')[0]] >= float(x.split('>=')[1])]
	elif '!=' in x:
		df = df[df[x.split('!=')[0]] != float(x.split('!=')[1])]
	elif '>' in x:
		df = df[df[x.split('>')[0]] > float(x.split('>')[1])]
	elif '<' in x:
		df = df[df[x.split('<')[0]] < float(x.split('<')[1])]
	else: raise Exception('Invalid filter')
	return df
