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

def tauC2rate(df:pd.DataFrame) -> pd.DataFrame:
	df['type'] = df['type'].astype('category')
	idx = df['type'] == 'tauC'
	df.loc[idx, 'value'] = (df.loc[idx,'stop'] - df.loc[idx,'start'])/df.loc[idx,'value']
	df['type'] = df['type'].cat.rename_categories({'tauC': 'rate'})
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