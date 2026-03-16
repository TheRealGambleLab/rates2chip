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
