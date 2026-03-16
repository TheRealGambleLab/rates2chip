import pyBigWig
import numpy as np
import pandas as pd

from pathlib import Path

from .utilities import subset_pandas
from .pro import calculateProIntervals

def getRegions(gene_df:pd.DataFrame,pro_df:pd.DataFrame|None=None,method:str='',filter_list:list[str]=['converged','replicated','valid_dependencies','~upper','~lower']) -> pd.DataFrame:
	local_df = gene_df.copy()
	if method == 'splice_total':
		pass
	else: # Elongation rate

		# Subset for only genes that pass threshold
		if pro_df is not None:
			local_df = calculateProIntervals(local_df,pro_df,False,pause_mask=500)
		else:
			local_df = subset_pandas(gene_df,filter_list)
			local_df = local_df[local_df['type'] == 'elongation']
	
	return local_df

def getCoverage(gene_df:pd.DataFrame, bw_path:Path):
	# Create output column
	col = bw_path.stem
	gene_df[col] = 0.0

	# Get bigwig coverage
	with pyBigWig.open(str(bw_path.absolute())) as bw:
		for i,row in gene_df.iterrows():
			gene_df.loc[i,col] = bw.stats(row['chromosome'],row['start'],row['stop'])

	gene_df = gene_df[['gene',col]]

	return gene_df