import logging
import argparse
import numpy as np
import pandas as pd

from typing import cast
from pathlib import Path

from rates2chip.parsing import getCoverage
from rates2chip.utilities import import_pandas, export_pandas, subset_pandas


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--rates',help='input rate path')
	parser.add_argument('--bigwig',nargs='+',default=[],help='path to bigwig files')
	parser.add_argument('--filter',type=str,default='converged,replicated,valid_dependencies,~upper,~lower',help='coma separated list of boolean columns to filter by. use ~ to filter by negation of column')
	parser.add_argument('--upstream',default=500,help='distance upstream of tss')
	parser.add_argument('--downstream',default=0,help='distance downstream of tss')
	parser.add_argument('--out',help='output path (if output path exists, will merge onto the dataframe)')
	parser.add_argument('--debug',action='store_true')
	args = parser.parse_args()

	logger = logging.getLogger(__name__)
	hdr = logging.StreamHandler()
	fmt = logging.Formatter('%(asctime)s\t%(message)s','%Y-%m-%d %H:%M:%S')
	hdr.setFormatter(fmt)
	logger.addHandler(hdr)
	if args.debug: logger.setLevel('DEBUG')
	else: logger.setLevel('INFO')

	out_path = Path(args.out)

	# Get gene properties
	rate_df = subset_pandas(import_pandas(Path(args.rates)),cast(str,args.filter).split(','))
	
	# Get 'rate'
	if 'tauC' in rate_df['type'].unique():
		pass
	rate_df = rate_df[rate_df['type'] == 'rate']

	# Get TSS
	idx = rate_df['strand'] == '+'
	rate_df['region_start'] = np.where(idx,rate_df['start']-args.upstream,rate_df['stop']-args.downstream)
	rate_df['region_stop'] = np.where(idx,rate_df['start']+args.downstream,rate_df['stop']+args.upstream)
	rate_df = rate_df[['gene','chromosome','strand','region_start','region_stop']]
	rate_df.rename(columns={'region_start':'start','region_stop':'stop'},inplace=True)

	# Get Bigwig
	for b in args.bigwig:
		bw_path = Path(b)
		logger.debug(bw_path.stem)
		rate_df[f'tss_{bw_path.stem}'] = getCoverage(rate_df[['gene','chromosome','strand','start','stop']].copy(),bw_path)
	rate_df.drop(columns=['start','stop'],inplace=True)
	# Check if out_path exists before exporting.
	export_pandas(rate_df,out_path)