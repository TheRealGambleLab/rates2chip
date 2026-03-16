import logging
import argparse
import numpy as np
import pandas as pd

from typing import cast
from pathlib import Path

from rates2chip.parsing import getRegions, getCoverage
from rates2chip.utilities import import_pandas, export_pandas, subset_pandas


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--rates',required=False,help='input rate path')
	parser.add_argument('--pro_path',required=False,help='path to PROseq dataframe')
	parser.add_argument('--bigwig',nargs='+',default=[],help='path to bigwig files')
	parser.add_argument('--filter',type=str,default='converged,replicated,valid_dependencies,~upper,~lower',help='coma separated list of boolean columns to filter by. use ~ to filter by negation of column')
	parser.add_argument('--method',default='',help='window method to use')
	parser.add_argument('--out',help='output path')
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
	if args.rates:
		rate_df = import_pandas(Path(args.rates))
		if args.pro_path: pro_df = import_pandas(Path(args.pro_path))
		else: pro_df = None
		out_df = getRegions(rate_df,pro_df,args.method,cast(str,args.filter).split(','))
	else: out_df = import_pandas(out_path)

	# Get Bigwig
	for b in args.bigwig:
		bw_path = Path(b)
		logger.debug(bw_path.stem)
		region_df = getCoverage(out_df[['gene','chromosome','strand','start','stop']].copy(),bw_path)
		out_df = pd.merge(out_df,region_df,how='outer',on=['gene','chromosome','strand','start','stop'])
	export_pandas(out_df,out_path)
	
