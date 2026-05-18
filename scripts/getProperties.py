import logging
import argparse
import numpy as np
import pandas as pd

from typing import cast
from pathlib import Path

from rates2chip.pro import getWindows
from rates2chip.parsing import getCoverage
from rates2chip.utilities import import_pandas, export_pandas, subset_pandas


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--rates',required=False,help='input rate path')
	parser.add_argument('--pro_path',required=False,help='path to PROseq dataframe')
	parser.add_argument('--bigwig',nargs='+',default=[],help='path to bigwig files')
	parser.add_argument('--filter',type=str,default='converged,replicated,valid_dependencies,~upper,~lower',help='coma separated list of boolean columns to filter by. use ~ to filter by negation of column')
	parser.add_argument('--method',default='',help='window method to use ()')
	parser.add_argument('--tss_offset',default=500,help='distance from tss to mask')
	parser.add_argument('--cps_offset',default=0,help='distance from CPS to mask')
	parser.add_argument('--min_len',default=50,help='minimum window length (nt)')
	parser.add_argument('--min_time',default=5,help='minimum window time (min)')
	parser.add_argument('--min_pro',default=1,help='minimum window polymerase count')
	parser.add_argument('--max_windows',required=False,help='maximum number of windows per gene')
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
		rate_df = subset_pandas(import_pandas(Path(args.rates)),cast(str,args.filter).split(','))
		if args.pro_path: 
			rate_df = rate_df[rate_df['type'] == 'rate']
			pro_df = import_pandas(Path(args.pro_path))
			out_df = getWindows(rate_df,pro_df,args.tss_offset,args.cps_offset,args.min_len,args.min_time,args.min_pro,args.max_windows)
		else: 
			out_df = rate_df[rate_df['type'] == 'elongation']
	
	else: out_df = import_pandas(out_path)

	# Get Bigwig
	for b in args.bigwig:
		bw_path = Path(b)
		logger.debug(bw_path.stem)
		out_df[bw_path.stem] = getCoverage(out_df[['gene','chromosome','strand','start','stop']].copy(),bw_path)
	export_pandas(out_df,out_path)
	
