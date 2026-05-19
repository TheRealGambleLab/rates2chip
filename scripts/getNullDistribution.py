import logging
import argparse
import numpy as np
import pandas as pd

from pathlib import Path

from rates2chip.parsing import getCoverage
from rates2chip.pro import getRandomWindows
from rates2chip.utilities import import_pandas, export_pandas

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--coordinates',required=False,help='input rate path')
	parser.add_argument('--bigwig',nargs='+',default=[],help='path to bigwig files')
	parser.add_argument('--seed',type=int,default=42,help='seed for random permutations')
	parser.add_argument('--n_random',type=int,default=1000,help='number of random intervals per window size')
	parser.add_argument('--max_window_size',type=int,default=2000,help='maximum window size for random intervals')
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

	# Assign paths
	input_path = Path(args.coordinates)
	out_path = Path(args.out)

	# Get gene properties
	input_df = import_pandas(input_path)
	seed:int = args.seed

	# Get random positiions 
	window_df = getRandomWindows(input_df,seed,args.n_random,args.max_window_size)

	# Get indexes for each window size
	window_sizes = np.sort(window_df['window_size'].unique())
	window_idx = {window_size: np.where(window_df['window_size'].to_numpy() == window_size)[0] for window_size in window_sizes}

	# Get Bigwig
	data = {'window_size':[],'mark':[],'mean':[],'std':[],'n':[],'log_offset':[],'log_mean':[],'log_std':[],'log_n':[]}
	for b in args.bigwig:
		bw_path = Path(b)
		logger.debug(bw_path.stem)

		# Get mark intensity for random windows
		coverage = getCoverage(window_df[['gene','chromosome','strand','start','stop']].copy(),bw_path).to_numpy()
		pc = np.quantile(coverage[coverage>0],0.01)
	
		# Compute per window mark average and standard deviation for each window
		for window_size,idx in window_idx.items():
			window_coverage = coverage[idx].astype(float)
			log_cov = np.log2(window_coverage+pc)
			data['window_size'].append(window_size)
			data['mark'].append(bw_path.stem)
			data['mean'].append(np.nanmean(window_coverage))
			data['std'].append(np.nanstd(window_coverage,ddof=1))
			data['n'].append(np.sum(np.isfinite(window_coverage)))
			data['log_offset'].append(pc)
			data['log_mean'].append(np.nanmean(log_cov))
			data['log_std'].append(np.nanstd(log_cov))
			data['log_n'].append(np.sum(np.isfinite(log_cov)))

	export_pandas(pd.DataFrame(data),out_path)
