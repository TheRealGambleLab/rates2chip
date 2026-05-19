import logging
import argparse
import numpy as np

from pathlib import Path

from rates2chip.utilities import import_pandas, export_pandas


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--properties',help='path to rate properties')
	parser.add_argument('--null',help='path to null distribution table')
	parser.add_argument('--log',action='store_true',help='return log2 scalled values')
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
	prop_path = Path(args.properties)
	null_path = Path(args.null)
	out_path = Path(args.out)

	# Import dataframes
	prop_df = import_pandas(prop_path)
	null_df = import_pandas(null_path)

	# Normalize each gene by window size
	prop_df['_length'] = prop_df['stop'] - prop_df['start']

	# Allocate output dataframe
	id_columns = ['chromosome','start','stop','strand','gene','value','_length']
	out_df = prop_df[id_columns].copy()

	# Define stats columns
	stat_prefix = 'log_' if args.log else ''
	mean_col = f'{stat_prefix}mean'
	std_col = f'{stat_prefix}std'

	# Find marks to normalize
	marks:list[str] = [m for m in null_df['mark'].unique() if m in prop_df.columns]
	for mark in marks:
		logger.debug(mark)
		# Get null distribution for each mark
		mark_null = null_df[null_df['mark'] == mark].drop_duplicates(subset=['window_size']).sort_values('window_size')
		
		# Use largest window size for very large windows
		max_null_size = mark_null['window_size'].max()
		lookup_size = prop_df['_length'].clip(upper=max_null_size)

		# Get parameters for each window size
		mean_by_size = mark_null.set_index('window_size')[mean_col]
		std_by_size = mark_null.set_index('window_size')[std_col]

		# Map means and std onto existing values
		expected = lookup_size.map(mean_by_size).to_numpy(dtype=float)
		std = lookup_size.map(std_by_size).to_numpy(dtype=float)
		values = prop_df[mark].to_numpy(dtype=float)

		if args.log:
			# Convert to log scale based on log offset
			pc_by_size = mark_null.set_index('window_size')['log_offset']
			pc = lookup_size.map(pc_by_size).to_numpy(dtype=float)
			values = np.log2(values + pc)

		# Normalize
		invalid = ~np.isfinite(expected) | ~np.isfinite(std) | (std <= 0)
		z = (values - expected) / std
		z[invalid] = np.nan
		out_df[mark] = z

	out_df.drop(columns=['_length'],inplace=True)
	export_pandas(out_df,out_path)
	
