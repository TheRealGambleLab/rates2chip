import argparse

from pathlib import Path

from rates2chip.utilities import import_pandas
from rates2chip.plotting import histone_histograms

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='input rate path')
	parser.add_argument('--out',help='output path')
	parser.add_argument('--column',help='target column')
	args = parser.parse_args()

	df = import_pandas(Path(args.input))

	histone_histograms(df,args.column,Path(args.out))