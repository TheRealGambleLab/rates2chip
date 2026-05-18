import argparse
import pandas as pd

from pathlib import Path

from rates2chip.utilities import import_pandas, export_pandas


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--x',help='dataframe #1')
	parser.add_argument('--y',help='dataframe #2')
	parser.add_argument('--columns',nargs='+',help='arguments to merge on')
	parser.add_argument('--out',help='output path')
	args = parser.parse_args()

	x_path = Path(args.x)
	y_path = Path(args.y)
	out_path = Path(args.out)

	x = import_pandas(x_path)
	y = import_pandas(y_path)

	z = pd.merge(x,y,how='inner',on=args.columns)

	export_pandas(z,out_path)