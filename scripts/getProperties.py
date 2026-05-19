import logging
import argparse

from pathlib import Path

from rates2chip.parsing import getCoverage
from rates2chip.utilities import import_pandas, export_pandas


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--coordinates',required=False,help='input rate path')
	parser.add_argument('--bigwig',nargs='+',default=[],help='path to bigwig files')
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

	coord_path = Path(args.coordinates)
	out_path = Path(args.out)

	# Get gene properties
	out_df = import_pandas(coord_path)

	# Get Bigwig
	for b in args.bigwig:
		bw_path = Path(b)
		logger.debug(bw_path.stem)
		out_df[bw_path.stem] = getCoverage(out_df[['gene','chromosome','strand','start','stop']].copy(),bw_path)
	export_pandas(out_df,out_path)
	
