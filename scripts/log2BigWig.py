import logging
import argparse
import pyBigWig
import numpy as np

from pathlib import Path

def log2_bigwig(input_path:Path, output_path:Path, pseudocount:float) -> None:
	# Set up logging
	logger = logging.getLogger('__main__')
	with pyBigWig.open(str(input_path)) as input_bw:
		chromosomes = input_bw.chroms()
		with pyBigWig.open(str(output_path), "w") as output_bw:
			output_bw.addHeader(list(chromosomes.items()))
			for chromosome in chromosomes:
				# Collect bigwig intervals for chromosome
				intervals = input_bw.intervals(chromosome)
				if intervals is None:
					continue

                # Calculate log values
				starts: list[int] = []
				ends: list[int] = []
				values: list[float] = []
				for start, end, value in intervals:
					adjusted_value = value + pseudocount
					if adjusted_value <= 0:
						raise ValueError(f"non-positive value at {chromosome}:{start}-{end}")
					starts.append(start)
					ends.append(end)
					values.append(np.log2(adjusted_value))

                # Export logged values
				output_bw.addEntries([chromosome] * len(starts), starts, ends=ends, values=values)
				logger.debug('Completed chromosome %s',chromosome)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--bigwig", required=True, type=Path, help="input BigWig")
	parser.add_argument("--out", required=True, type=Path, help="output BigWig")
	parser.add_argument("--pseudocount",type=float,default=0.01,help="value added before log2 (default: 0; exact log2 transform)")
	parser.add_argument("--debug", action="store_true")
	args = parser.parse_args()

    # Set up logger
	logger = logging.getLogger('__main__')
	hdr = logging.StreamHandler()
	fmt = logging.Formatter('%(asctime)s\t%(message)s','%Y-%m-%d %H:%M:%S')
	hdr.setFormatter(fmt)
	logger.addHandler(hdr)
	if args.debug:logger.setLevel('DEBUG')
	else: logger.setLevel('INFO')
	logger.info("Writing %s", args.out)

    # Write bigwig
	log2_bigwig(args.bigwig, args.out, args.pseudocount)

