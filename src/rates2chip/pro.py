import numpy as np
import pandas as pd

from typing import cast

from optimizeTT.compiled import tau

def getPRO(df:pd.DataFrame,tss:int,stop:int|float=np.nan,pos_strand:bool=True) -> np.ndarray:
	# Get numpy arrays of data
	if pos_strand: pos = df['pos'].to_numpy() - tss
	else:pos = tss - df['pos'].to_numpy()
	val = df['value'].to_numpy()

	# Allocate output array
	if np.isnan(stop):
		pro = np.zeros(max(pos)+1)
		mask = pos > 0
	else:
		pro = np.zeros(cast(int,stop))
		mask = np.logical_and(pos > 0,pos < stop)
	
	# Assign coverage
	pro[pos[mask]] = val[mask]
	return pro

def getTime2nt(pro_df:pd.DataFrame,cleave_df:pd.DataFrame,tss:int,rate:float,alt_time:bool=False) -> np.ndarray:
	# Determine which strand the gene is on
	strand = cleave_df['strand'].unique()[0]

	# Find cleave positions
	if strand == '+':
		pos = True
		cleave_pos = cleave_df['stop'].to_numpy() - tss
	else:
		pos = False
		cleave_pos = tss - cleave_df['start'].to_numpy()

	# Get PRO coverage
	min_pro_length = int(max(cleave_pos)+1)
	pro = getPRO(pro_df,tss,min_pro_length,pos)
	
	# Sort cleavage positions and values
	idx = np.argsort(cleave_pos)
	cleave_pos = cleave_pos[idx]
	cleave_value = cleave_df['value'].to_numpy()[idx]

	# Calculate time to nucleotide
	if alt_time:
		rate = cleave_pos[0] / (np.sum(pro[:cleave_pos[0]]) * rate)

	time2nt = tau(np.array([0]+[i for i in cleave_value] + [rate]),pro,np.array([(x,i) for i,x in enumerate(cleave_pos)]))
	return time2nt

def addValidColumn(rate_df:pd.DataFrame,row_filter:str) -> pd.DataFrame:
	req_true_col = []
	req_false_col = []
	for x in row_filter.split(','):
		if x[0] == '~': req_false_col.append(x[1:])
		else: req_true_col.append(x)
	rate_df['valid'] = rate_df[req_true_col].all(axis=1) & ~rate_df[req_false_col].any(axis=1)
	return rate_df

def calculateProIntervals(rate_df:pd.DataFrame,pro_df:pd.DataFrame,use_elongation_index:bool=False,pause_mask:int=200,cps_mask:int=0,bin_size_metric:float=5,min_bin_nt:int=50,space_bins:bool=False)-> pd.DataFrame:
	# Check rate validity
	rate_df = addValidColumn(rate_df,'converged,replicated,~upper,~lower,valid_dependencies')
	
	# Group dfs by gene
	rate_group = rate_df.groupby('gene',observed=True)
	pro_group = pro_df.groupby('gene',observed=True)

	data = {x:[] for x in ['gene','chromosome','strand','start','stop','rate']}
	for gene,df in rate_group:
		row = df.loc[df['type'] == 'rate'].iloc[0]
		g = df.loc[df['type'] == 'g'].iloc[0]['value']
		gene_start = row['start']
		gene_stop = row['stop']
		# Skip rows that fail filters
		if not row['valid']: continue

		# Get PROseq dataframe
		try: p = pro_group.get_group(gene)
		except KeyError: continue
		
		# Get CPS dataframe
		cleave_df = df.loc[df['type'] == 'cleave']

		if row['strand'] == '+': tss = gene_start
		else: tss = gene_stop

		if use_elongation_index: time2nt = getTime2nt(p,cleave_df,tss,1/g,alt_time=True)
		else: time2nt = getTime2nt(p,cleave_df,tss,row['value'])

		# Drop masked regions
		if cps_mask > 0: time2nt = time2nt[pause_mask:-1*cps_mask]
		else: time2nt = time2nt[pause_mask:]

		# Skip genes smaller than minimum bin size
		if len(time2nt) < min_bin_nt: continue

		# Find bins
		if space_bins:
			# Calculate number of bins (at LEAST bin_size_metric large) required for gene
			num_bins = np.ceil(len(time2nt) // bin_size_metric)

			# Skip genes that are too short:
			if num_bins < 1: continue

			# Identify bin locations
			boundaries = np.linspace(0, len(time2nt)-1 , num_bins + 1)

			# Round to nearest integer
			boundaries = np.rint(boundaries).astype(int)

		else:
			b = [0]
			for i in range(len(time2nt)):
				# Skip bins too short (in nt or time)
				if (i - b[-1]) < min_bin_nt: continue
				elif (time2nt[i] - time2nt[b[-1]]) < bin_size_metric: continue

				# Check if next bin is invalid
				if (len(time2nt) - i) < min_bin_nt: break
				elif (time2nt[-1] - time2nt[i]) < bin_size_metric: break
				
				# Add index i to array b
				b.append(i)
			b.append(len(time2nt)-1)
			boundaries = np.array(b)
			
		# Save data
		for i in range(len(boundaries)-1):
			start = boundaries[i]
			stop = boundaries[i+1]
			delta = time2nt[stop] - time2nt[start]
			if delta == 0: rate = np.nan
			else: rate = (stop - start) / delta

			data['gene'].append(gene)
			data['chromosome'].append(row['chromosome'])
			data['strand'].append(row['strand'])
			if row['strand'] == '+':
				data['start'].append(start + gene_start + pause_mask)
				data['stop'].append(stop + gene_start + pause_mask)
			else:
				data['start'].append(gene_stop - pause_mask - stop)
				data['stop'].append(gene_stop - pause_mask - start)
			data['rate'].append(rate)
	return pd.DataFrame(data)