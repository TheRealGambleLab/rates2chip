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

def getWindows(rate_df:pd.DataFrame,pro_df:pd.DataFrame,tss_offset:int=500,cps_offset:int=0,min_win_size:int=100,min_win_time:float=5.0, min_win_pol:int=1, max_win:int|None = None, merge_remainder:bool=False) -> pd.DataFrame:
	# Group dataframes
	rate_group = rate_df.groupby('gene',observed=True)
	pro_group = pro_df.groupby('gene',observed=True)

	# Allocate output dataframe
	data = {'chromosome':[],'start':[],'stop':[],'strand':[],'gene':[],'value':[]}
	for gene,df in rate_group:
		# Get sequence properties
		row = df[df['type'] == 'rate'].iloc[0]
		gene_length = row['stop'] - row['start']

		# Confirm elongation rate is converged and valid
		if row['upper'] or row['lower'] or not row['converged'] or not row['replicated'] or not row['valid_dependencies']: continue
		elif gene_length <= tss_offset + cps_offset: continue
		elif not np.isfinite(row['value']) or row['value'] <= 0: continue

		# Get PROseq dataframe
		try: p_df = pro_group.get_group(gene)
		except KeyError: continue

		# Get TSS normalized polymerase counts
		if row['strand'] == '+': pro_pos = p_df['pos'].to_numpy() - row['start']
		else: pro_pos = row['stop'] - p_df['pos'].to_numpy()
		pro_val = p_df['value'].to_numpy()
		pro = np.zeros(gene_length + 1)
		pro_mask = np.logical_and(pro_pos >= 0,pro_pos <= gene_length)
		pro[pro_pos[pro_mask].astype(int)] = pro_val[pro_mask]
		if np.nansum(pro) <= 0: continue

		# Calculate time2nt with the same single-CPS setup used by optimizeTT.
		time2nt = tau(np.array([0,row['value']]),pro,np.array([(gene_length,0)]))

		# Get transcript-relative window breakpoints. Windows cover [tss_offset, gene_length - cps_offset) without leaving a short tail.
		window_start = tss_offset
		window_stop = gene_length - cps_offset
		if window_start >= window_stop: continue
		pairs = []
		window_rates = []
		for i in range(window_start + 1,window_stop + 1):
			# Window is too small
			if i - window_start < min_win_size: continue

			# Calculate time change between current position and window start
			delta = time2nt[i] - time2nt[window_start]

			# Check if has enough time or polymerases
			if delta <= 0 or delta < min_win_time: continue
			elif np.nansum(pro[window_start:i]) < min_win_pol: continue

			# Add window
			pairs.append((window_start,i))
			window_rates.append((i - window_start) / delta)
			
			# Begin new window
			window_start = i

		# Assign any terminal sequence that cannot satisfy the minimums
		# to the last valid window.
		if merge_remainder:
			if window_start < window_stop:
				if len(pairs) == 0: continue
				pairs[-1] = (pairs[-1][0],window_stop)
				delta = time2nt[pairs[-1][1]] - time2nt[pairs[-1][0]]
				if delta <= 0: continue
				window_rates[-1] = (pairs[-1][1] - pairs[-1][0]) / delta

		# Get absolute genomic positions
		if max_win:
			if len(pairs) > max_win:
				pairs = pairs[:max_win]
				window_rates = window_rates[:max_win]
		if len(pairs) == 0: continue
		pairs = np.array(pairs)
		if row['strand'] == '+': pairs = pairs + row['start']
		else: pairs = row['stop'] - pairs[:,[1,0]]

		# Add to data
		data['chromosome'].extend([row['chromosome'] for _ in range(len(pairs))])
		data['strand'].extend([row['strand'] for _ in range(len(pairs))])
		data['gene'].extend([row['gene'] for _ in range(len(pairs))])
		data['start'].extend(pairs[:,0])
		data['stop'].extend(pairs[:,1])
		data['value'].extend(window_rates)

	# Convert to dataframe
	data = pd.DataFrame(data)
	data.drop_duplicates(inplace=True)
	data.reset_index(inplace=True,drop=True)
	return data

def getRandomWindows(input_df:pd.DataFrame,seed:int,n_random:int,max_window_size:int) -> pd.DataFrame:
	rng = np.random.default_rng(seed)
	input_df = input_df.copy()
	input_df['window_size'] = input_df['stop'] - input_df['start']

	# Random intervals are drawn from gene spans represented in the input table.
	# This keeps the null distribution intragenic instead of drifting into intergenic sequence.
	gene_df = input_df.groupby(['gene','chromosome'],observed=True).agg(min_start=('start','min'),max_stop=('stop','max')).reset_index()
	gene_df['available'] = gene_df['max_stop'] - gene_df['min_start']

	window_sizes:np.ndarray = np.sort(input_df.loc[(input_df['window_size'] > 0) & (input_df['window_size'] <= max_window_size),'window_size'].unique())

	data = {'gene':[],'chromosome':[],'strand':[],'start':[],'stop':[],'window_size':[]}
	for window_size in window_sizes:
		valid_genes = gene_df[gene_df['available'] >= cast(float,window_size)].copy()
		if valid_genes.empty:
			continue

		valid_genes['n_start_positions'] = valid_genes['available'] - window_size + 1
		weights = valid_genes['n_start_positions'].to_numpy(dtype=float)
		weights = weights / weights.sum()
		gene_idx = rng.choice(valid_genes.index.to_numpy(),size=n_random,replace=True,p=weights)

		for idx in gene_idx:
			row = valid_genes.loc[idx]
			low = int(row['min_start'])
			high = int(row['max_stop'] - window_size)
			start = int(rng.integers(low,high + 1))
			stop = start + int(window_size)

			data['gene'].append(row['gene'])
			data['chromosome'].append(row['chromosome'])
			data['strand'].append('+')
			data['start'].append(start)
			data['stop'].append(stop)
			data['window_size'].append(window_size)

	return pd.DataFrame(data)