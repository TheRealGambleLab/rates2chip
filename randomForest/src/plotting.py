import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

from pathlib import Path

def plot_bar(df:pd.DataFrame,x_col:str,y_col:str,hue_col:str,out_path:Path) -> None:
	plt.figure(figsize=(10,4))
	sb.barplot(df,x=x_col,y=y_col,hue=hue_col)
	plt.xticks(rotation=45, ha='right')
	#plt.legend(loc='lower left')
	plt.tight_layout()
	plt.savefig(out_path,dpi=200)
	plt.close()

def histone_histograms(df:pd.DataFrame,data_column:str,out_path:Path) -> None:
	histones = ['H2A','H2B','H3','H4']
	df['Histone'] = ''
	for h in histones:
		idx = df['Feature'].str.contains(h)
		df.loc[idx,'Histone'] = h

	df.sort_values(by=['Histone','Feature'], inplace=True)

	plot_bar(df,'Feature',data_column,'Histone',out_path)