import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("All_Clusters.dat", sep=" ", header=None)
df.columns=['time','size']
data = df['size'].values

d = np.diff(np.unique(data)).min()
left_of_first_bin = data.min() - float(d)/2
right_of_last_bin = data.max() + float(d)/2

n, bins, patches = plt.hist(data, 
	np.arange(left_of_first_bin, right_of_last_bin + d, d),  
	weights=data, 
	alpha=0.5,
	normed = True)
plt.savefig("histo.pdf")
print np.sum(n)
