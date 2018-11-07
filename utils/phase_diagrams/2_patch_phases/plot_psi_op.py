
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from sys import argv
sns.set(style="white")
prog, file1, file2 = argv


data1=np.loadtxt(file1)
data2=np.loadtxt(file2)

fig, ax = plt.subplots()
plt.tick_params(axis='both', which='major', labelsize=24)

X=np.linspace(2.4,4.1,10)
Yerr1=[0.09, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.07,0.07]
Y=[ 0.4 for i in range(len(X))]
#plt.errorbar(x, y, xerr=0.2, yerr=0.4)
print data1[:,2]
plt.errorbar(data1[:,0], data1[:,1], yerr=data1[:,2], c='#FF9933', fmt="o", ms=10, capthick=2 )
plt.errorbar(data2[:,0], data2[:,1], yerr=Yerr1, c='#0099FF', lw=3, fmt="-o", ms=10, alpha=0.5, capthick=2)
plt.plot(X,Y, 'k--', alpha=0.4, lw=3)
plt.xlabel("$\\alpha$", size=40)
plt.ylabel("$\\psi$", size=40)
plt.xlim([2.4,4.1])
plt.ylim([-0.2,0.95])
plt.text(3.0, 0.6, r' ordered', fontsize=30, color="k", alpha=0.4)
plt.text(3.5, 0.1, r' random ', fontsize=30, color="k", alpha=0.4)
plt.tight_layout()
plt.savefig("psi_op_comparision.pdf")
plt.show()
