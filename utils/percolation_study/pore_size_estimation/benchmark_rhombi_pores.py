import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

# read data 
histos = []

lx=1
alpha = np.pi/3.
sigma = lx*np.sqrt(2+2*np.cos(alpha))

box_all = np.fromfile("Box.bin")
blx = box_all[3]
total_area = blx*blx 

phi_particles =0.125    
phi_all_pores = 1 - phi_particles 

area_rhombus = lx*lx * np.sin(alpha)

Nhisto=10
color="#391e3a"
for vs in range(1,Nhisto+1):
	a=pd.read_csv("histogram_vs_{}.dat".format(vs),header=None).values
	histos.append(a[:,0])


# pore packing bench mark 
fig,ax=plt.subplots()
arr=np.zeros((Nhisto,2))
for vs in range(1,Nhisto+1):
	arr[vs-1] = [vs, np.sum(histos[vs-1])/total_area]

plt.plot(sigma/arr[:,0], phi_all_pores*np.ones(10), lw=2, color='r', linestyle='--')
plt.plot(sigma/arr[:,0],arr[:,1], marker='o',c=color,lw=2)
plt.xlabel("$l_{p} / l_{r}$", size=15)
plt.ylabel("$\phi_{pores}$", size=15)
plt.savefig("pore_packing.pdf")

# large void area bench mark 
fig,ax=plt.subplots()
arr=np.zeros((Nhisto,2))
for vs in range(1,Nhisto+1):
	arr[vs-1] = [vs, np.max(histos[vs-1]/total_area)]

plt.plot(sigma/arr[:,0],arr[:,1], marker='o',c=color,lw=2)
plt.xlabel("$l_{p} / l_{r}$", size=15)
plt.ylabel("$\phi_{void}$", size=15)
plt.savefig("void_packing.pdf")

# histograms 
fig,ax=plt.subplots()
for vs in [2,5,7,10]:
	arr_p = np.sort(histos[vs-1])[:-1]
	hist, bin_edges = np.histogram(arr_p, bins=range(100),density=True)
	start=bin_edges[0]
	x=(bin_edges[:-1] + (bin_edges[1]-bin_edges[2])/2)/area_rhombus
	plt.plot(x,hist, lw=1, label='$l_p/l_r = {}$'.format(np.round(sigma/vs,3)))

print(area_rhombus)
plt.legend()
plt.xlabel("area pore  / area rhomus", size=15)
plt.ylabel("P", size=15)
plt.savefig("area_histograms.pdf")


# large void area bench mark 
fig,ax=plt.subplots()
arr  =np.zeros((Nhisto-1,3))
for vs in range(2,Nhisto+1):
	arr_p = np.sort(histos[vs-1])[:-1]/area_rhombus
	arr[vs-2] = [vs, np.mean(arr_p), np.std(arr_p)]

plt.errorbar(sigma/arr[:,0],arr[:,1],yerr=arr[:,2],capsize=10,c=color, lw=2)

plt.xlabel("$l_{p} / l_{r}$", size=15)
plt.ylabel("mean area pore / area rhombus", size=15)
plt.savefig("mean_pore_area.pdf")


# largest pore  
fig,ax=plt.subplots()
arr=np.zeros((Nhisto-1,2))
for vs in range(2,Nhisto+1):
	arr_p = np.sort(histos[vs-1])[:-1]
	arr[vs-2] = [vs, np.max(arr_p)]

plt.plot(sigma/arr[:,0],arr[:,1]/area_rhombus, marker='o',c=color,lw=2)
plt.xlabel("$l_{p} / l_{r}$", size=15)
plt.ylabel("area largest pore / area rhombus", size=15)
plt.savefig("largest_pore_area.pdf")

