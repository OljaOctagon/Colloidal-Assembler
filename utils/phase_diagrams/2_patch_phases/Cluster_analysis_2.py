import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import seaborn as sns
from matplotlib import cm

Chemical_Potentials = [0.3, 0.4]
Energies = [ 5.2, 6.2 ,7.2 ,8.2 ]
Patch_Toplogies = ['symm', 'Asymm']
Patch_Positions = [0.2, 0.3 ,0.4 ,0.5, 0.6, 0.7, 0.8] 

max_length = 100
Nruns=1
N_vals = np.zeros(max_length)

NVAL = np.zeros((len(Energies),len(N_vals)))
N_vals_single = np.zeros((len(N_vals),Nruns))
N_VALS_std = np.zeros((len(Energies), len(N_vals)))
Multimers = np.zeros((len(Energies),max_length))

phi_mean=np.zeros(len(Energies))
phi_std = np.zeros(len(Energies))


f = open("cluster_distribution.csv", 'w')
f.write("chemical_potential,energy,phi_mean, phi_std, patch_topoloy,patch_position,N1,N2,N3,N4,>N4\n")

for patch_kind in Patch_Toplogies:

	if patch_kind == "symm":
		patch_kind_num = 0

	else:
		patch_kind_num = 1


	for patch_pos in Patch_Positions:
		for mu in Chemical_Potentials:
			epsi=0
			for energy in Energies:
				data=np.zeros(1)
				Phi=np.zeros(1)
				for nrun in xrange(1,Nruns+1):
					directory = "mu_"+str(mu)+"Energy_"+str(energy)+patch_kind+"_patchpos_"+str(patch_pos)+"_"+str(nrun)
					fname = directory+"/All_Clusters.dat"
					print fname
					df = pd.read_csv(fname, delim_whitespace=True, header=None)                                                                                                                                                    
					df.columns=['time','size']
					
					data_single = df['size'].values
					data = np.concatenate((data, data_single))
					
					d = np.diff(np.unique(data_single)).min()
					left_of_first_bin = data_single.min() - float(d)/2
					right_of_last_bin = data_single.max() + float(d)/2

					N_single, bins, patches = plt.hist(data_single, 
					np.arange(left_of_first_bin, right_of_last_bin + d, d),
					weights=data_single,
					alpha=0.5,
					normed = True)



					for nl in xrange(min(max_length-1,len(N_single))):
						N_vals_single[nl,nrun-1] = N_single[nl]

					N_vals_single[max_length-1,nrun-1] = np.sum(N_single[max_length-1:])

					fname2=directory+"/NPT_OUT.txt"
					dg = pd.read_csv(fname2, delim_whitespace=True ,header=None)
					dg.columns=['time', 'crap', 'N', 'crap2', 'phi']
					Phi=np.concatenate((Phi, dg['phi'].values))

					
				for nl in xrange(max_length):
					N_VALS_std[epsi,nl] = np.std(N_vals_single[nl])	


				data = np.delete(data,0)
				Phi = np.delete(Phi,0)

				phi_mean[epsi] = np.mean(Phi)
				phi_std[epsi] = np.std(Phi)

				d = np.diff(np.unique(data)).min()
				left_of_first_bin = data.min() - float(d)/2
				right_of_last_bin = data.max() + float(d)/2
				N_cluster, bins, patches = plt.hist(data, 
					np.arange(left_of_first_bin, right_of_last_bin + d, d),
					weights=data,
					alpha=0.5,
					normed = True)
				for nl in xrange(min(max_length,len(N_cluster))):
					N_vals[nl] = N_cluster[nl]

				N_vals[max_length-1] = np.sum(N_cluster[max_length-1:])
				NVAL[epsi] = N_vals
				state_vars = np.array([energy, phi_mean[epsi], phi_std[epsi], patch_kind_num, patch_pos])
				Output = np.insert(N_vals,0,state_vars)
				np.savetxt("cluster_data.csv", Output, delimiter=",",newline="\n")
				epsi=epsi+1


			phi_mean = np.around(phi_mean, decimals=3)
			phi_std  = np.around(phi_std, decimals=3)

			N_epsilon = len(Energies)
			
			for epsi  in  xrange(N_epsilon):
				for mer in xrange(max_length):
					Multimers[epsi,mer]= NVAL[epsi][mer]	 

			print "plot data"
			error_config = dict(ecolor='black', lw=2, capsize=6, capthick=2, alpha=0.6)

			sns.set(style='white')
			fig, ax =  plt.subplots(figsize=(15,11))
			sns.set(style='white')
			width = 1.
			ind = np.arange(N_epsilon) + 5.2 - width/2.
			plt.xticks(ind+width/2.)
			plt.xlim([5.2-width/2.,10.2+width/2.])
			plt.ylim((0,1))
			plt.tick_params(axis='both', which='major', labelsize=20)
			plt.xlabel("$\\epsilon [k_{B}T]$", size=30)
			plt.ylabel("$\\phi [\\%]$", size=30)

			plt.bar(ind, Multimers[:,0], width, alpha=0.7 )

			for epsi in xrange(N_epsilon):
				print np.argmax(Multimers[epsi]) + 1

			MER = np.zeros(len(Multimers[:,0]))

			for tm in xrange(1,max_length):
				MER = MER + Multimers[:,tm-1]
				#plt.bar(ind, Multimer[:,tm], width, yerr=N_VALS_std[:,tm], error_kw=error_config, bottom=MER, alpha=0.7)
				plt.bar(ind, Multimers[:,tm], width, color=cm.jet(1.*tm/len(ind)), bottom=MER, alpha=0.7)
			
			for i in xrange(len(ind)):
				ax.text(ind[i], 1.07, "$\\psi=$"+str(phi_mean[i]), fontsize=20)
				ax.text(ind[i], 1.03, "$\pm$"+str(phi_std[i]), fontsize=20)


			plt.savefig("mu_"+str(mu)+patch_kind+"_patchpos_"+str(patch_pos)+".pdf")

			######## Equilibrium length distribution

			fig, ax =  plt.subplots(figsize=(15,11))
			
			sns.set(style='white')
			plt.xlabel("$\\langle l \\rangle$", size=35)
			plt.ylabel("percentage of occurence ", size=30)
			plt.tick_params(axis='both', which='major', labelsize=25)

			epsi=0


			for energy in Energies:
				Xm = np.linspace(1,len(Multimers[epsi]),len(Multimers[epsi]))
				plt.plot(Xm, Multimers[epsi], lw=5, marker='o',ms="10", label="$\\epsilon = $"+str(energy)+"$k_{B}T$")
				lims = np.argmax(Multimers[epsi]) + 1
				plt.xlim((1,50))
				plt.legend(loc=1, prop={'size': 30})
			
				epsi += 1

			plt.savefig("equilibrium_length_mu_"+str(mu)+patch_kind+"_patchpos_"+str(patch_pos)+".pdf")
			
