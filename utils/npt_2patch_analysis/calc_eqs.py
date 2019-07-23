import numpy as np 
import pandas as pd 

Pressures = [5,6,8,10,11,12, 14, 16,17, 18, 19, 20, 25, 30, 35, 40, 50, 60, 80,100]
patch_positions = [0.3]
Nruns = 8 
phi = [ [] for i in range(Nruns) ]

filetrace = 'mu_0.3Energy_8.2Asymm_patchpos_'
for pos in patch_positions:
    eqs = []
    for p in Pressures:
        for nrun in range(1, Nruns+1):
            run_trace = filetrace+str(pos)+"_Pressure_"+str(p)+"_"+str(nrun)
            file = "/NPT_OUT.txt"
            dfile = run_trace+file
            phi = pd.read_csv(dfile, header=None, delim_whitespace=True).values[:,4]

            '''
            file = "/cluster_center_of_mass.dat"
            dfile = run_trace+file
            arr = pd.read_csv(dfile, header=None, delim_whitespace=True).values[:,1:4]
           
            arr_unique = np.unique(arr, axis=0)
            
            Np = arr_unique[-30:,0]
            Lx = arr_unique[-30:,1]
            Ly = arr_unique[-30:,2]

            area = Lx*Ly 

            # vol of box
            rlx = 1.
            alpha= (60*np.pi)/180.0;
            rh = rlx * np.sin(alpha)
            Vp= 3*rlx*rh 
       
            # packing fraction 
            phi[nrun-1] = (Np*Vp)/area
            '''

        sphi = np.array(phi).flatten()
        print(sphi)
        mean_phi = np.mean(sphi)
        std_phi = np.std(sphi)
        eqs.append([p, mean_phi, std_phi])

    np.savetxt('eqs_'+str(pos)+'.dat', eqs, fmt='%.4f')




