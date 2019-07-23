# calculate equation of state 
python calc_eqs.py
# plot equation of state
python plot_eqs.py
# calculate radial distribution function 
python gr_2D.py -dirs mu_0.3Energy_8.2Asymm_patchpos_0.3_Pressure_100_  -fname cluster_center_of_mass.dat

#asign color value
asign_color_to_cluster.py

# hexatic 
python hexatic_op.py 
