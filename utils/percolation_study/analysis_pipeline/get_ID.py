
import numpy as np 

ptype = ['double_manta_asymm_1','double_mouse_asymm_1', 'double_mouse_symm_1', 'double_mouse_symm_1']
delta = [0.2,0.3,0.4,0.5]
temperature = np.append([0.01,0.04,0.05,0.07],np.linspace(0.09,0.16,8))
phi = [0.01,0.05,0.075, 0.1, 0.125,0.15,0.175,0.2, 0.225,0.227,0.25, 0.3,0.35,0.375, 0.4, 0.425,0.45.0.475, 0.5, 0.525]

def id_gen():
	ID=-1
	for ptype_i in ptype:
		for delta_i in delta:
			for phi_i in phi:
				for temperature_i in temperature:
					ID = ID + 1 
					yield ptype_i, delta_i, phi_i, temperature_i, ID


def get_id():
	g = id_gen()

	IDS
	for item in g: 


def get_vars():



if __name__ == '__main__':

	# write id file 
	g = id_gen()
	with open("id_file.dat", 'w') as fh:
	for item in g:
		ptype, delta, phi, temperature, ID = item
		fh.write({},{},{},{},{}).format(ptype, delta, phi, temperature, ID)



