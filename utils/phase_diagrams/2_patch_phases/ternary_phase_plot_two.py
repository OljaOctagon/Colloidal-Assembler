import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import colors 

def get_color(a,b,c,d):
			print(a,b,c,d)
			# 3 unique phases
			#pink 
			if c > 2/3.:
				cv = 7
			# purple 
			elif d > 2/3.:
				cv = 9
			#yellow
			elif b > 2/3.:
				cv = 8
			#blue 
			elif a > 2/3.:
				cv = 0
			# 6 mixtures
			# pink mixtures 
			# pink + yellow = salmon
			elif (c > 1/3. and b > 1/3.) or (d > 1/3. and e > 1/3.):
				cv = 2
			# pink + blue = pale purple
			elif (c > 1/3. and f > 1/3.)  or ( d>1/3. and a > 1/3.):
				cv = 3
			# pink + rest = pale pink
			elif (c > 1/3. and f < 1/3. and b < 1/3.) or ( d > 1/3. and a < 1/3. and e < 1/3.) :
				cv = 1
				print(c,d)
			# blue mixtures 
			# blue + yellow = pale green
			elif (f > 1/3. and b > 1/3.) or ( a>1/3. and e> 1/3.):
				cv = 4
			# blue + rest = pale blue 
			elif ( f > 1/3. and b < 1/3. and c < 1/3. ) or ( a <1/3. and e < 1/3. and d < 1/3. ):
				cv = 5
			# yellow mixtures 
			# yellow + rest = pale yellow
			elif ( b > 1/3. and c < 1/3.  and f < 1/3.) and (e > 1/3. and d < 1/3.  and a < 1/3.):
				cv = 6

			return cv 


df = pd.read_csv("cluster_distribution.csv")
df = df.fillna(value=0)

# value matrix for ternary color plot: 
# row dimension is number of energies
# rolumn dimension is number of positions
# the color vector is 3 dimensional: 
# 0: N_cluster < 3
# 1: N_cluster == 3 
# 2: N_cluster > 3 
print(df.patch_position)
Epsi = df.energy.unique()
N_Epsi = len(Epsi)
Pos = df.patch_position.unique()
N_Pos  = len(Pos)
vals = np.zeros((N_Epsi, N_Pos, 4))
zvals = np.zeros((N_Epsi, N_Pos))


for ei, nei in zip(Epsi, range(len(Epsi))):
	for pi, npi in zip(Pos, range(len(Pos))):
		df_sub = df[(df.energy == ei) & (df.patch_position == pi)] 

		a = df_sub.Nleq5.values[0]
		b = df_sub.Ngeq6.values[0]
		c = df_sub.N6.values[0]
		d = df_sub.N5.values[0]

		# bigger than 5 
		e =  df_sub.Ngeq6.values[0] + df_sub.N6.values[0]
		
		# smaller than 6 
		f =  df_sub.Nleq5.values[0] + df_sub.N5.values[0]

		print (ei, pi, a,b,c,d)
		zvals[nei,npi] = get_color(a,b,c,d)	
		vals[nei,npi] = np.array([a,b,c,d])

####

import matplotlib as mpl 
from matplotlib import pyplot
import numpy as np
from matplotlib.ticker import MultipleLocator

# make values from -5 to 5, for this example
#zvals = np.random.randint(9, size=(N_Epsi,N_Pos))

# make a color map of fixed colors
# 0: blue #236AB9
# 1: pale pink #ff80b3
# 2: salmon #ff9999
# 3: pale purple #ccccff
# 4: pale green: #92C591
# 5: pale blue: #67AFCB
# 6: pale yellow: #F1F791
# 7: pink #0099cc #ff0066
# 8: yellow #ffff00
# 9: purple "#9900cc"


color_list = ['#236AB9', '#ff80b3', '#ff9999',
'#ccccff', '#92C591', '#67AFCB', 
'#F1F791', '#ff0066', '#ffff00', "#9900cc"]

fig, ax = plt.subplots()

cmap = mpl.colors.ListedColormap(color_list)
#bounds=[-6,-2,2,6]
#bounds = [0,1,2,3,4,5,6,7,8,9,10]
bounds = np.arange(11) - 0.1 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# tell imshow about color map so that only set colors are used
img = pyplot.imshow(zvals,interpolation='nearest',
                    cmap = cmap,norm=norm)


x_og = np.arange(-.5,6.5,1)
y_og = np.arange(-.5,5.5,1)
#x_og = [0.5,1.5,2.5,3.5,4.5,5.5,6.5]
ax.set_xticks(x_og)
ax.set_yticks(y_og)


plt.grid(b=True,
	which='major', 
	axis='both', 
	linestyle='-', 
	color='k', linewidth=2)

xlabels=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
plt.xticks(x_og, xlabels, rotation='horizontal')
ylabels=[5.2,6.2,7.2,8.2,9.2,10.2]
plt.yticks(y_og, ylabels, rotation='horizontal')

plt.xlabel("$\\Delta$", size=25)
plt.ylabel("$\\epsilon [ k_{B}T ]$", size=25)


#from PIL import Image 

#size= 120,120
#im = Image.open('ternary_color_code_vs2.png')
#im.thumbnail(size, Image.ANTIALIAS)
#height = im.size[1]
#width = im.size[0]
#im = np.array(im).astype(np.float) / 255
#fig.figimage(im,fig.bbox.xmax-width, fig.bbox.ymax - height, zorder=10)
plt.tight_layout()
plt.savefig("star_phase_diagram_2.pdf")
plt.show()

#fig.bbox.xmax - height*2, fig.bbox.ymax - 2*height
