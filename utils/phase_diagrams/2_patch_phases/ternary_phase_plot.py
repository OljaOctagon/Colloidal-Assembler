import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import colors 
import matplotlib as mpl 
from matplotlib import pyplot
import numpy as np
from matplotlib.ticker import MultipleLocator
import argparse

def get_color_box(a,b,c):
	# 3 unique phases
	#pink
    if c > 2/3.:
        cv = 7
    #yellow
    elif b > 2/3.:
        cv = 8
    #blue
    elif a > 2/3.:
        cv = 0
    # 6 mixtures
    # pink mixtures
    # pink + yellow = salmon
    elif c > 1/3. and b > 1/3.:
        cv = 2
    # pink + blue = pale purple
    elif c > 1/3. and a > 1/3.:
       cv = 3
    # pink + rest = pale pink
    elif c > 1/3.  and a < 1/3. and b < 1/3.:
        cv = 1
    # blue mixtures
	# blue + yellow = pale green
    elif a > 1/3. and b > 1/3.:
        cv = 4
    # blue + rest = pale blue
    elif a > 1/3. and b < 1/3. and c < 1/3.:
        cv = 5
    # yellow mixtures
    # yellow + rest = pale yellow
    elif b > 1/3. and c < 1/3. and a < 1/3.:
        cv = 6

    return cv 

def get_color_star(a,b,c,d,e,f):
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

def star_vars(df_sub):
    a = df_sub.Nleq5.values[0]
    b = df_sub.Ngeq6.values[0]
    c = df_sub.N6.values[0]
    d = df_sub.N5.values[0]
    # bigger than 5
    e =  df_sub.Ngeq6.values[0] + df_sub.N6.values[0]
    # smaller than 6
    f =  df_sub.Nleq5.values[0] + df_sub.N5.values[0]

    return a,b,c,d,e,f

def box_vars(df_sub):
    a = df_sub.N1.values[0] + df_sub.N2.values[0]
    b = df_sub.N4.values[0] + df_sub.btN4.values[0]
    c = df_sub.N3.values[0]
    return a,b,c

def asymm_box_vars(df_sub):
    a = df_sub.Nl3.values[0] 
    b = df_sub.Nb3.values[0]
    c = df_sub.N3.values[0]
    return a,b,c

def ternary(df,Epsi,Pos,cluster_type):
    N_Epsi = len(Epsi)
    N_Pos = len(Pos)
    vals = np.zeros((N_Epsi, N_Pos, 4))
    zvals = np.zeros((N_Epsi, N_Pos))

    for ei, nei in zip(Epsi, range(N_Epsi)):
        for pi, npi in zip(Pos, range(N_Pos)):
            df_sub = df[(df.energy == ei) & (df.patch_position == pi)] 

            if cluster_type == 'asymm_box':
                a,b,c = asymm_box_vars(df_sub)
                d=0
                zvals[nei,npi] = get_color_box(a,b,c)

            if cluster_type == 'box':
                a,b,c = box_vars(df_sub)
                d=0
                zvals[nei,npi] = get_color_box(a,b,c)
            if cluster_type == 'star':
                a,b,c,d,e,f = star_vars(df_sub)
                zvals[nei,npi] = get_color_star(a,b,c,d,e,f)

            vals[nei,npi] = np.array([a,b,c,d])

    zvals = zvals[::-1,:]
    return zvals, vals

'''
Row dimension is number of energies
Column dimension is number of positions

The color vector is 3 dimensional: 
0: N_cluster < x
1: N_cluster == x 
2: N_cluster > x

colors:
0: blue #236AB9
1: pale pink #ff80b3
2: salmon #ff9999
3: pale purple #ccccff
4: pale green: #92C591
5: pale blue: #67AFCB
6: pale yellow: #F1F791
7: pink #0099cc #ff0066
8: yellow #ffff00
9: purple "#9900cc"

'''

#------------- parse file ----------

parser = argparse.ArgumentParser()
parser.add_argument("-type", help="available types: star, box, asymm_box",
                    type=str)

args = parser.parse_args()


#------------ read and process file 
df = pd.read_csv("cluster_distribution.csv")
df = df.fillna(value=0)

Epsi = df.energy.unique()
Pos = df.patch_position.unique()
zvals, vals = ternary(df, Epsi, Pos, args.type)

#------------ plot teneray phase diagram 

color_list = ['#236AB9', '#ff80b3', '#ff9999',
'#ccccff', '#92C591', '#67AFCB', 
'#F1F791', '#ff0066', '#ffff00', "#9900cc"]


fig, ax = plt.subplots()

cmap = mpl.colors.ListedColormap(color_list)
bounds = np.arange(11) - 0.1 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

img = pyplot.imshow(zvals,interpolation='nearest',
                    cmap = cmap,norm=norm)

x_og = np.arange(0,len(Pos),1)
y_og = np.arange(0,len(Epsi),1)
ax.set_xticks(x_og)
ax.set_yticks(y_og)

xlabels=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
ax.set_xticklabels(xlabels, minor=False)

ylabels=[5.2,6.2,7.2,8.2,9.2,10.2]
ax.set_yticklabels(ylabels[::-1], minor=False)

# minor ticks
ax.set_xticks(x_og-0.5, minor=True);
ax.set_yticks(y_og-0.5, minor=True);

ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())


plt.grid(b=True,
         which='minor',
         axis='both',
         linestyle='-',
         color='k', linewidth=1)

plt.xlabel("$\\Delta$", size=22)
plt.ylabel("$\\epsilon [ k_{B}T ]$", size=22)


ax.set_xticklabels(xlabels, minor=False)
ax.set_yticklabels(ylabels[::-1], minor=False)

ax.tick_params(axis='both',labelsize=16)

plt.tight_layout()
plt.savefig("phase_diagram_"+str(args.type)+".pdf")
plt.show()

