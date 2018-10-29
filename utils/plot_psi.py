import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
sns.set(style='white')

# ncluster > 400
names = ['double manta', 'double mouse', 'whitlam et al', 'comparison']
marker=['o', 'X', 's', 'd']
color=['r', 'b', 'g', 'y']
x = np.array([[0.028080306781609193, 0.04677896872550968],
              [0.5183065888324874, 0.021978519622216603],
              [-0.082800232105,	 0.05],
              [0.024920415017857137, 0.035504098095303006]])

fig,ax = plt.subplots()
ax.set_aspect(0.99)
plt.ylim([-1,1])
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#ax.get_xaxis().set_visible(False)
#ax.yaxis.set_ticks_position('left')

# plot single compairsion lines at center position
for i in range(len(x)):
    ax.errorbar(0, x[i,0], yerr=x[i,1], color=color[i], elinewidth=3, capsize=2, label = names[i], marker=marker[i], ms=8)

def parse_json(filen):
    with open(filen) as fhandle:
        data = json.load(fhandle)

    pdata = [ value.append(key) for key, value in zip(data.keys(), data.values()) ]
    print(pdata)
    return data
# plot assym manta
ad_manta = parse_json("psi_mean_assym_manta.json")

# plot assym mouse
ad_mouse = parse_json("psi_mean_assym_mouse.json")

plt.legend(loc='best')
plt.savefig("whitlam_comparision.png", bbox_inches='tight', dpi=300)
plt.show()

n_points = np.array([264,591,-1, 84])

