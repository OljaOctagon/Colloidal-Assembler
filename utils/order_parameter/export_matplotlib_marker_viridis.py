import matplotlib.pyplot as plt
from cycler import cycler 
import numpy as np

cmap=plt.cm.viridis
c = cycler('color', cmap(np.linspace(0,1,5)) )
colors = [item['color'] for item in c]
point = [0.5]
for i, color in enumerate(colors): 
    fig,ax = plt.subplots()
    plt.axis('off')
    plt.plot(point, c=color, lw=0, marker='s', ms=200)
    plt.savefig('marker_{}.png'.format(i))
    
