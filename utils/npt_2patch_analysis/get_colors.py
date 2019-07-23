import matplotlib.pyplot as plt 
import numpy as np 
cmap  = plt.cm.cubehelix(np.linspace(0,1,15))
x= np.arange(5)

for c,i in zip(cmap, range(1,16)):
    plt.plot(x,i*x, lw=3, c=c)

print(cmap)

plt.show()
