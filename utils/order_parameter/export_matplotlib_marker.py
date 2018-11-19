import matplotlib.pyplot as plt

markers = ['o','X','s','d','o','X','D','P','o']
colors = ['r','b','k','g','#ff6600', '#6600ff', 'y', '#33cc33', '#cc0066']

point = [0.5]
for color, marker in zip(colors, markers):
    fig,ax = plt.subplots()
    plt.axis('off')
    plt.plot(point, lw=0, marker=marker, c=color, ms=200)
    plt.savefig('marker_{}_{}.png'.format(marker, color))

