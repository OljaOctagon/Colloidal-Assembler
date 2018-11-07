import numpy as np

class disk:
	def __init__(self, x,y, radius):
		self.r = radius
		self.x = x
		self.y = y
		self.vol = self.r*self.r*np.pi

disk1 = disk(0,0,1)

class box:
	def __init__(self, x,y, length):
		self.x = x
		self.y = y
		self.l = length
		self.vol = self.l*self.l

box1 = box(0,0,10)

# Brownian Motion with periodic boundaries
T=100
x = np.zeros(T)
y = np.zeros(T)

for time in range(1,T):
	disk1.x = disk1.x + np.random.normal(0,0.5)
	disk1.x = disk1.x - box1.x*np.rint(disk1.x/box1.l)
	x[time] = disk1.x

	disk1.y = disk1.y + np.random.normal(0,0.5)
	disk1.y = disk1.y - box1.y*np.rint(disk1.y/box1.l)
	y[time] = disk1.y


import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go

particle_shape = {
	'type': 'circle',
    'x0': x[time] - disk1.r,
    'y0': y[time] - disk1.r,
    'x1': x[time] + disk1.r, 
    'y1': y[time] + disk1.r,
    'line': {
        'color': 'rgba(0,0,0,1)',
        'width': 3
    },
}

_shapes = [particle_shape]

trace0 = go.Scatter(
	x=[x[0]],
	y=[y[0]]
)
data = [trace0]

layout=dict(xaxis=dict(range=[-12, 12], autorange=False, zeroline=False),
            yaxis=dict(range=[-12, 12], autorange=False, zeroline=False),
            title='Brownian Motion of a disk', hovermode='closest',
            shapes = _shapes,
            updatemenus= [{'type': 'buttons',
                           'buttons': [{'label': 'Play',
                                        'method': 'animate',
                                        'args': [None, {'frame': {'duration': 0, 'redraw': True},
                                        'transition': {'duration': 0 }}]
                                        }]}])

frames=[dict(data=[dict(x=[x[time]], 
                        y=[y[time]], 
                        mode='markers', 
                        marker=dict(color='red', size=50),
                        duration = 0
                        )
                  ]) for time in range(T)]    
          
figure1=dict(data=data, layout=layout, frames=frames)          

plotly.offline.plot(figure1, filename='brownian.html')
