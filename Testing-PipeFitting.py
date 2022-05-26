# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 10:19:48 2022

@author: M. McCrea

A simple example of using the Pipesfitting and path.py to re-reoute wires around a cylindrical pipe.
"""

import numpy as np
import matplotlib as plt

from patch import *
from PipesFitting import *

from datetime import datetime

#to apply a function to a list of points:
#print("checking point list in cube = ", np.apply_along_axis(PointInCube, 1, point_list, minimum = cmin, maximum=cmax))

print("RerouteSimplified::Coil Creation Start: " , datetime.now().strftime("%H:%M:%S"))
coils=coilset()
mypipes=pipelist()

#load pre-existing coil path files
# coils.input_csv("get_field-OutputFiles/CoilPathOutput-s0_005")
# coils.rotate(np.pi/2,0,0)
# print("RerouteSimplified::Coil Import End: " , datetime.now().strftime("%H:%M:%S"))

#minimal coil loop example
#side loops
loopa = np.array([[1.2,1.2,0.],[1.2,-1.2,0],[1,-1.2,0],[1,1.2,0],[1.2,1.2,0.]])
loopb = np.array([[-1.2,1.20,0.],[-1.2,-1.2,0],[-1,-1.2,0],[-1,1.20,0]])
loopc = loopa+np.array([0,0,-0.105])
loopd = loopb+np.array([0,0,0.1])
loope = loopa+np.array([0,0,-0.1])
loopf = loopb+np.array([0,0,-0.1])
loopg = loopb+np.array([0,0,0.5])
looph = loopa+np.array([0,0,-0.005])
loopi = loopa+np.array([0,0,+0.005])
loopj = loopa+np.array([0,0,+0.0075])

mypipes.add_pipe(0.101,0,0.0,'y','z',-1.21,1.21)

# #floor to ceiling loops
# loopz = np.array([[0,1.2,1.2],[0,-1.2,1.2],[0,-1.2,-1.2],[0,1.2,-1.2]])
# loopy = loopz+np.array([0+0.05,0,0])
# loopx = loopz+np.array([0+0.15,0,0])
# loopv = loopz+np.array([0+0.25,0,0])
# loopu = loopz+np.array([0+0.38,0,0])

coils.add_coil(loopa)
coils.add_coil(loopb)
coils.add_coil(loopc)
coils.add_coil(loopd)
coils.add_coil(loope)
coils.add_coil(loopf)
coils.add_coil(loopg)
coils.add_coil(looph)
coils.add_coil(loopi)
coils.add_coil(loopj)

# coils.add_coil(loopz)
# coils.add_coil(loopy)
# coils.add_coil(loopx)
# coils.add_coil(loopv)
# coils.add_coil(loopu)
print('Length %f'%coils.length())

fig = plt.figure(figsize=(9.5,9.5))
ax3=[]
ax3.append(fig.add_subplot(2, 2, 3)) #lower left, xy-plane
ax3.append(fig.add_subplot(2, 2, 4)) #lower right,yz-plane
ax3.append(fig.add_subplot(2, 2, 1)) #upper left, xz-plane
ax3.append(fig.add_subplot(2, 2, 2, projection='3d')) #upper right isometric view
# print(dir(coils))

#setting direction arrow and point position labels for all wire layout plots.
arrow=False   #adding arrows to show sgment directions can be slow for large sets
poslabel=False

print("Draw Original Layout Start:" , datetime.now().strftime("%H:%M:%S"))
coils.draw_layout(fig,arrow=arrow,poslabel=poslabel)
ax3[3].set_xlim3d(-1.4, 1.4)
ax3[3].set_ylim3d(-1.4, 1.4)
ax3[3].set_zlim3d(-1.4, 1.4)
print("Draw Original Layout End: " , datetime.now().strftime("%H:%M:%S"))
# plt.show()

pipe_density=14 #4# number of points to inscribe around the pipe

print("mypipes = ", mypipes)


print("ReRoutePipes Start: " , datetime.now().strftime("%H:%M:%S"))
mypipes.reroute_wires(coils,pipe_density=pipe_density)
print("ReRoutePipes End: " , datetime.now().strftime("%H:%M:%S"))

print("Draw rerouted Layout Start: " , datetime.now().strftime("%H:%M:%S"))
fig1 = plt.figure(figsize=(9.5,9.5))

ax4 = coils.draw_layout(fig1,arrow=arrow,poslabel=False)
ax4[3].set_xlim3d(-1.4, 1.4)
ax4[3].set_ylim3d(-1.4, 1.4)
ax4[3].set_zlim3d(-1.4, 1.4)

print("Draw Pipes Start: " , datetime.now().strftime("%H:%M:%S"))
mypipes.draw_layout(fig1,div_rad=40,color="b")
print("Draw Pipes End: " , datetime.now().strftime("%YYYY-%MM-%DD%H:%M:%S"))

fig.tight_layout()
plt.show()
