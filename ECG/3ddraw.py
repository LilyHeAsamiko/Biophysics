# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 09:15:40 2020

@author: LilyHeAsamiko
"""

from numpy import *
import pylab as p
#import matplotlib.axes3d as p3
import mpl_toolkits.mplot3d.axes3d as p3

# u and v are parametric variables.
# u is an array from 0 to 2*pi, with 100 elements
u=np.linspace(0, 100j, int(round(100/(2*np.pi))))
# v is an array from 0 to 2*pi, with 100 elements
v=np.linspace(0, 100j, int(round(100/np.pi)))
# x, y, and z are the coordinates of the points for plotting
# each is arranged in a 100x100 array
xu=10*np.outer(np.cos(u),np.sin(v))
yv=10*np.outer(np.sin(u),np.sin(v))
#zuv=10*np.outer(np.ones((np.size(u))),np.cos(v))


#!python
# this connects each of the points with lines
fig=plt.figure()
ax = mpl_toolkits.mplot3d.axes3d.Axes3D(fig)
# plot3D requires a 1D array for x, y, and z
# ravel() converts the 100x100 array into a 1x10000 array
ax.plot3D(np.ravel(xu),np.ravel(yv),np.ravel(zuv))
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
fig.add_axes(ax)
p.show()


#!python
fig=p.figure()
ax = p3.Axes3D(fig)
# scatter3D requires a 1D array for x, y, and z
# ravel() converts the 100x100 array into a 1x10000 array
ax.scatter3D(ravel(x),ravel(y),ravel(z))
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()

#!python
fig=p.figure()
ax = p3.Axes3D(fig)
# x, y, and z are 100x100 arrays
ax.plot_surface(x,y,z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()


#!python
delta = 0.025
x = arange(-3.0, 3.0, delta)
y = arange(-2.0, 2.0, delta)
X, Y = p.meshgrid(x, y)
Z1 = p.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = p.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = 10.0 * (Z2 - Z1)
fig=p.figure()
ax = p3.Axes3D(fig)
ax.contour3D(X,Y,Z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()

#!python
delta = 0.025
x = arange(-3.0, 3.0, delta)
y = arange(-2.0, 2.0, delta)
X, Y = p.meshgrid(x, y)
Z1 = p.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = p.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = 10.0 * (Z2 - Z1)
fig=p.figure()
ax = p3.Axes3D(fig)
ax.contour3D(X,Y,Z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()

#!python
x=r_[-10:10:100j]
y=r_[-10:10:100j]
z= add.outer(x*x, y*y)
### Contour plot of z = x**2 + y**2
p.contour(x,y,z)
### ContourF plot of z = x**2 + y**2
p.figure()
p.contourf(x,y,z)
p.show()


#!python
# note that for the following to work you have to modify the test funcitons in your site-packages/matplotlib/axes3d.py like this:
#def test_xxxx():
#    import pylab
#    ax = Axes3D(pylab.figure())
#    ....
#    ....
#    pylab.show()
# the following then work on 0.87.5
p3.test_bar2D()
p3.test_contour()
p3.test_scatter()
p3.test_scatter2D()
p3.test_surface()
# the following fail on 0.87.5
p3.test_plot()
p3.test_polys()
p3.test_wir