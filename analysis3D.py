import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt_1
from numpy import random, histogramdd, diff
from scipy.interpolate import interpn
#from scipy.interpolate import interpn
#from scipy import interpolate
#from mipylib.numeric import interpolate
#from scipy.interpolate import RectBivariateSpline
#from scipy.interpolate import SmoothBivariateSpline
#import scipy as sp
#import scipy.interpolate
from scipy.interpolate import Rbf



#atplotlib.use('TkAgg')

def read_info(path):
    d = {}
    with open(path) as myFile:
        line1 = next(myFile)
        l1 = line1.strip().split(",")
        line2 = next(myFile)
        l2 = line2.strip().split(",")
        for name, elem in zip(l1, l2):
            try:
                d[name] = int(elem)
            except Exception:
                try:
                    d[name] = float(elem)
                except Exception:
                    d[name] = elem
    return d

def printInfo(info):
    print("Coupling constant alpha_s\t = ", info["alpha_s"])
    print("Mass of top quark\t\t = ", info["top_quark_mass"], "GeV")
    print("Energy of single proton\t\t = ", info["energy"], "GeV")

filename = "lol.txt"

info = read_info(filename)
printInfo(info)

data = pd.read_csv(filename, skiprows=3, skip_blank_lines=False)

x=data["x"].astype(np.float32)
y=data["y"].astype(np.float32)
z=data["z"].astype(np.float32)
w=data["w"].astype(np.float32)

xnew = np.arange(-350, 350, 20)
ynew = np.arange(-350, 350, 20)
znew = np.arange(200, 900, 20)

spline = Rbf(x, y, z, w,function='thin_plate',smooth=5, epsilon=5)
weightsnew =spline (xnew,ynew,znew)

pd.DataFrame({'col1':weightsnew,'col2':xnew,'col3':xnew}).to_csv('out.txt', index=False, header=None, sep='\t')

#fig, ax = plt.subplots(figsize =(10, 7))
#H, xedges, yedges, zedges, image=plt.hist3d(data["x"], data["y"], data["z"], weights=data["w"], bins =[35, 35, 35])

#nbins = 35
#H, [xedges, yedges, zedges] = histogramdd((data["x"], data["y"], data["z"]) , bins=35, weights=data["w"])

#def centers(edges):
 #  return edges[:-1] + diff(edges[:2])/2

#xcenters = centers(xedges)
#ycenters = centers(yedges)
#zcenters = centers(zedges)

#x = np.arange(-350, 350, 20)
#y = np.arange(-350, 350, 20)
#z = np.arange(200, 900, 20)
#xx, yy, zz = np.meshgrid(xcenters, ycenters, zcenters)

#weights = interpn((xx, yy, zz), H, (x,y,z))

#xx, yy, zz = np.meshgrid(xcenters, ycenters, zcenters)
#print (xx)
#pd.DataFrame({'col1':xx,'col2':yy,'col3':zz}).to_csv('out.txt', index=False, header=None, sep='\t')
# Construct interpolator
#pdf = interpn((xcenters, ycenters, zcenters), H, interpolate_grid)

# define range of x values

# define range of z values

#x = np.arange(-350, 350, 20)
#y = np.arange(-350, 350, 20)
#z = np.arange(200, 900, 20)
#interpolate_grid = np.array(x,y,z)
# create a meshgrid
#xx, yy = np.meshgrid(x, y)
# these weights would be given to you in your histogram - here is an example function to generate weights
#pdf = interpn((xcenters, ycenters, zcenters), H, (x,y,z),'cubic')
#weights = interpn((xcenters, ycenters, zcenters), H, (x,y,z))
# create interpolation object
#f= interpolate.interp2d(x, y, weights, kind='cubic')
#xnew = np.arange(-350, 350, 20)
#ynew = np.arange(-350, 350, 20)
#znew = np.arange(200, 900, 20)
#f= RectBivariateSpline(x, y, weights)

#sparse_points = np.array([data["x"], data["y"], data["z"]]) 
#dense_points = np.array([xnew, ynew, znew]) 

#spline = Rbf(xcenters, ycenters, zcenters,H,function='thin_plate',smooth=5, epsilon=5)
#spline = Rbf(x, y, z, xnew,function='thin_plate',smooth=5, epsilon=5)
#weightsnew =spline (xnew,ynew,znew)
#zfun_smooth_rbf = RBFInterpolator((x,y), z, smoothing=0, kernel='thin_plate_spline')  

#weightsnew = RBFInterpolator((data["x"], data["y"], data["z"]),data["w"],smoothing=0.0, kernel='thin_plate_spline') 

# generate new ranges of x and z values

# interpolate 
#weightsnew = zfun_smooth_rbf(dense_points)
# Adding color bar
#plt_1.colorbar()


#print((weights))





#plt_1.imshow(weightsnew)
#plt_1.colorbar()
#plt_1.show()

  
#ax.set_xlabel('Z-axis') 
#ax.set_ylabel('X-axis') 
  
# show plot
#plt.tight_layout() 
#plt.show()

##########spline##########

# define range of x values
#x = np.arange(0, 1350, 38)
# define range of z values
#y = np.arange(-350, 350, 20)
# create a meshgrid
#xx, yy = np.meshgrid(x, y)
# these weights would be given to you in your histogram - here is an example function to generate weights
#weights = np.sin(xx**2+zz**2)
# create interpolation object
#interpolate.interp2d(x, y, data["w"], kind='cubic')


# generate new ranges of x and z values
#x = np.arange(0, 1350, 38)
#y = np.arange(-350, 350, 20)
# interpolate 
#weightsnew = f(xnew, ynew)

# plot
#plt_1.imshow(weightsnew)
#plt_1.show()




