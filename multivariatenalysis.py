import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt_1
from numpy import random, histogram2d, diff
from scipy.interpolate import interp2d
#from scipy.interpolate import interpn
from scipy import interpolate
#from mipylib.numeric import interpolate
#from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import SmoothBivariateSpline


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

filename = "z-x-w.txt"

info = read_info(filename)
printInfo(info)

data = pd.read_csv(filename, skiprows=3, skip_blank_lines=False)

#fig, ax = plt.subplots(figsize =(10, 7))
H, xedges, yedges, image=plt.hist2d(data["z"], data["x"], weights=data["w"], bins =[35, 35])

#nbins = 35
#H, xedges, yedges = histogram2d(data["z"], data["x"] , bins =[35, 35], weights=data["w"])

def centers(edges):
    return edges[:-1] + diff(edges[:2])/2

xcenters = centers(xedges)
ycenters = centers(yedges)

# Construct interpolator
pdf = interp2d(xcenters, ycenters, H)

# define range of x values
x = np.arange(250, 950, 20)
# define range of z values
y = np.arange(-350, 350, 20)
# create a meshgrid
#xx, yy = np.meshgrid(x, y)
# these weights would be given to you in your histogram - here is an example function to generate weights
weights = pdf(x,y)
# create interpolation object
#f= interpolate.interp2d(x, y, weights, kind='cubic')

#f= RectBivariateSpline(x, y, weights)
f= SmoothBivariateSpline(x, y, weights)
# generate new ranges of x and z values
xnew = np.arange(250, 950, 20)
ynew = np.arange(-350, 350, 20)
# interpolate 
weightsnew = f(xnew, ynew)
# Adding color bar
#plt_1.colorbar()

print(weightsnew)
a_file = open("test.txt", "w")
for row in weightsnew:
    np.savetxt(a_file, row)

a_file.close()

prova= f(400,50)
print (prova)


plt_1.imshow(weightsnew)
plt_1.colorbar()
plt_1.show()

  
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




