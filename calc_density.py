import numpy
import os
from numpy import *
from random import randint
from matplotlib.pyplot import *
import matplotlib.patches as patches
import matplotlib.colors as colors
from pylab import *
#from mpl_toolkits.basemap import Basemap, cm
import string
import warnings
warnings.filterwarnings("ignore")
#import paramiko
import string
import webbrowser
import os
from progress.bar import Bar

def read_data(fname):
#    print fname
    array1 = loadtxt(fname)
    return array1

fname1 = 'snorm_dunes.dat'
fname2 = 'santi_dunes.dat'
snd = read_data(fname1)
sad = read_data(fname2)
fname3 = 'unorm_dunes.dat'
fname4 = 'uanti_dunes.dat'
und = read_data(fname3)
uad = read_data(fname4)


arr1 = snd
arr2 = sad
arr3 = und
arr4 = uad
len_arr1 = len(arr1[:,3])
len_arr2 = len(arr2[:,3])
len_arr3 = len(arr3[:,3])
len_arr4 = len(arr4[:,3])

nx = 100
ny = 50
bar = Bar('Processing Data', max=nx*ny, suffix='%(index)d/%(max)d - %(percent).1f%% - %(eta)ds')
x = linspace(0.5, 3.0,nx)
y = linspace(0.1,2.5,ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
density1 = (nx,ny)
density1 = zeros(density1)
density2 = (nx,ny)
density2 = zeros(density2)
density3 = (nx,ny)
density3 = zeros(density3)
density4 = (nx,ny)
density4 = zeros(density4)
for i in range(0,nx-1):
    #print i
    for j in range(0,ny-1):
        m1 = 0
        for k in range(len_arr1):
            if arr1[k,3] > x[i]  and arr1[k,3]<= x[i+1]:
                if  arr1[k,1] >= y[j] and arr1[k,1]< y[j+1]:
                    m1 += 1
        density1[i,j] = m1
        m2 = 0
        for k in range(len_arr2):
            if arr2[k,3] > x[i]  and arr2[k,3]<= x[i+1]:
                if  arr2[k,1] >= y[j] and arr2[k,1]< y[j+1]:
                    m2 += 1
        density2[i,j] = m2
        m3 = 0
        for k in range(len_arr3):
            if arr3[k,3] > x[i]  and arr3[k,3]<= x[i+1]:
                if  arr3[k,1] >= y[j] and arr3[k,1]< y[j+1]:
                    m3 += 1
        density3[i,j] = m3
        m4 = 0
        for k in range(len_arr4):
            if arr4[k,3] > x[i]  and arr4[k,3]<= x[i+1]:
                if  arr4[k,1] >= y[j] and arr4[k,1]< y[j+1]:
                    m4 += 1
        density4[i,j] = m4
        bar.next()
bar.finish()
figure(figsize = (5,20))
subplot(4,1,1)
pcolormesh(x,y,density1.T)

subplot(4,1,2)
pcolormesh(x,y,density2.T)

subplot(4,1,3)
pcolormesh(x,y,density3.T)

subplot(4,1,4)
pcolormesh(x,y,density4.T)
show()
