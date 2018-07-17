import numpy
import os
from numpy import *
from random import randint
from matplotlib.pyplot import *
from matplotlib import gridspec
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

def tanh(z):
    return (exp(z) - exp(-z))/(exp(z) + exp(-z))

def find_velo(fl,array):
    nxx = len(array)
    for i in range(nxx):
        if array[i]>=fl:
            index = i-1
            break
    return index


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['axes.labelsize'] = 17
matplotlib.rcParams['xtick.labelsize'] = 15
matplotlib.rcParams['ytick.labelsize'] = 15
matplotlib.rcParams['legend.fontsize'] = 15
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

nx = 1000
ny = 50
bar = Bar('Processing Data', max=nx*ny, suffix='%(index)d/%(max)d - %(percent).1f%% - %(eta)ds')
x = linspace(0.1, 3.0,nx)
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
# for i in range(0,nx-1):
#     #print i
#     for j in range(0,ny-1):
#         m1 = 0
#         for k in range(len_arr1):
#             if arr1[k,3] > x[i]  and arr1[k,3]<= x[i+1]:
#                 if  arr1[k,1] >= y[j] and arr1[k,1]< y[j+1]:
#                     m1 += 1
#         density1[i,j] = m1
#         m2 = 0
#         for k in range(len_arr2):
#             if arr2[k,3] > x[i]  and arr2[k,3]<= x[i+1]:
#                 if  arr2[k,1] >= y[j] and arr2[k,1]< y[j+1]:
#                     m2 += 1
#         density2[i,j] = m2
#         #m3 = 0
#         for k in range(len_arr3):
#             if arr3[k,3] > x[i]  and arr3[k,3]<= x[i+1]:
#                 if  arr3[k,1] >= y[j] and arr3[k,1]< y[j+1]:
#                     m1 += 1
#         density1[i,j] = m1
#         m4 = 0
#         for k in range(len_arr4):
#             if arr4[k,3] > x[i]  and arr4[k,3]<= x[i+1]:
#                 if  arr4[k,1] >= y[j] and arr4[k,1]< y[j+1]:
#                     m4 += 1
#         density4[i,j] = m4
#         bar.next()
# bar.finish()
# figure(1,figsize = (10,5))
# subplot(1,2,1)
# density1[density1<0.5] = 'NaN'
# density1 = ma.masked_invalid(density1)
# cmap = get_cmap("Blues")
# pcolor(x,y,density1.T,cmap=cmap)
# xlabel(r"$ kd$")
# ylabel(r"$\mathrm{Froude\;Number}$")
# xlim(x[0],x[-1])
# #colorbar()
# subplot(1,2,2)
# density4[density4<0.5] = 'NaN'
# density4 = ma.masked_invalid(density4)
# cmap = get_cmap("Blues")
# pcolor(x,y,density4.T,cmap=cmap)
# xlabel(r"$ kd$")
# # ylabel(r"$\mathrm{Froude\;Number}$")
# xlim(x[0],x[-1])
# savefig('density.png',dpi=501, bbox_inches="tight")
#
# fig = figure(2,figsize = (10,5))
# gs = gridspec.GridSpec(2,2)
# fac1 = 0.9
# fac2 = fac1
# ax1= fig.add_subplot(gs[0,:])
# delta_s = 120.0
# fac3 = delta_s / pi
# ax1.plot(x,sqrt(fac1*x* tanh(fac2*x))**(-1.0)-0.04,'k-')
# ax1.plot(x, sqrt(tanh(x)/x),'k-')
# ax1.plot(x, sqrt(1/x),'k-')
# ax1.pcolor(x,y,density4.T,cmap=cmap)
# ax1.pcolor(x,y,density1.T,cmap=cmap)
# xlabel(r'$kd$')
# ylabel(r"$\textrm{Froude Number [']}$")
# #title_text = '$\delta_s = {test1:3.1f}$'.format(test1=delta_s)
# #print title_text
# #title(title_text,size =20)
# ylim(0.1,2.5)
# text(0.75,0.5,r'$\textrm{Ripples/Dunes}$', color='black',size=15)
# annotate(r'$\textrm{Transition Zone}$', xy=(0.5, 1.2), xytext=(1.25, 1.3),size =15,
#             arrowprops=dict(arrowstyle="->"),
#             )
#
# annotate(r'$\textrm{Antidunes}$', xy=(0.51, 1.6), xytext=(0.75, 1.9), size=15,
#             arrowprops=dict(arrowstyle="->"),
#             )
# text(2.0,1.9,r'$\textrm{Wash Out}$',size = 15)
# ax2= fig.add_subplot(gs[1,0])
# delta_s = 120.0
# fac3 = delta_s / pi
# ax2.plot(x*fac3,sqrt(fac1*x* tanh(fac2*x))**(-1.0)-0.04,'k-')
# ax2.plot(x*fac3, sqrt(tanh(x)/x),'k-')
# ax2.plot(x*fac3, sqrt(1/x),'k-')
# ax2.pcolor(x*fac3,y,density4.T,cmap=cmap)
# ax2.pcolor(x*fac3,y,density1.T,cmap=cmap)
# xlabel(r'$\textrm{Flow depth [m]} = \frac{kd \delta_s}{\pi}$')
# ylabel(r"$\textrm{Froude Number [']}$")
# title_text = '$\delta_s = {test1:3.1f} $'.format(test1=delta_s)
# print title_text
# title(title_text,size =20)
# xlim(3.75,20)
# ylim(0.1,2.5)
#
# ax3= fig.add_subplot(gs[1,1])
# delta_s = 50.0
# fac3 = delta_s / pi
# ax3.plot(x*fac3,sqrt(fac1*x* tanh(fac2*x))**(-1.0)-0.04,'k-')
# ax3.plot(x*fac3, sqrt(tanh(x)/x),'k-')
# ax3.plot(x*fac3, sqrt(1/x),'k-')
# ax3.pcolor(x*fac3,y,density4.T,cmap=cmap)
# ax3.pcolor(x*fac3,y,density1.T,cmap=cmap)
# xlabel(r'$\textrm{Flow depth [m]} = \frac{kd \delta_s}{\pi}$')
# ylabel(r"$\textrm{Froude Number [']}$")
# title_text = '$\delta_s = {test1:3.1f} $'.format(test1=delta_s)
# gs.update(wspace=0.25, hspace = 0.75)
# print title_text
# title(title_text,size =20)
# xlim(2,15)
# ylim(0.1,2.5)
# fig.savefig('FroudevsScoureD.png',dpi=501, bbox_inches="tight")
# #xlim(5,15)
# #figure(2,figsize = (10,5))
#


# figure(3,figsize = (10,5))
# delta_s = 11.0
# flow_depth =3.0
#
# fac1 = 0.9
# fac2 = fac1
# fac3 = delta_s / pi
# fd_index = find_velo(flow_depth,x*fac3)
#
# print fd_index
# cmap = get_cmap("Blues")
# density1[density1<0.5] = 'NaN'
# density1 = ma.masked_invalid(density1)
# plot(x*fac3,sqrt(fac1*x* tanh(fac2*x))**(-1.0)-0.04,'k-')
# plot(x*fac3, sqrt(tanh(x)/x),'k-')
# plot(x*fac3, sqrt(1/x),'k-')
# pcolor(x*fac3,y,density4.T,cmap=cmap)
# pcolor(x*fac3,y,density1.T,cmap=cmap)
# xlabel(r'$\textrm{Flow depth [m]} = \frac{kd \delta_s}{\pi}$')
# ylabel(r"$\textrm{Froude Number [']}$")
# title_text = '$\delta_s = {test1:3.1f} $'.format(test1=delta_s)
# text1=r"$\textrm{Dunes:}$" +" $< {test2:4.2f}\;m/s$".format(test2= sqrt(tanh(x[fd_index])/x[fd_index])*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# text2 = r"$ \textrm{Antidunes: }$"+"$ {test2:4.2f}\;m/s > u > {test3:4.2f} \;m/s$".format(test2= sqrt(1/x[fd_index])*sqrt(9.81*x[fd_index]*fac3),test3=(sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)-0.04)*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# print "Flow depth = {test1:3.1f}".format(test1=x[fd_index]*fac3)
# print "Flow speed:"
# print "\t Dunes: < {test2:4.2f} m/s".format(test2= sqrt(tanh(x[fd_index])/x[fd_index])*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# print "\t Antidunes: {test2:4.2f} m/s > u > {test3:4.2f}".format(test2= sqrt(1/x[fd_index])*sqrt(9.81*x[fd_index]*fac3),test3=(sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)-0.04)*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# #text(4,0.5,text1,size = 20)
# #text(6,1.5,text2,size = 20)
# plot([3.0,3.0],[0.1,sqrt(tanh(x[fd_index])/x[fd_index])],'r-')
# plot([3.0,3.0],[sqrt(1/x[fd_index]),sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)],'r-')
# annotate(text1, xy=(3., 0.5), xytext=(2.0, 0.65),size =20,
#             arrowprops=dict(arrowstyle="->"),
#             )
# annotate(text2, xy=(3., 1.3), xytext=(1., 2.0),size =20,
#             arrowprops=dict(arrowstyle="->"),
#             )
# #fd = u / sqrt(g h) -> u = fd * sqrt(g h)
# xlim(0.0,5)
# ylim(0.1,2.5)
# #savefig('flowspeed.png',dpi=501, bbox_inches="tight")
# show()

figure(3,figsize = (10,5))
delta_s = 30.0
flow_depth = linspace(0.5,10,100)

fac1 = 0.9
fac2 = fac1
fac3 = delta_s / pi
v_dune = []
v_antidune_l = []
v_antidune_h =[]
for i in range(100):
    fd_index = find_velo(flow_depth[i],x*fac3)
    v_dune.append(sqrt(tanh(x[fd_index])/x[fd_index])*sqrt(9.81*x[fd_index]*fac3))
    v_antidune_l.append(sqrt(1/x[fd_index])*sqrt(9.81*x[fd_index]*fac3))
    v_antidune_h.append((sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)-0.04)*sqrt(9.81*x[fd_index]*fac3))
plot(flow_depth,v_dune,'b-',alpha=0.5)
plot(flow_depth,v_antidune_l,'r-',alpha=0.5)
plot(flow_depth,v_antidune_h,'r-',alpha = 0.5)
fill_between(flow_depth, v_antidune_l,v_antidune_h,color='r',alpha=0.5)
fill_between(flow_depth, 0,v_dune,color='b',alpha=0.5)
xlim(0.5,5.0)
ylim(2,18)
xlabel(r'$\textrm{Flow depth,}\; m$')
ylabel(r'$\textrm{Flow Speed,}\; m/s$')
title_text = '$\delta_s = {test1:3.1f} $'.format(test1=delta_s)
title(title_text,size =20)
text(3.0,3.0,r'$\textrm{Dunes}$', size=20)
text(0.6,8.0,r'$\textrm{Antidunes}$', size=20)

# text1=r"$\textrm{Dunes:}$" +" $< {test2:4.2f}\;m/s$".format(test2= sqrt(tanh(x[fd_index])/x[fd_index])*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# text2 = r"$ \textrm{Antidunes: }$"+"$ {test2:4.2f}\;m/s > u > {test3:4.2f} \;m/s$".format(test2= sqrt(1/x[fd_index])*sqrt(9.81*x[fd_index]*fac3),test3=(sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)-0.04)*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# print "Flow depth = {test1:3.1f}".format(test1=x[fd_index]*fac3)
# print "Flow speed:"
# print "\t Dunes: < {test2:4.2f} m/s".format(test2= sqrt(tanh(x[fd_index])/x[fd_index])*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# print "\t Antidunes: {test2:4.2f} m/s > u > {test3:4.2f}".format(test2= sqrt(1/x[fd_index])*sqrt(9.81*x[fd_index]*fac3),test3=(sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)-0.04)*sqrt(9.81*x[fd_index]*fac3))
# title(title_text,size =20)
# #text(4,0.5,text1,size = 20)
# #text(6,1.5,text2,size = 20)
# plot([3.0,3.0],[0.1,sqrt(tanh(x[fd_index])/x[fd_index])],'r-')
# plot([3.0,3.0],[sqrt(1/x[fd_index]),sqrt(fac1*x[fd_index]* tanh(fac2*x[fd_index]))**(-1.0)],'r-')
# annotate(text1, xy=(3., 0.5), xytext=(2.0, 0.65),size =20,
#             arrowprops=dict(arrowstyle="->"),
#             )
# annotate(text2, xy=(3., 1.3), xytext=(1., 2.0),size =20,
#             arrowprops=dict(arrowstyle="->"),
#             )
# #fd = u / sqrt(g h) -> u = fd * sqrt(g h)
# xlim(0.0,5)
# ylim(0.1,2.5)
#savefig('flowspeed.png',dpi=501, bbox_inches="tight")
show()
