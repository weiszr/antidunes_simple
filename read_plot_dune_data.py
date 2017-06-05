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


def read_data(fname):
#    print fname
    array1 = loadtxt(fname)
    return array1


def gauss(a,b,c,xx):
    return a * exp(-(xx-b)**2. / (2.0 * c**2.0))

fname1 = 'snorm_dunes.dat'
fname2 = 'santi_dunes.dat'
snd = read_data(fname1)
sad = read_data(fname2)
fname3 = 'unorm_dunes.dat'
fname4 = 'uanti_dunes.dat'
und = read_data(fname3)
uad = read_data(fname4)
#print shape(snd),shape(sad),shape(und),shape(uad)
figure(1)
# subplot(141)
# plot(nd[:,0],'ko')
# plot(ad[:,0],'ro')
#
# subplot(142)
# plot(nd[:,1],'ko')
# plot(ad[:,1],'ro')
#
# subplot(143)
x = linspace(0.5,3.0,100)
# #plot(x,0.57+gauss(5.,-1.,1.1,x))
#plot(x/(pi),0.2+(1/x)**0.88,'k-')
# plot(x/(pi),0.55+gauss(0.65,0.,1.15,x),'k-')
# plot(x/(pi),0.55+gauss(0.42,0.,1.18,x),'k-')
#for i in range(len(und))
scatter(und[:,3],und[:,1],s=2,color='yellow',alpha=1.0)
scatter(uad[:,3],uad[:,1],s=2,color='red',alpha=0.2)
scatter(snd[:,3],snd[:,1],s=2,color='blue',alpha=0.3)
scatter(sad[:,3]/pi,sad[:,1],s=2,color='green',alpha=0.5)
plot(x, sqrt(1/x),'k-')
plot(x, sqrt((x*tanh(x))**(-1.0)),'k-')
plot(x, sqrt(tanh(x)/x),'k-')
xlabel(r"Water/Scour Length")
ylabel("Froude Number")
xlim(x[0],x[-1])
text(1.7,1.6,"Wash out")
text(0.75,0.6,"Ripples/Dunes", color='white')
annotate('Transition Zone', xy=(0.3, 0.95), xytext=(0.4, 1.4),
            arrowprops=dict(arrowstyle="->"),
            )

annotate('Antidunes', xy=(0.6, 1.4), xytext=(0.3, 1.8),
            arrowprops=dict(arrowstyle="->"),
            )
savefig('PhasePlotBedforms.png',dpi=501, bbox_inches="tight")
# plot(snd[:,3],snd[:,1],'k-')
# plot(sad[:,3],sad[:,1],'k-')
#
#
# subplot(144)
# plot(nd[:,3],'ko')
# plot(ad[:,3],'ro')



# subplot(331)
# plot(nd[:,0],nd[:,1],'ko')
# plot(ad[:,0],ad[:,1],'ro')
#
# subplot(332)
# plot(nd[:,0],nd[:,2],'ko')
# plot(ad[:,0],ad[:,2],'ro')
# subplot(333)
# plot(nd[:,0],nd[:,3],'ko')
# plot(ad[:,0],ad[:,3],'ro')
# subplot(334)
# plot(nd[:,1],nd[:,2],'ko')
# plot(ad[:,1],ad[:,2],'ro')
#
# subplot(335)
# plot(nd[:,1],nd[:,3],'ko')
# plot(ad[:,1],ad[:,3],'ro')
#
# subplot(336)
# plot(nd[:,2],nd[:,3],'ko')
# plot(ad[:,2],ad[:,3],'ro')

# figure(2)
# plot(nd[:,1]/nd[:,2],nd[:,0]*nd[:,2],'ko')
# plot(ad[:,1]/ad[:,2],ad[:,0]*ad[:,2],'ro')

show()
