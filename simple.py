import pylab
from pylab import *
from numpy import *
import matplotlib

### Functions #####

####################fall velocity##################################
def ferg(se):
    rhos=2650.0
    rho=1000.0
    g=9.81
    nu=1.0004E-6
    gamma=g*(rhos-rho)
    r=gamma/(rho*g)
    c1=18.0
    c2=0.4
    ws=r*g*se**2.0/(c1*nu+sqrt(0.75*c2*r*g*se**3.0))
    return ws
##################################################################
###################critical shear velocity########################
def usstarcr(dg):
    uscr=0.00001
    nu=1.0004E-6
    rhos=2650.0
    rho=1000.0
    g=9.81
    gamma=g*(rhos-rho)
    for i in range(20):
        xr=uscr*dg/nu
        yr=0.15*xr**(-1.0)+0.05*exp(-8.0*xr**(-0.9))
        uscr=sqrt(dg*gamma*yr/rho)
    return uscr
###################################################################


def tanh(z):
    return (exp(z) - exp(-z))/(exp(z) + exp(-z))

def coth(z):
    return (exp(z) + exp(-z))/(exp(z) - exp(-z))


def cosh(z):
    return 0.5 * (exp(z) + exp(-z))

# def calc_u_b(u1,u1_c,F1,k1,depth1,delta1,m1,n1):
#     term1 = (1. - F1**2.0 * k1 * (depth1) *
#             coth(k1*(depth1)))/(coth(k1*(depth1)) - F1**2. * k1 *depth1)
#     #print term1, F1**2.0 * k1 * depth1 * tanh(k1*depth1)
#     vel_exec = u1/(u1 - u1_c)
#     return calc_t_bar(u1,u1_c,m1,n1) * vel_exec*n1 * k1 * term1 * cos(k1*delta1)

def calc_u_b(u1,u1_c,F1,k1,depth1,delta1,m1,n1):
    term1 = (1. - F1**2.0 * k1 * (depth1) * tanh(k1 *depth1))/(tanh(k1 * depth1)  - F1**2. * k1 * depth1)
    #print term1, F1**2.0 * k1 * depth1 * tanh(k1*depth1)
    vel_exec = u1/(u1 - u1_c)
    return calc_t_bar(u1,u1_c,m1,n1) * vel_exec*n1 * k1 * term1 * cos(k1*delta1)


def calc_ub_old(u1,u1_c,n1,m1,depth1,totald1,delta1,k1):
    #depth1=abs(depth1)
    return -(n1*calc_t_bar(u1,u1_c,n1,m1)*k1/(9.81*1440.))*coth(totald1-depth1)*cos(k1*delta1)


def calc_t_bar(u1,u1_c,m1,n1):
    return m1 * ( u1 - u1_c)**n1

def calc_a(t1,u1,u1_c,F1,k1,depth1,delta1,m1,n1,A1_o):
    vel_exec = u1/(u1 - u1_c)
    term1 = A1_o /(F1**2.0 * k1 * depth1)
    term2 = cosh(k1*depth1) * (F1**2. * k1 * depth1 - tanh(k1*depth1))
    term3 = (1. - F1**2.0 * k1 * depth1 * tanh(k1*depth1))/(F1**2.0 * k1 * depth1 -tanh(k1*depth1))
    term4 = exp(t1*n1*calc_t_bar(u1,u1_c,m1,n1)*k1**2.0*vel_exec*term3*sin(k1*delta1))
    return term1 * term2 * term4


def message():
    print " ******************"
    print " * Simple bedform *"
    print " *     code       *"
    print " * Kennedy's work *"
    print " ******************"

#### Main Program ######
#F=linspace(0.5,2.0,100)

F=1.1
nt=100
nx=500
A_o=0.01
#vel = 0.5
dg = 1.0 * 1.0e-3
w_s = ferg(dg)
vel_c = usstarcr(dg)*10.0
print "vel", vel_c
length = 5.1809
depth = 2.2
vel=F*sqrt(9.81*depth)
#vel[:]=F[:]*sqrt(10.81*depth)
delta = 0.91
delta = delta *length#/(2.0*pi)
delta = 0.77 * vel_c**2.0 * depth/(w_s**2.0)
delta = delta *length/(2.0*pi)
print "Delta", delta
# totaldepth=10.0
#delta = delta *2.0*pi
m = 0.05
n = 0.64
k =  2. * pi / length
# F = vel / sqrt(depth * 9.81)
# print F,delta/depth
#uu=zeros([100])
#print delta/depth*depth*k/pi
#for  ii in range(100):
u_b = calc_u_b(vel,vel_c,F,k,depth,delta,m,n)
print u_b
#    uu[ii]=u_b
#plot(F,uu)
#print min(u_b),u_b[-1]
#show()
# u_b_old = calc_ub_old(vel,vel_c,n,m,depth, totaldepth,delta,k)
# print u_b,u_b_old
t = linspace(0,20,nt)
a =  calc_a(t,vel,vel_c,F,k,depth,delta,m,n,A_o)
x = linspace(0,4.0*length,nx)
eta = zeros([nt,nx])
xi = zeros([nt,nx])
maxeta = []
maxxi = []
for tt in range(nt):
    print a[tt],u_b,F,k*delta,delta/depth,F**2.0*k*depth*tanh(k*depth)
    for ii in range(nx):
        eta[tt,ii] = -depth + a[tt] * sin(k*(x[ii]-u_b*t[tt]))
        xi[tt,ii] =  sin(k*(x[ii]-u_b*t[tt]))*a[tt] / ((1.0 - tanh(k * depth)/(F**2.0*k*depth))*cosh(k *depth))
    maxeta.append(max(eta[tt,:]))
    maxxi.append(max(xi[tt,:]))

print shape(eta)
# plot(x,eta[50,:]/max(eta[50,:]))
# plot(x,eta[0,:]/max(eta[0,:]),'k-')
figure(1)
plot(x,eta[0,:]/max(eta[0,:]),'k--')
# plot(x,eta[1,:]/max(eta[1,:]),'k-')
plot(x,xi[0,:]/max(xi[0,:]),'b--')
# plot(x,xi[1,:]/max(xi[1,:]),'b-')
# figure(2)
# plot(maxeta)
show()
