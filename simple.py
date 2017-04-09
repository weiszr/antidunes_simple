import pylab
from pylab import *
from numpy import *
import matplotlib

### Functions #####

def tanh(z):
    return (exp(z) - exp(-z))/(exp(z) + exp(-z))

def coth(z):
    return (exp(z) + exp(-z))/(exp(z) - exp(-z))


def cosh(z):
    return 0.5 * (exp(z) + exp(-z))

def calc_u_b(u1,u1_c,F1,k1,depth1,delta1,m1,n1):
    term1 = (1. - F1**2.0 * k1 * (depth1) *
            coth(k1*(depth1)))/(coth(k1*(depth1)) - F1**2. * k1 *depth1)
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
    term3 = (1. - F1**2.0 * k1 * depth1 * tanh(k1*depth1))/(tanh(k1*depth1) - F1**2. * k1 *depth1)
    term4 = exp(t1*n1*calc_t_bar(u1,u1_c,m1,n1)*k1**2.0*u1/(u1 -u1_c)*term3*sin(k1*delta1))
    return term1 * term2 * term4


def message():
    print " ******************"
    print " * Simple bedform *"
    print " *     code       *"
    print " * Kennedy's work *"
    print " ******************"

#### Main Program ######
#F=linspace(0.5,2.0,100)

F=1.0
nt=100
nx=500
A_o=0.05
#vel = 0.5

vel_c = 0.1
length = 4.0
depth = .25
vel=F*sqrt(9.81*depth)
#vel[:]=F[:]*sqrt(10.81*depth)
delta = 0.0
delta = delta *length#/(2.0*pi)
# totaldepth=10.0
#delta = delta *2.0*pi
m = 1.05
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
for tt in range(nt):
    print a[tt],u_b,F,k*delta,delta/depth,F**2.0*k*depth*tanh(k*depth)
    for ii in range(nx):
        eta[tt,ii] = a[tt] * sin(k*(x[ii]-u_b*t[tt]))
print shape(eta)
# plot(x,eta[50,:]/max(eta[50,:]))
# plot(x,eta[0,:]/max(eta[0,:]),'k-')
plot(x,eta[2,:])
plot(x,eta[0,:],'k-')
show()
