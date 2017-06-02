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

def random_floats(low, high, size):
    # returns an array of random numbers of size 'size' between 'low' and 'hihg'.
    return [random.uniform(low, high) for _ in xrange(size)]


def message():
    print " ******************"
    print " * Simple bedform *"
    print " *     code       *"
    print " * Kennedy's work *"
    print " ******************"

#### Main Program ######
#F=linspace(0.5,2.0,100)

#F=0.8
mn_n = 3000
Froude = random_floats(0.5,2.0,mn_n)
kd = random_floats(4.0,4.0,mn_n)
j = random_floats(0.8,0.8,mn_n)
dd = random_floats(2.0,2.0,mn_n)
Froude = linspace(0.5,2.5,mn_n)
kd = linspace(0.5,3.0,mn_n)
nt=100
nx=500
A_o=0.6
#kd = 0.1
#vel = 0.5
dg = 1.0 * 1.0e-3
w_s = ferg(dg)
vel_c = usstarcr(dg)*10.0
print "vel", vel_c
#j = 2.0
depth = 1.0

# vel=F*sqrt(9.81*depth)

#vel[:]=F[:]*sqrt(10.81*depth)
#delta =0.8
# delta = delta *length#/(2.0*pi)
# delta = 0.77 * vel_c**2.0 * depth/(w_s**2.0)
# delta = j* depth * 2.* pi / k
# print "Delta", delta

m = 0.5
n = 2.64
sad_x = []
sad_y = []
sad_z = []
sad_v = []
snd_x = []
snd_y = []
snd_z = []
snd_v = []

uad_x = []
uad_y = []
uad_z = []
uad_v = []
und_x = []
und_y = []
und_z = []
und_v = []

for kk1 in range(mn_n):
    for kk2 in range(mn_n):
        F=Froude[kk1]
        depth=dd[kk1]
        vel=F*sqrt(9.81*depth)
        length = 2.0*pi*depth/kd[kk2]
        k =  2. * pi / length
        delta = j[kk1]* length
        u_b = calc_u_b(vel,vel_c,F,k,depth,delta,m,n)
        print u_b,kd[kk1],j[kk1],Froude[kk1]

        t = linspace(0,30.0,nt)
        dt = t[1]-t[0]
        a= zeros(nt)
        a[0] = A_o/1.
        for i in range(1,nt):
            term1 = (tanh(k * depth) - (F**2. * k * depth)**(-1.)) / (1. - tanh(k *depth)/(F**2.0 * k * depth))
            a[i]=a[i-1]- dt*(a[i-1] * calc_t_bar(vel,vel_c,m,n) * k**2.0 * vel / (vel - vel_c) *term1 * sin(k * delta))

    # a =  calc_a(t,vel,vel_c,F,k,depth,delta,m,n,A_o)
        # x = linspace(0,4.0*length,nx)
        # eta = zeros([nt,nx])
        # xi = zeros([nt,nx])
        # maxeta = []
        # maxxi = []
        # maxdistance_eta = []
        # maxdistance_xi = []
        # for tt in range(nt):
        # #print a[tt],u_b,F,k*delta,delta/depth,F**2.0*k*depth*tanh(k*depth)
        #     for ii in range(nx):
        #         eta[tt,ii] = -depth + a[tt] * sin(k*(x[ii]-u_b*t[tt]))
        #         xi[tt,ii] = sin(k*(x[ii]-u_b*t[tt]))*a[tt] / ((1.0 - tanh(k * depth)/(F**2.0*k*depth))*cosh(k *depth))
        # noinspection PyTypeChecker
        if 0.90 < a[10]/a[70] < 1.1:
            if u_b < 0.0:
                sad_x.append(dd[kk1])
                sad_y.append(Froude[kk1])
                sad_z.append(j[kk1])
                sad_v.append(kd[kk2])
            if u_b >= 0.0:
                snd_x.append(dd[kk1])
                snd_y.append(Froude[kk1])
                snd_z.append(j[kk1])
                snd_v.append(kd[kk2])
        # noinspection PyTypeChecker
        else:
            if u_b < 0.0:
                uad_x.append(dd[kk1])
                uad_y.append(Froude[kk1])
                uad_z.append(j[kk1])
                uad_v.append(kd[kk2])
            if u_b >= 0.0:
                und_x.append(dd[kk1])
                und_y.append(Froude[kk1])
                und_z.append(j[kk1])
                und_v.append(kd[kk2])

und_x = array(und_x)
und_y = array(und_y)
und_z = array(und_z)
und_v = array(und_v)

uad_x = array(uad_x)
uad_y = array(uad_y)
uad_z = array(uad_z)
uad_v = array(uad_v)

snd_x = array(snd_x)
snd_y = array(snd_y)
snd_z = array(snd_z)
snd_v = array(snd_v)

sad_x = array(sad_x)
sad_y = array(sad_y)
sad_z = array(sad_z)
sad_v = array(sad_v)

#

ar3=column_stack((snd_x,snd_y,snd_z,snd_v))
ar4=column_stack((sad_x,sad_y,sad_z,sad_v))

ar5=column_stack((und_x,und_y,und_z,und_v))
ar6=column_stack((uad_x,uad_y,uad_z,uad_v))

outfile1=open('snorm_dunes.dat','w')
outfile2=open('santi_dunes.dat','w')
outfile3=open('unorm_dunes.dat','w')
outfile4=open('uanti_dunes.dat','w')
savetxt(outfile1, ar3,fmt='%6.3e')
savetxt(outfile2, ar4,fmt='%6.3e')
savetxt(outfile3, ar5,fmt='%6.3e')
savetxt(outfile4, ar6,fmt='%6.3e')
outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()

print 'Length Stable Antidune', len(sad_v)
print "Length Unstable Antidudes",len(uad_v)

print 'Length Stable Dunes', len(snd_v)
print 'Length Unstable Dunes',len(und_v)
# for tt in range(nt):
#     for ii in range(0,nx-2):
#         if x[ii]>0.25*length and eta[tt,ii+1]< eta[tt,ii]:
#             maxeta.append(eta[tt,ii])
#             break
# for tt in range(nt):
#     for ii in range(0, nx-2):
#         if x[ii]>0.25*length and xi[tt,ii+1] < xi[tt,ii]:
#             maxxi.append(xi[tt,ii])
#             break
# for tt in range(nt):
#     for ii in range(nx):
#         if x[ii]>0.0:
#             if eta[tt,ii] == maxeta[tt]:
#                 maxdistance_eta.append(x[ii])
#                 break
# for tt in range(nt):
#     for ii in range(nx):
#         if x[ii]>0.0:
#             if xi[tt,ii] == maxxi[tt]:
#                 maxdistance_xi.append(x[ii])
#                 break
# print shape(eta)
# # plot(x,eta[50,:]/max(eta[50,:]))
# # plot(x,eta[0,:]/max(eta[0,:]),'k-')
# for tt in range(200):
#     figure(1)
#     inn = tt
#     plot(x,eta[0,:],'r-')
#     plot(x,eta[inn,:],'k--')
# # plot(x,eta[1,:]/max(eta[1,:]),'k-')
#     plot(x,xi[inn,:]/15.,'b--')
#     title(str('t={test1:4.3f}'.format(test1=t[tt])))
#     #vlines(x=maxdistance_eta[inn], ymin=min(eta[inn,:]), ymax=2*max(xi[inn,:]), color='green', zorder=2)
#     #vlines(x=maxdistance_xi[inn], ymin=min(eta[inn,:]), ymax=2*max(xi[inn,:]), color='red', zorder=2)
#     savefig(str('name_{test:04d}.png'.format(test=tt)),dpi=101, bbox_inches="tight")
#     close()
# plot(x,xi[1,:]/max(xi[1,:]),'b-')
# figure(2)
# plot(maxeta)
# plot(ad_x,ad_y,'bo')
# plot(nd_x,nd_y,'ro')
# show()
# print len(ad_x)
