import pylab
from pylab import *
from numpy import *
import matplotlib
from scipy.signal import correlate
from numpy.fft import fft, ifft, fftshift
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
from scipy import signal, fftpack
#import click


import warnings
warnings.filterwarnings("ignore")
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
    return [random.uniform(low, high) for _ in range(size)]

def stable(tt, Tb , kk, vel_e, fr,dd,de, nn):
    term1 = tt * nn * vel_e
    term2 = (1.0 - fr**2.0*kk*dd*tanh(kk*dd))/(fr**2.0*kk*dd - tanh(kk*dd))
    return term1 * term2 * sin(kk*de)

def cross_correlation_using_fft(x, y):
    f1 = fft(x)

    # flip the signal of y
    f2 = fft(np.flipud(y))
    cc = np.real(ifft(f1 * f2))

    return fftshift(cc)

def compute_shift(x, y):
    # we make sure the length of the two signals are the same
    assert len(x) == len(y)
    c = cross_correlation_using_fft(x, y)
    assert len(c) == len(x)
    zero_index = int(len(x) / 2) - 1
    shift = zero_index - np.argmax(c)
    return shift

def err_func(p):
    return interp1d(X,Y)(X[1:-1]+p[0]) - Y_shifted[1:-1]
def get_max_correlation(original, match):
    z = signal.fftconvolve(original, match[::-1])
    lags = np.arange(z.size) - (match.size - 1)
    return ( lags[np.argmax(np.abs(z))] )


def find_phase(arr1,arr2,locs,le):
    nnx = len(arr1)
    mm1 = max(arr1)
    mm2 = max(arr2)
    # print mm1,mm2
    for i in range(nnx):
        # if locs[i]>locs[10]:
        if arr1[i] == mm1:
            index1 = i
            # print index1, locs[index1]
            break
    for i in range(nnx):
        # if locs[i]>locs[10]:
        if arr2[i] == mm2:
            index2 = i
            # print index2, locs[index2]
            break
    return index1,index2


def message():
    print()
    print("\t **************************************")
    print("\t * \t  \033[1m Creating bed forms \033[0m        *")
    print("\t *          Monte Carlo code          *")
    print("\t *   Based on work from J.F. Kennedy  *")
    print("\t *    Copyright and License to        *")
    print("\t *  CoastalHazards@VT (Robert Weiss)  *")
    print("\t *                                    *")
    print("\t *          \033[1mDeveloper Mode \033[0m           * ")
    print("\t **************************************")
    print()
    print()
    print()


def load_model(mn_n,nprocess):
    import multiprocessing
    from multiprocessing import Pool
    from numpy import array, loadtxt,concatenate
    import tqdm
    global Froude
    global kd
    global j
    global dd
    global vel_c,m,n,nt,A_o,nx
    Froude = random_floats(0.1,5.57,mn_n)
    kd = random_floats(0.1,4.0,mn_n)
    j = random_floats(0.1,0.9,mn_n)
    dd = random_floats(0.5,5.0,mn_n)
    nt=100
    nx=500
    A_o=2.0
    dg = 1.25 * 1.0e-3
    w_s = ferg(dg)
    vel_c = usstarcr(dg)*1.0

    depth = 1.0
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
    result1=[]
    pool = Pool(processes=nprocess)
    inputs = range(mn_n)
    for dummy in tqdm.tqdm(pool.imap_unordered(func1, inputs), total=mn_n):
        result1.append(dummy)
        #print(shape(dummy))
        #result1= concatenate([result1,dummy])
        pass
    print(shape(result1))
#     result1 = pool.map(func1, inputs)
    pool.close()
    pool.terminate()
    array_n = len(result1)
    print(len(result1),array_n)

    for array_i in range(array_n):
        if result1[array_i] == 1:
            sad_x.append(dd[array_i])
            sad_y.append(Froude[array_i])
            sad_z.append(array_i)
            sad_v.append(kd[array_i])
        if result1[array_i] == 2:
            snd_x.append(dd[array_i])
            snd_y.append(Froude[array_i])
            snd_z.append(array_i)
            snd_v.append(kd[array_i])
        if result1[array_i] == 3:
            uad_x.append(dd[array_i])
            uad_y.append(Froude[array_i])
            uad_z.append(array_i)
            uad_v.append(kd[array_i])
        if result1[array_i] == 4:
            und_x.append(dd[array_i])
            und_y.append(Froude[array_i])
            und_z.append(array_i)
            und_v.append(kd[array_i])
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
    return result1

def analyze_results(array):
    array_n = len(array)
    sad_value = 0
    uad_value = 0
    snd_value = 0
    und_value = 0
    for array_i in range(array_n):
        if array[array_i] == 1:
            sad_value += 1
        if array[array_i] == 2:
            snd_value += 1
        if array[array_i] == 3:
            uad_value += 1
        if array[array_i] == 4:
            und_value += 1
    print('Length Stable Antidune', sad_value)
    print("Length Unstable Antidudes", uad_value)

    print('Length Stable Dunes', snd_value)
    print('Length Unstable Dunes', und_value)

    print('Total', sad_value+uad_value+snd_value+und_value)
    print(sad_value, snd_value,uad_value,und_value)

def func1(kk):
    res = 0
    F=Froude[kk]
    depth=dd[kk]
    vel=F*sqrt(9.81*depth)
    length = 2.0*pi*depth/kd[kk]
    k =  2. * pi / length
    delta = j[kk]* length
    u_b = calc_u_b(vel,vel_c,F,k,depth,delta,m,n)
    t = linspace(0,30.0,nt)
    dt = t[1]-t[0]
    a= zeros(nt)
    t1 = 1/(F**2.0 *k * depth)*cosh(k*depth) - sinh(k*depth)
    a[0] = A_o*n*calc_t_bar(vel,vel_c,m,n)*k**2.0* vel / (vel - vel_c)*t1*sin(k*delta)
    for i in range(1,nt):
        term1 = (tanh(k * depth) - (F**2. * k * depth)**(-1.)) / (1. - tanh(k *depth)/(F**2.0 * k * depth))
        a[i]=a[i-1]- dt*(a[i-1] * calc_t_bar(vel,vel_c,m,n) * k**2.0 * vel / (vel - vel_c) *term1 * sin(k * delta))
    stab = stable(t[-1], calc_t_bar(vel,vel_c,m,n) , k, vel/(vel - vel_c),
    F, depth,delta, n)
# a =  calc_a(t,vel,vel_c,F,k,depth,delta,m,n,A_o)
    x = linspace(0,4.0*length,nx)
    eta = zeros([nt,nx])
    xi = zeros([nt,nx])
    maxeta = []
    maxxi = []
    maxdistance_eta = []
    maxdistance_xi = []
    max_eta = []
    max_xi = []
    for tt in range(nt):
        eta[tt,:] = a[tt] * sin(k*(x[:]-u_b*t[tt]))
        xi[tt,:] = sin(k*(x[:]-u_b*t[tt]))*a[tt] / ((1.0 - tanh(k * depth)/(F**2.0*k*depth))*cosh(k *depth))
    for tt in range(nt):
        max_eta.append(max(eta[tt,:]))
        max_xi.append(max(xi[tt,:]))
    eta_max_loc0,xi_max_loc0=find_phase(eta[0,:],xi[0,:],x,length)
    eta_max_loc1,xi_max_loc1=find_phase(eta[1,:],xi[1,:],x,length)
    shift1 = abs(x[eta_max_loc0] -x[xi_max_loc0])
    shift2 = abs(x[eta_max_loc1] -x[xi_max_loc1])
    if stab <= 0.0 and abs(u_b) <= abs(vel):
#            if u_b < 0.0 and F**2.0 * k *depth  > tanh(k *depth) and F**2.0 * k *depth *tanh(k *depth)<1.0:# and delta<0.5*length:
        if u_b < 0.0 and shift1 == 0.0 and shift1 == shift2 and max_eta[0] <= max_eta[-1]:# F**2.0 * k *depth *tanh(k *depth)>1.0:# and delta<0.5*length:
            res = 1
        if u_b > 0.0 and shift1>=0.5*length and shift1 == shift2 and max_eta[0] <= max_eta[-1]:#F**2.0 *k *depth < tanh(k * depth):# and delta>=0.5*length and delta <length:
            res = 2
    if stab > 0.0 and abs(u_b) <= abs(vel):
#            if u_b < 0.0 and F**2.0 * k *depth  > tanh(k *depth) and F**2.0 * k *depth *tanh(k *depth)<1.0:# and delta < 0.5*depth:
        if u_b < 0.0 and shift1 == 0.0 and shift1 == shift2 and max_eta[0] <= max_eta[-1]:#F**2.0 * k *depth *tanh(k *depth)<1.0:# and delta < 0.5*depth:
            res = 3
        if u_b >= 0.0 and shift1>=0.5*length and shift1 == shift2 and max_eta[0]<=max_eta[-1]:#F**2.0 *k *depth < tanh(k * depth):#  and delta>=0.5*length and delta <length:
            res = 4
    return res

def default():
    #F=0.8


    # vel=F*sqrt(9.81*depth)

    #vel[:]=F[:]*sqrt(10.81*depth)
    #delta =0.8
    # delta = delta *length#/(2.0*pi)
    # delta = 0.77 * vel_c**2.0 * depth/(w_s**2.0)
    # delta = j* depth * 2.* pi / k
    # print "Delta", delta



    if mode == 'd':
        bar = Bar('Processing', max=mn_n, suffix='%(index)d/%(max)d - %(percent).1f%% - %(eta)ds')

    for kk in range(mn_n):
        F=Froude[kk]
        depth=dd[kk]
        vel=F*sqrt(9.81*depth)
        length = 2.0*pi*depth/kd[kk]
        k =  2. * pi / length
        delta = j[kk]* length
        u_b = calc_u_b(vel,vel_c,F,k,depth,delta,m,n)


        t = linspace(0,30.0,nt)
        dt = t[1]-t[0]
        a= zeros(nt)
        t1 = 1/(F**2.0 *k * depth)*cosh(k*depth) - sinh(k*depth)
        a[0] = A_o*n*calc_t_bar(vel,vel_c,m,n)*k**2.0* vel / (vel - vel_c)*t1*sin(k*delta)

        for i in range(1,nt):
            term1 = (tanh(k * depth) - (F**2. * k * depth)**(-1.)) / (1. - tanh(k *depth)/(F**2.0 * k * depth))
            a[i]=a[i-1]- dt*(a[i-1] * calc_t_bar(vel,vel_c,m,n) * k**2.0 * vel / (vel - vel_c) *term1 * sin(k * delta))
        stab = stable(t[-1], calc_t_bar(vel,vel_c,m,n) , k, vel/(vel - vel_c),
        F, depth,delta, n)
        #print kk,u_b,k*delta,vel/(usstarcr(dg)*F), stab
    # a =  calc_a(t,vel,vel_c,F,k,depth,delta,m,n,A_o)
        x = linspace(0,4.0*length,nx)
        eta = zeros([nt,nx])
        xi = zeros([nt,nx])
        maxeta = []
        maxxi = []
        maxdistance_eta = []
        maxdistance_xi = []
        max_eta = []
        max_xi = []
        for tt in range(nt):
        #print a[tt],u_b,F,k*delta,delta/depth,F**2.0*k*depth*tanh(k*depth)
        #    for ii in range(nx):
            eta[tt,:] = a[tt] * sin(k*(x[:]-u_b*t[tt]))
            xi[tt,:] = sin(k*(x[:]-u_b*t[tt]))*a[tt] / ((1.0 - tanh(k * depth)/(F**2.0*k*depth))*cosh(k *depth))
        for tt in range(nt):
            max_eta.append(max(eta[tt,:]))
            max_xi.append(max(xi[tt,:]))

        #print len(max_eta)
        # noinspection PyTypeChecker
        eta_max_loc0,xi_max_loc0=find_phase(eta[0,:],xi[0,:],x,length)
        eta_max_loc1,xi_max_loc1=find_phase(eta[1,:],xi[1,:],x,length)
        shift1 = abs(x[eta_max_loc0] -x[xi_max_loc0])
        shift2 = abs(x[eta_max_loc1] -x[xi_max_loc1])
        # print 'shift', shift
        # plot(x, eta[0,:])
        # plot(x, xi[0,:])
        # plot(x[eta_max_loc], eta[0,eta_max_loc],'ro')
        # plot(x[xi_max_loc], xi[0,xi_max_loc],'ko')
        #
        # show()
        if stab <= 0.0 and abs(u_b) <= abs(vel):
    #            if u_b < 0.0 and F**2.0 * k *depth  > tanh(k *depth) and F**2.0 * k *depth *tanh(k *depth)<1.0:# and delta<0.5*length:
            if u_b < 0.0 and shift1 == 0.0 and shift1 == shift2 and max_eta[0] <= max_eta[-1]:# F**2.0 * k *depth *tanh(k *depth)>1.0:# and delta<0.5*length:
                sad_x.append(dd[kk])
                sad_y.append(Froude[kk])
                sad_z.append(k)
                sad_v.append(kd[kk])
            if u_b > 0.0 and shift1>=0.5*length and shift1 == shift2 and max_eta[0] <= max_eta[-1]:#F**2.0 *k *depth < tanh(k * depth):# and delta>=0.5*length and delta <length:
                snd_x.append(dd[kk])
                snd_y.append(Froude[kk])
                snd_z.append(k)
                snd_v.append(kd[kk])
        # noinspection PyTypeChecker
        if stab > 0.0 and abs(u_b) <= abs(vel):
    #            if u_b < 0.0 and F**2.0 * k *depth  > tanh(k *depth) and F**2.0 * k *depth *tanh(k *depth)<1.0:# and delta < 0.5*depth:
            if u_b < 0.0 and shift1 == 0.0 and shift1 == shift2 and max_eta[0] <= max_eta[-1]:#F**2.0 * k *depth *tanh(k *depth)<1.0:# and delta < 0.5*depth:
                uad_x.append(dd[kk])
                uad_y.append(Froude[kk])
                uad_z.append(k)
                uad_v.append(kd[kk])
            if u_b >= 0.0 and shift1>=0.5*length and shift1 == shift2 and max_eta[0]<=max_eta[-1]:#F**2.0 *k *depth < tanh(k * depth):#  and delta>=0.5*length and delta <length:
                und_x.append(dd[kk])
                und_y.append(Froude[kk])
                und_z.append(k)
                und_v.append(kd[kk])
        if mode == 'd':
            bar.next()
    if mode == 'd':
        bar.finish()
    # und_x = array(und_x)
    # und_y = array(und_y)
    # und_z = array(und_z)
    # und_v = array(und_v)
    #
    # uad_x = array(uad_x)
    # uad_y = array(uad_y)
    # uad_z = array(uad_z)
    # uad_v = array(uad_v)

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

    if mode == 'd':
        print('Length Stable Antidune', len(sad_v))
        print("Length Unstable Antidudes",len(uad_v))

        print('Length Stable Dunes', len(snd_v))
        print('Length Unstable Dunes',len(und_v))

        print('Total', len(sad_v)+len(uad_v)+len(snd_v)+len(und_v))
