from cmath import cos
import random 
import math 
import numpy as np
import matplotlib.pyplot as plt


# Constants
pi = math.pi
m = 0.000510 # gev 
alpha = 1/132.507 # alpha
w = .001 # Mev
w = w*1E-3 # Gev
s = 2*m*w + m**2
delta = 2.0
convert = 3.894E8 # GeV^-2 to pb

def dsigma(costh):
    pre_fact = pi*alpha*alpha/(m**2)
    pre_fact *= (1+(w*(1-costh)/m))**2
    sum_1 = (1+(w*(1-costh)/m))
    sum_1 += w*(1-costh)/m
    sum_1 += costh**2
    return pre_fact*sum_1

def xsec(w_max):
    costh = -1 + random.random()*delta
    w_ii = dsigma(costh)*delta
    if w_max < w_ii:
        w_max = w_ii
    return w_ii, costh, w_max

def plot_photon(var,var_name):
    plt.figure()
    plt.hist(var*1000,bins=100)
    string = "Photon " + var_name + " keV"
    plt.xlabel(string)
    plt.savefig(string+".png")

def plot_elec(var,var_name):
    plt.figure()
    plt.hist(var*1000,bins=100)
    string = "Electron " + var_name + " keV"
    plt.xlabel(string)
    plt.savefig(string+".png")

def gen_events(Nevent,w_max):
    i = 0
    ph_px = np.zeros(Nevent)
    ph_py = np.zeros(Nevent)
    ph_pz = np.zeros(Nevent)
    ph_e = np.zeros(Nevent)

    el_px = np.zeros(Nevent)
    el_py = np.zeros(Nevent)
    el_pz = np.zeros(Nevent)
    el_e = np.zeros(Nevent)

    while i < Nevent:
        costh = -1 + random.random()*delta 
        w_ii =  dsigma(costh) * delta
        prob = w_ii/w_max
        random_num = random.random()
        if random_num<prob: 
            sinth = math.sqrt(1 - costh**2)
            ph_px[i]  = w * sinth 
            ph_py[i] = 0 
            ph_pz[i] = w * costh 
            ph_e[i] = w 

            el_px[i]  = - w * sinth 
            el_py[i] = 0 
            el_pz[i] = w  - w * costh 
            el_e[i] = math.sqrt(2) * w *sinth 
            i = i + 1

    plot_photon(ph_px,"Px")
    plot_photon(ph_py,"Py")
    plot_photon(ph_pz,"Pz")
    plot_photon(ph_e,"E")

    plot_elec(el_px,"Px")
    plot_elec(el_py,"Py")
    plot_elec(el_pz,"Pz")
    plot_elec(el_e,"E")




def main():
    N = 1000000
    w_sum = 0
    w_max = 0
    w_square = 0
    costh = np.zeros(N)
    for i in range(1,N):
       w, angle, w_max = xsec(w_max)
       w_sum += w
       w_square += w*w
       costh[i] = angle
    variance = math.sqrt(w_square/N - (w_sum/N)**2)
    gen_events(100,w_max)
    print("cross section in barn ", w_sum*convert/(N*1E12), " barn")
    print("error in cross section in barn ", variance*convert/(math.sqrt(N)*1E12), " barn")
    plt.figure()
    plt.hist(costh)
    plt.show()

main()