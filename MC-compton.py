import random 
import math 
import numpy as np
import matplotlib.pyplot as plt

# Lab Frame

# Constants
pi = math.pi
m = 0.000510 # mass of electron in gev  
alpha = 1/132.507 # alpha
E = .1 # initial photon energy in MeV
E = E*1E-3 # initial photon energy in Gev 
delta = 2.0
w = E #(E - m**2)/(2*m)
w = w/m
convert = 3.894E8 # GeV^-2 to pb

def dsigma(costh):
    pre_fact = pi*(alpha/m)**2
    pre_fact *= 1/(1 +  w*(1-costh))**2
    sum_1 = 1/(1 + w*(1-costh))
    sum_1 += w*(1-costh)
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
    plt.hist(var*1E6,bins=40)
    string = "Photon " + var_name + " keV"
    plt.xlabel(string)
    plt.ylabel("Counts")
    string_new = string.replace(" ","_")
    plt.savefig(string_new+".png")


def plot_elec(var,var_name):
    plt.figure()
    plt.hist(var*1E6,bins=40)
    string = "Electron " + var_name + " keV"
    plt.xlabel(string)
    plt.ylabel("Counts")
    string_new = string.replace(" ","_")
    plt.savefig(string_new+".png")

def plot_2d(var_1,var_2,var_1_name,var_2_name):
    plt.figure()
    plt.scatter(var_1,var_2*1E6)
    plt.xlabel(var_1_name)
    plt.ylabel(var_2_name)
    plt.savefig(var_1_name+"_"+var_2_name.replace(" ","_")+".png")

def gen_events(Nevent,w_max):
    i = 0

    #sum_ = np.zeros(Nevent)
    m_costh = np.zeros(Nevent)
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
            m_costh[i] = math.degrees(math.acos(costh))
            sinth = math.sqrt(1 - costh**2)
            ph_e[i]  = m*E/(E*(1-costh) + m)   #E/(1+((E/m)*(1-costh))) 
            ph_px[i] = 0
            ph_py[i] = ph_e[i] * sinth 
            ph_pz[i] = ph_e[i] * costh 


            el_px[i] = - ph_px[i] 
            el_py[i] = - ph_py[i] 
            el_pz[i] = E - (ph_e[i] * costh) #math.sqrt(el_e[i]**2 - m**2  - el_px[i]**2)
            el_e[i]  = E - ph_e[i] # math.sqrt(m**2 + el_px[i]**2 + el_py[i]**2 + el_pz[i]**2)
            i = i + 1

    plot_photon(ph_px,"Px")
    plot_photon(ph_py,"Py")
    plot_photon(ph_pz,"Pz")
    plot_photon(ph_e,"E")


    plot_elec(el_px,"Px")
    plot_elec(el_py,"Py")
    plot_elec(el_pz,"Pz")
    plot_elec(el_e,"E")

    plot_2d(m_costh,ph_e,"Angle","Photon Energy")
    plot_2d(m_costh,el_e,"Angle","Electron Energy")
    #plot_photon(sum_,"E_sum")


def main():
    N = 100000
    Nevent = 10000
    w_sum = 0
    w_max = 0
    w_square = 0
    costh = np.zeros(N)
    for i in range(1,N):
       w_i, angle, w_max = xsec(w_max)
       w_sum += w_i
       w_square += w_i*w_i
       costh[i] = angle
    variance = math.sqrt(w_square/N - (w_sum/N)**2)
    gen_events(Nevent,w_max)
    sigma_x = w_sum*convert/(N*1E12)
    print("Number of Events ", Nevent)
    print("Energy  ", E*1000, "MeV")
    print("Effective Luminsoity ", Nevent/sigma_x, "barn^-1")
    print("cross section in barn ", sigma_x, " barn")
    print("error in cross section in barn ", variance*convert/(math.sqrt(N)*1E12), " barn")


main()


'''
pre_fact = pi*alpha*alpha/(m**2)
pre_fact *= (1+(E*(1-costh)/m))**2
sum_1 = (1+(E*(1-costh)/m))
sum_1 += E * (1-costh)/m
sum_1 += costh**2
'''