import random 
import math 
import numpy as np
import matplotlib.pyplot as plt

# Constants
pi = math.pi
m = 0.000510 # gev 
alpha = 1/132.507 # alpha
w = 0.000001 # gev
w = w*1E-6
s = 2*m*w + m**2
delta = 2.0
convert = 3.894E8

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
    print(w_max)
    print("cross section in barn ", w_sum*convert/(N*1E12), " barn")
    print("error in cross section in barn ", variance*convert/(math.sqrt(N)*1E12), " barn")
    plt.figure()
    plt.hist(costh)
    plt.show()

main()