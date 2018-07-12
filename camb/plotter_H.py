import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

mod = np.loadtxt('test_H.dat')
std = np.loadtxt('test_H_nomin.dat')

plt.plot(mod[:,0], mod[:,1], ls='-', color='#8E001C' , label=r'minimized $H^{\rm DHOST}(z)$')
plt.plot(mod[:,0], mod[:,2], ls='-', color='#FFB300', label=r'$H^{\rm DHOST}(z)$')

plt.plot(std[:,0], std[:,1], ls='--', color='#8E001C' , label=r'not-minimized $H^{\rm DHOST}(z)$')

plt.legend(loc='upper left', fontsize='small')
plt.xlabel(r'$z$')
plt.ylabel(r'$H(z)$')
plt.xlim([0,2])
plt.ylim([0,200])
plt.savefig('hubble.pdf')
