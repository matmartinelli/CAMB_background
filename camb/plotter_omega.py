import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

mod = np.loadtxt('test_modomega.dat')
std = np.loadtxt('test_standard.dat')

plt.plot(mod[:,0], mod[:,1], ls='-', color='#8E001C' , label=r'$\Omega^{\rm DHOST}_m(z)$')
plt.plot(mod[:,0],mod[:,2], ls='-', color='#FFB300', label=r'$\Omega^{\rm DHOST}_\Lambda(z)$')

plt.plot(std[:,0], std[:,1], ls='--', color='#8E001C' , label=r'$\Omega^{\rm std}_m(z)$')
plt.plot(std[:,0], std[:,2], ls='--', color='#FFB300', label=r'$\Omega^{\rm std}_\Lambda(z)$')

plt.legend(loc='upper right', fontsize='small')
plt.xlabel(r'$z$')
plt.ylabel(r'$\Omega(z)$')
plt.xlim([0,2])
plt.savefig('omegas.pdf')
