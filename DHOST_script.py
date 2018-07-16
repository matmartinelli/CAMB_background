##matplotlib inline
import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
print('Using CAMB installed at %s'%(os.path.realpath(os.path.join(os.getcwd(),'..'))))
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
import camb
from camb import model, initialpower
print('CAMB version: %s '%camb.__version__)

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.122, mnu=0.06, omk=0, tau=0.06, c2_dhost=3, c3_dhost=5, c4_dhost=1, beta_dhost=-5.3, inired=20.0)
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

print pars
#calculate results for these parameters
results = camb.get_results(pars)


z  = np.linspace(0,4,100)
DA = results.angular_diameter_distance(z)
plt.plot(z, DA)
plt.xlabel('$z$')
plt.ylabel(r'$D_A /\rm{Mpc}$')
plt.title('Angular diameter distance')
plt.ylim([0,2000]);

plt.show()
