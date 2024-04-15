import numpy as np
import time
import sys
from halotools.mock_observables import mean_radial_velocity_vs_r

nthreads = int(sys.argv[1])

if len(sys.argv) == 2 :
    n_halos = 100_000
else :
    n_halos = int(sys.argv[2])

if len(sys.argv) == 4 :
    density = 10**5 / 250**3
    boxsize = (n_halos / density)**(1/3)
else : 
    boxsize = 250.

rbins = [0.,2.,4.,6.,8.,10.]
Lbox = [boxsize,boxsize, boxsize]

if n_halos == 100_000 :
    positions = np.loadtxt("./data/florpi100k_positions.dat")
    velocities = np.loadtxt("./data/florpi100k_velocities.dat")
else :
    positions = boxsize*np.random.random((n_halos,3))
    velocities = np.random.random((n_halos,3))

print("------------------------------------------------")
print("nthreads = ", nthreads)
print("volume = ", boxsize ** 3)
print("n, dim = ", positions.shape)
t0 = time.time()
v_12_halotools = mean_radial_velocity_vs_r(
    positions, 
    velocities, 
    rbins_absolute=rbins,
    period=Lbox,
    num_threads=nthreads,
)
t_halotools = time.time()
print("Mean radial velocity = ", v_12_halotools)
print(f"Halotools took {t_halotools-t0} seconds")



