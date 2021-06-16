import numpy as np
import time
import sys
from halotools.mock_observables import mean_radial_velocity_vs_r

if len(sys.argv) == 1 :
    n_halos = 100_000
else :
    n_halos = int(sys.argv[1])

if len(sys.argv) == 3 :
    density = 10**5 / 250**3
    boxsize = (n_halos / density)**(1/3)
else : 
    boxsize = 250.

n_threads = 8
rbins = [0.,2.,4.,6.,8.,10.]
Lbox = [boxsize,boxsize, boxsize]

if n_halos == 100_000 :
    positions = np.loadtxt("./florpi100k_positions.dat")
    velocities = np.loadtxt("./florpi100k_velocities.dat")
else :
    positions = boxsize*np.random.random((n_halos,3))
    velocities = np.random.random((n_halos,3))

print(positions.shape)
t0 = time.time()
v_12_halotools = mean_radial_velocity_vs_r(
    positions, 
    velocities, 
    rbins_absolute=rbins,
    period=Lbox,
    num_threads=8,
)
t_halotools = time.time()
print("Mean radial velocity = ", v_12_halotools)
print(f"Halotools took {t_halotools-t0} seconds")


