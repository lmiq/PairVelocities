import numpy as np
import time
from halotools.mock_observables import mean_radial_velocity_vs_r

boxsize = 250.

#density = 10**4 / 250**3
#boxsize = (n_halos / density)**(1/3)

n_halos = 100_000

n_threads = 1
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
    num_threads=1,
)
t_halotools = time.time()
print("Mean radial velocity = ", v_12_halotools)
print(f"Halotools took {t_halotools-t0} seconds")


