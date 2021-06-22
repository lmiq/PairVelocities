import numpy as np
import time
import sys
from halotools.mock_observables import mean_radial_velocity_vs_r

n_threads = 8
rbins = list(np.arange(0.,201.,20.))
print(rbins)
boxsize = 2000.
Lbox = [boxsize,boxsize, boxsize]

N = 8_000_000

print("Reading positions...")
positions = np.loadtxt("~/Downloads/data_enrique/pos.dat",max_rows=N)
print("Reading velocities...")
velocities = np.loadtxt("~/Downloads/data_enrique/vel.dat",max_rows=N)

print("Computing histogram...")
t0 = time.time()
v_12_halotools = mean_radial_velocity_vs_r(
    positions, 
    velocities, 
    rbins_absolute=rbins,
    period=Lbox,
    num_threads=n_threads,
)
t_halotools = time.time()
for v in v_12_halotools :
  print(v)
#print("Mean radial velocity = ", v_12_halotools)

print(f"Halotools took {t_halotools-t0} seconds")


