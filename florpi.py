import numpy as np 
import time
import sys
from halotools.mock_observables import mean_radial_velocity_vs_r

ns = [  10000  ,
        50000  ,
        100000 ,
        200000 ,
        300000 ,
        400000 ,
        500000 ,
        600000 ,
        700000 ,
        800000 ,
        900000 ,
        1000000,
        1500000,
        2000000,
        3000000,
        4000000,
        5000000,
        6000000,
        7000000,
        8000000,
        9000000,
        10000000 ]


rbins = [0.,2.,4.,6.,8.,10.]

for type in [ "cv", "cd" ] :

    if type == "cd" :
        density = 10**5 / 250**3
        boxsize = (n_halos / density)**(1/3)
        _ns = ns
    else : 
        _ns = ns[0:15]
        boxsize = 250.
    
    Lbox = [boxsize,boxsize, boxsize]

    for nthreads in [ 8, 1 ] :
    
        for n_halos in _ns :
        
            if n_halos == 100_000 :
                positions = np.loadtxt("./data/florpi100k_positions.dat")
                velocities = np.loadtxt("./data/florpi100k_velocities.dat")
            else :
                positions = boxsize*np.random.random((n_halos,3))
                velocities = np.random.random((n_halos,3))
            
            t0 = time.time()
            v_12_halotools = mean_radial_velocity_vs_r(
                positions, 
                velocities, 
                rbins_absolute=rbins,
                period=Lbox,
                num_threads=nthreads,
            )
            t_halotools = time.time()
            print(type, nthreads, n_halos, t_halotools-t0)
        
        
