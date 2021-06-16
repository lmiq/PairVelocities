# PairVelocities

Code to compute pairwise velocity distributions for cosmological simulations.

Current conclusions:

1. Currently `CellListMap.jl` scales linearly with the number of atoms of the system, for constant density. `halotools` scales more than linearly, and becomes slow for *very large* number of particles (more than a million).  

2. The algorithm for pairwise distances of `halotools` runs over less pairs than the `CellListMap.jl`, thus it scales better if the density increases. 


![image](https://user-images.githubusercontent.com/31046348/122289175-c9d6a300-cec8-11eb-86af-dbca257cb3a3.png)

![image](https://user-images.githubusercontent.com/31046348/122289339-fab6d800-cec8-11eb-8b2e-b3d38456fb4e.png)


