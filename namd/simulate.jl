import Chemfiles
using CellListMap
using FastPow
using StaticArrays
using Printf
using Base.Threads
using Parameters
using Statistics: mean
using LinearAlgebra: norm_sqr

#
# Simulation setup
#
@with_kw struct Params{V,N,T,M,UnitCellType}
    x0::V = getcoor("o6.dcd")  
    temperature::T = 300.
    nsteps::Int = 10_000
    dt::T = 2.0 # fs
    ibath::Int = 10
    print_energy::Int = 50 
    print_traj::Int = 100
    trajfile::String = "traj.xyz"
    cutoff::T = 10.
    box::Box{UnitCellType,N,T,M} = Box([ 80., 80., 80. ], cutoff)
    # Parameters for Neon
    mass::T = 20.17900 # g/mol 
    ε::T = 0.0441795 # kcal/mol
    σ::T = 2*1.64009 # Å
    kB::T = 0.001985875 # Boltzmann constant kcal / mol K
end

# Structure that carries both the energy and force vector
@with_kw struct UF{T,V}
    u::T
    f::V
end

# function that computes energy and forces
@fastpow function energy_and_force(x,y,i,j,d2,ε,σ,uf::UF{T,V}) where {T,V}
    @unpack u, f = uf
    d = sqrt(d2)
    u += ε*( σ^12/d^12 - 2*σ^6/d^6 )
    r = y - x
    dudr = -12*ε*(σ^12/d^13 - σ^6/d^7)*(r/d)
    f[i] = f[i] + dudr
    f[j] = f[j] - dudr
    return UF(u,f)
end

# reduction rule for the (u,f) tuple (u, and f are zeroed outside)
function reduceuf(uf,uf_threaded)
    @unpack u, f = uf
    for it in 1:nthreads()
        u += uf_threaded[it].u
        @. f += uf_threaded[it].f
    end
    return UF(u,f)
end

# Kinetic energy and temperature 
compute_kinetic(v::AbstractVector,m) = (m/2)*sum(x -> norm_sqr(x), v)
compute_temp(kinetic,kB,n) = 2*kinetic/(3*kB*n)
compute_temp(v::AbstractVector,m,kB) = 2*compute_kinetic(v,m)/(3*kB*length(v))

# Remove drift from velocities
function remove_drift!(v)
    vmean = mean(v)
    v .= v .- Ref(vmean)
end

# Function to print output data
function print_data(istep,x,params,u,kinetic,trajfile)
    @unpack print_energy, print_traj, kB = params
    if istep%print_energy == 0
        temp = compute_temp(kinetic,kB,length(x))
        @printf(
            "STEP = %8i U = %12.5f K = %12.5f TOT = %12.5f TEMP = %12.5f\n", 
            istep, u, kinetic, u+kinetic, temp
        )
    end
    if istep%print_traj == 0 && istep > 0
        println(trajfile,length(x))
        println(trajfile," step = ", istep)
        for i in 1:length(x)
           @printf(trajfile,"Ne %12.5f %12.5f %12.5f\n", ntuple(j -> x[i][j], 3)...)
        end
    end
    return nothing
end

# Read coordinates from NAMD-DCD file
function getcoor(file)
    traj = redirect_stdout(() -> Chemfiles.Trajectory(file), devnull)
    frame = Chemfiles.read_step(traj,0)
    return copy(reinterpret(reshape,SVector{3,Float64},Chemfiles.positions(frame)))
end

#
# Simulation
#
function simulate(params::Params{V,N,T,UnitCellType}) where {V,N,T,UnitCellType}
    @unpack x0, temperature, nsteps, box, dt, ε, σ, mass, kB = params
    trajfile = open(params.trajfile,"w")

    # To use coordinates in Angstroms, dt must be in 10ps. Usually packages
    # use ps and nm internally (thus multiply coordinates by 10 and divide
    # the timestep given in fs by 1000)
    dt = dt/100

    # Initial arrays
    x = copy(x0)
    f = similar(x)
    flast = similar(x)

    # Initial velocities
    v = randn(eltype(x),size(x))
    remove_drift!(v)
    # Adjust average to desidred temperature
    t0 = compute_temp(v,mass,kB) 
    @. v = v * sqrt(temperature/t0)
    # Remove drift

    # Build cell lists for the first time
    cl = CellList(x,box)

    # preallocate threaded output, since it contains the forces vector
    f .= Ref(zeros(eltype(f)))
    uf_threaded = [ UF(0.,deepcopy(f)) for _ in 1:nthreads() ]
    aux = CellListMap.AuxThreaded(cl)

    # Compute energy and forces at initial point
    uf = UF(0.,f)
    uf = map_pairwise!( 
        (x,y,i,j,d2,output) -> energy_and_force(x,y,i,j,d2,ε,σ,output),
        uf, box, cl, parallel=true,
        reduce=reduceuf,
        output_threaded=uf_threaded
    ) 
    u = uf.u
    kinetic = compute_kinetic(v,mass)
    print_data(0,x,params,u,kinetic,trajfile)

    # Simulate
    for istep in 1:nsteps

        # Update positions (velocity-verlet)
        @. x = x + v*dt + 0.5*(f/mass)*dt^2

        # Reset energy and forces
        flast .= f
        u = 0.
        f .= Ref(zeros(eltype(f)))
        uf = UF(u,f)
        @threads for it in 1:nthreads()
           ft = uf_threaded[it].f 
           ft .= Ref(zeros(eltype(f)))
           uf_threaded[it] = UF(0.,ft)
        end

        # Compute energy and forces
        uf = map_pairwise!( 
            (x,y,i,j,d2,output) -> energy_and_force(x,y,i,j,d2,ε,σ,output),
            uf, box, cl, parallel=true,
            reduce=reduceuf,
            output_threaded=uf_threaded
        ) 
        u = uf.u
        if u/length(x) > 1e10
            println("Simulation is unstable. ")
            return nothing
        end
         
        # Update velocities
        @. v = v + 0.5*((flast + f)/mass)*dt 

        # Print data and output file
        kinetic = compute_kinetic(v,mass)
        print_data(istep,x,params,u,kinetic,trajfile)

        # Isokinetic bath
        if istep%params.ibath == 0
            temp = compute_temp(kinetic,kB,length(v))
            remove_drift!(v)
            @. v = v * sqrt(temperature/temp)
        end

        # Update cell lists
        cl = UpdateCellList!(x,box,cl,aux)

   end
    close(trajfile)

end

#params = Params()
#simulate(params)





