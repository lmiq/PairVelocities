import Chemfiles
using CellListMap
using FastPow
using StaticArrays
using Printf
using Base.Threads
using Parameters
using Statistics: mean
using LinearAlgebra: norm_sqr

# Simulation parameters
@with_kw struct Params{V,N,T,M,UnitCellType}
    x0::V = getcoor("o6.dcd")  
    temperature::T = 300.
    nsteps::Int = 10_000
    dt::T = 0.002
    ibath::Int = 10
    print_energy::Int = 50 
    save_traj::Int = 100
    trajfile::String = "o1.xyz"
    cutoff::T = 10.
    box::Box{UnitCellType,N,T,M} = Box([ 80.      0.      0.
                                          0.     80.      0. 
                                          0.      0.     80. ], cutoff)
end

# function that computes energy and forces
@fastpow function energy_and_force(x,y,i,j,d2,ε,σ,uf::UF{T,V}) where {T,V}
    @unpack u, f = uf
    σ6 = σ^6
    σ12 = σ6^2
    d6 = d2^3
    d12 = d6^2
    σ12d12 = σ12/d12
    σ6d6 = σ6/d6 
    u += ε*( σ12d12 - 2*σ6d6 )
    r = y - x
    ∂u∂r = -12*ε*(σ12d12 - σ6d6)*(r/d2)
    f[i] = f[i] + ∂u∂r 
    f[j] = f[j] - ∂u∂r 
    return UF{T,V}(u,f)
end

# Structure that carries both the energy and force vector
@with_kw struct UF{T,V}
    u::T
    f::V
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
compute_temp(kinetic,k,n) = (2/(3k))*kinetic/n
compute_temp(v::AbstractVector,m,k) = (2/(3k))*compute_kinetic(v,m)/length(v)

# Remove drift from velocities
function remove_drift!(v)
    vmean = mean(v)
    v .= v .- Ref(vmean)
end

# Read coordinates from NAMD-DCD file
function getcoor(file)
    traj = Chemfiles.Trajectory(file)
    frame = Chemfiles.read_step(traj,0)
    return copy(reinterpret(reshape,SVector{3,Float64},Chemfiles.positions(frame)))
end

#
# Simulation
#
function simulate(params::Params{V,N,T,UnitCellType}) where {V,N,T,UnitCellType}
    @unpack x0, temperature, nsteps, box, dt = params
    trajfile = open(params.trajfile,"w")

    # Parameters for Neon
    ε = 0.0441795
    σ = 2*1.64009
    mass_Ne = 20.17900 # kg/mol

    x = copy(x0)
    f = similar(x)
    flast = similar(x)

    # Initial velocities
    v = randn(eltype(x),size(x))
    # Adjust average to desidred temperature
    k = 0.001985875 # Boltzmann constant kcal / mol K
    t0 = compute_temp(v,mass_Ne,k) 
    @. v = v * sqrt(temperature/t0)
    # Remove drift
    remove_drift!(v)

    # preallocate threaded output, since it contains the forces vector
    uf_threaded = [ UF{T,V}(0.,deepcopy(f)) for _ in 1:nthreads() ]

    # Build cell lists for the first time
    cl = CellList(x,box)

    println("Initial temperature = ", compute_temp(v,mass_Ne,k))

    f .= Ref(zeros(eltype(f)))
    uf = UF{T,V}(0.,f)
    uf = map_pairwise!( 
        (x,y,i,j,d2,output) -> energy_and_force(x,y,i,j,d2,ε,σ,output),
        uf, box, cl, parallel=true,
        reduce=reduceuf,
        output_threaded=uf_threaded
    )
    println("Energy at initial point = ", uf.u)

    for istep in 1:nsteps

        # Reset energy and forces
        flast .= f
        u = 0.
        f .= Ref(zeros(eltype(f)))
        uf = UF{T,V}(u,f)
        @threads for it in 1:nthreads()
           ft = uf_threaded[it].f 
           ft .= Ref(zeros(eltype(f)))
           uf_threaded[it] = UF{T,V}(0.,ft)
        end

        # Compute energy and forces
        uf = map_pairwise!( 
            (x,y,i,j,d2,output) -> energy_and_force(x,y,i,j,d2,ε,σ,output),
            uf, box, cl, parallel=true,
            reduce=reduceuf,
            output_threaded=uf_threaded
        ) 
        u = uf.u
        kinetic = compute_kinetic(v,mass_Ne)

        if istep%params.print_energy == 0
            temp = compute_temp(kinetic,k,length(v))
            @printf(
                "STEP = %8i U = %12.5f K = %12.5f TOT = %12.5f TEMP = %12.5f\n", 
                istep, u, kinetic, u+kinetic, temp
            )
        end

        if u/length(x) > 1e10
            println("Simulation is unstable. ")
            return
        end

        if istep%params.save_traj == 0
            println(trajfile,length(x))
            println(trajfile," step = ", istep)
            for i in 1:length(x)
               @printf(trajfile,"Ne %12.5f %12.5f %12.5f\n", ntuple(j -> x[i][j], N)...)
            end
        end

        # Update positions and velocities (velocity-verlet)
        @. x = x + v*dt + 0.5*(f/mass_Ne)*dt^2
        @. v = v + 0.5*(f + flast)*dt 

        # Isokinetic bath
        if istep%params.ibath == 0
            temp = compute_temp(kinetic,k,length(v))
            @. v = v * (temperature/temp)
            remove_drift!(v)
        end

        # Update cell lists
        cl = UpdateCellList!(x,box,cl)

   end
    close(trajfile)

end

#params = Params()
#simulate(params)


















