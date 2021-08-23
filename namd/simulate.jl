import Chemfiles
using CellListMap
using FastPow
using StaticArrays
using Printf
using Threads
using LinearAlgebra: norm_sq

# Simulation parameters
@with_kw struct Params{N,T,UnitCellType}
    x0::Vector{SVector{N,T}} = getcoor("o1.dcd")  
    temperature::T = 300.
    nsteps::Int = 10_000
    print_energy::Int = 100
    save_traj::Int = 1000
    trajfile::String = "o1.xyz"
    cutoff::T = 10.
    box::Box{UnitCellType,N,T} = Box([ 50.      0.      0.
                                        0.     50.      0. 
                                        0.      0.     50. ], cutoff)
end

# function that computes energy and forces
function uf(i,j,x,y,d2,ε,σ,uf) 
    u, f = uf
    @fastpow sd6 = (σ^2/d2)^3
    sd12 = sd6^2
    u += ε*( sd12 - 2*sd6 )
    dx = y - x
    df = -12*ε*sd6/d2*( sd6 - 1 )*dx
    f[i] = f[i] + df
    f[j] = f[j] - df
    return (u,f)
end

# reduction rule for the (u,f) tuple
function reduceuf(uf,uf_threaded)
    u, f = uf_threaded[1]
    for it in 2:nthreads()
        ut, ft = uf_threaded[it]
        u += ut
        f .+ ft
    end
    return (u,f)
end

# Remove drift from velocities
function remove_drift!(v)
    vmean = mean(v)
    for i in 1:length(v)
      v[i] -= vmean
    end
end

# Read coordinates from NAMD-DCD file
function getcoor(file)
    traj = Chemfiles.Trajectory(file)
    frame = Chemfiles.read_step(traj,0)
    return reinterpret(reshape,SVector{3,Float64},Chemfiles.positions(frame))
end

#
# Simulation
#
function simulate(params::Params{N,T}) where {N,T}
    @unpack x0, temperature, nsteps, box = params
    trajfile = open(params.trajfile,"w")

    # Parameters for Neon
    ε = 0.0441795
    σ = 2*1.64009
    mass_Ne = 20.17900e-3 # kg/mol

    x = copy(x)
    f = similar(x)
    flast = similar(x)

    # Initial velocities
    v = randn(eltype(x),size(x))
    # Remove drift
    # Adjust average to desidred temperature  for i in 1:n
    k = 1.38e-23*1e20 # k in Å
    kT = temperature*k
    @. v = sqrt(3*kT/2)*v
    remove_drift!(v)

    # preallocate threaded output, since it contains the forces vector
    uf_threaded = [ (0.,deepcopy(f)) for it in 1:nthreads() ]

    # Build cell lists for the first time
    cl = CellList(x,box)

    for i in 1:nsteps

        # Reset energy and forces
        @threads for it in 1:nthreads()
           ut, ft = uf_threaded[it]
           ut = zero(T)
           ft .= zero(SVector{3,T}) 
           uf_threaded[it] = (ut,ft)
        end
        # Compute energy and forces
        u, f = map_pairwise( 
            (i,j,x,y,d2,output) -> (i,j,x,y,d2,ε,σ,output),
            box, cl, (u,f),
            reduce=reduceuf,
            output_threaded=uf_threaded
        ) 

        if i%print_energy == 0
            @printf("STEP = %8i ENERGY = %12.5f TEMPERATURE = %12.5f)", i, u, temp)
        end
        if i%save_traj == 0
            println(trajfile,length(x))
            println(trajfile,"Step = ", i)
            for i in 1:length(x)
                if N == 3 
                    @printf(trajfile,"Ne %12.5f %12.5f %12.5f", ntuple(j -> x[i][j], N)...)
                elseif N == 2
                    @printf(trajfile,"Ne %12.5f %12.5f %12.5f", ntuple(j -> x[i][j], N)..., 0.)
                end
            end
        end

        # Update positions and velocities (velocity-verlet)
        @. x = x + v*dt + 0.5*(f/mass_Ne)*dt^2
        @. flast = f
        @. v = v + 0.5*(f + flast)*dt 

        # Isokinetic bath
        if istep <= iequil && mod(istep,ibath) == 0
            temp = (one(T)/2)*mean(el -> norm_sq(el), v)
            @. v = v * temperature/temp 
            remove_drift!(v)
        end

        # Update cell lists
        cl = UpdateCellList!(x,box,cl)

    end
    close(trajfile)

end

params = Params()
simulate(params)


















