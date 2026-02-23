import Chemfiles
using CellListMap
using FastPow
using StaticArrays
using Printf
using Base.Threads
using Parameters
using Statistics: mean
using LinearAlgebra: norm_sqr
using LoopVectorization

#
# Simulation setup
#
@with_kw struct Params{V,N,T,M,UnitCellType}
    x0::V = getcoor("./ne10k_initial.pdb")  
    temperature::T = 300.
    nsteps::Int = 10_000
    dt::T = 1.0 # fs
    ibath::Int = 10
    print_energy::Int = 50 
    print_traj::Int = 100
    trajfile::String = "ne10k_traj.xyz"
    cutoff::T = 12.
    box::Box{UnitCellType,N,T,M} = Box([ 46.37, 46.37, 46.37 ], cutoff)
    # Parameters for Neon
    mass::T = 20.17900 # g/mol 
    ε::T = 0.0441795 # kcal/mol
    σ::T = 2*1.64009 # Å
    kB::T = 0.001985875 # Boltzmann constant kcal / mol K
end

@inline function potential_energy(d2,ffpair,u)
    @fastpow u += ( ffpair.p12/d2^6 - 2*ffpair.p6/d2^3 )
    return u
end

@inline function forces(x,y,i,j,d2,ffpair,f)
    r = y - x
    @fastpow dudr = -12*(ffpair.p12/d2^7 - ffpair.p6/d2^4)*r
    @inbounds f[i] += dudr
    @inbounds f[j] -= dudr
    return f
end

struct FFPair{T}
    p12::T
    p6::T
end

# Kinetic energy and temperature 
function compute_kinetic(v::AbstractVector{SVector{N,T}},m::T) where {N,T}
    kinetic = zero(T)
    @tturbo for i in eachindex(v)
        kinetic += norm_sqr(v[i]) 
    end
    kinetic = (m/2) * kinetic
    return kinetic
end
compute_temp(kinetic,kB,n) = 2*kinetic/(3*kB*n)
compute_temp(v::AbstractVector,m,kB) = 2*compute_kinetic(v,m)/(3*kB*length(v))

# Remove drift from velocities
function remove_drift!(v::AbstractVector{SVector{N,T}}) where {N,T}
    vmean = zero(SVector{N,T}) 
    @tturbo for i in eachindex(v)
        vmean += v[i]
    end
    vmean = vmean / length(v)
    @tturbo for i in eachindex(v)
        v[i] = v[i] - vmean
    end
end

# Function to print output data
function print_data(istep,x,params,ffpair,cl,kinetic,trajfile)
    @unpack print_energy, print_traj, kB, box, ε, σ = params
    if istep%print_energy == 0
        u = map_pairwise!( 
            (x,y,i,j,d2,output) -> potential_energy(d2,ffpair,output),
            0., box, cl, parallel=true,
        ) 
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
    Chemfiles.close(traj)
    return copy(reinterpret(reshape,SVector{3,Float64},Chemfiles.positions(frame)))
end

#
# Simulation
#
function simulate_fast(params::Params{V,N,T,UnitCellType}) where {V,N,T,UnitCellType}
    @unpack x0, temperature, nsteps, box, dt, ε, σ, mass, kB = params
    trajfile = open(params.trajfile,"w")

    # To use coordinates in Angstroms, dt must be in 10ps. Usually packages
    # use ps and nm internally (thus multiply coordinates by 10 and divide
    # the timestep given in fs by 1000)
    dt = dt/100

    ffpair = FFPair(ε*σ^12,ε*σ^6)

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

    # Build cell lists for the first time
    cl = CellList(x,box)

    # preallocate threaded output, since it contains the forces vector
    f .= Ref(zeros(eltype(f)))
    f_threaded = [ deepcopy(f) for _ in 1:nthreads() ]
    aux = CellListMap.AuxThreaded(cl)

    # Print data at initial point
    kinetic = compute_kinetic(v,mass)
    print_data(0,x,params,ffpair,cl,kinetic,trajfile)

    # Simulate
    for istep in 1:nsteps

        # Update positions (velocity-verlet)
        @tturbo for i in eachindex(x)
            x[i] = x[i] + v[i]*dt + 0.5*(f[i]/mass)*dt^2
        end

        # Reset forces
        @tturbo for i in eachindex(f)
            flast[i] = f[i]
        end
        fill!(f,zeros(eltype(f)))
        @threads for it in 1:nthreads()
            fill!(f_threaded[it],zeros(eltype(f)))
        end

        # Update forces
        map_pairwise!( 
            (x,y,i,j,d2,output) -> forces(x,y,i,j,d2,ffpair,output),
            f, box, cl, parallel=true,
            output_threaded=f_threaded
        ) 
         
        # Update velocities
        @tturbo for i in eachindex(v)
            v[i] = v[i] + 0.5*((flast[i] + f[i])/mass)*dt 
        end

        # Print data and output file
        kinetic = compute_kinetic(v,mass)
        print_data(istep,x,params,ffpair,cl,kinetic,trajfile)

        # Isokinetic bath
        if istep%params.ibath == 0
            remove_drift!(v)
            temp = compute_temp(kinetic,kB,length(v))
            @tturbo for i in eachindex(v)
                v[i] = v[i] * sqrt(temperature/temp)
            end
        end

        # Update cell lists
        cl = UpdateCellList!(x,box,cl,aux)

   end
    close(trajfile)

end

#params = Params()
#simulate(params)





