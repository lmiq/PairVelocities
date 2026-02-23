import Chemfiles
using CellListMap
using FastPow
using StaticArrays
using Printf
using Base.Threads
using Statistics: mean

#
# Simulation setup
#
@kwdef struct Params{V,N,T}
    x0::V = getcoor("./ne10k_initial.pdb")  
    temperature::T = 300.0
    nsteps::Int = 10_000
    dt::T = 1.0 # fs
    ibath::Int = 10
    print_energy::Int = 50 
    print_traj::Int = 100
    trajfile::String = "ne10k_traj.xyz"
    cutoff::T = 12.0
    unitcell::SVector{N,T} = SVector(46.37, 46.37, 46.37)
    # Parameters for Neon
    mass::T = 20.17900 # g/mol 
    ε::T = 0.0441795 # kcal/mol
    σ::T = 2*1.64009 # Å
    kB::T = 0.001985875 # Boltzmann constant kcal / mol K
end

# Kinetic energy and temperature 
norm_sqr(v::SVector) = sum(abs2, v)
compute_kinetic(v::AbstractVector,m) = (m/2)*sum(norm_sqr, v)
compute_temp(kinetic,kB,n) = 2*kinetic/(3*kB*n)
compute_temp(v::AbstractVector,m,kB) = 2*compute_kinetic(v,m)/(3*kB*length(v))

# Remove drift from velocities
function remove_drift!(v)
    vmean = mean(v)
    v .= v .- Ref(vmean)
end

# Function to print output data
function print_data(istep,x,params,sys,kinetic,trajfile)
    (; print_energy, print_traj, kB, ε, σ) = params
    if istep%print_energy == 0
        u = sys.energy_and_forces.u
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

mutable struct EnergyAndForces{N,T}
    u::T
    f::Vector{SVector{N,T}}
end
CellListMap.copy_output(uf::EnergyAndForces) = EnergyAndForces(uf.u, copy(uf.f))
function CellListMap.reset_output!(uf::EnergyAndForces) 
    uf.u = 0
    fill!(uf.f, zeros(eltype(uf.f)))
    return uf
end
function CellListMap.reducer!(uf1::EnergyAndForces, uf2::EnergyAndForces)
    uf1.u += uf2.u
    uf1.f .+= uf2.f
    return uf1
end

function compute_energy_and_forces(pair,ε,σ,uf::EnergyAndForces)
    (; x, y, i, j, d2) = pair
    r = y - x
    @fastpow u = ε*( σ^12/d2^6 - 2*σ^6/d2^3 )
    @fastpow dudr = -12*ε*(σ^12/d2^7 - σ^6/d2^4)*r
    uf.u += u
    uf.f[i] = uf.f[i] + dudr
    uf.f[j] = uf.f[j] - dudr
    return uf
end

#
# Simulation
#
function simulate(params::Params)
    (; x0, temperature, nsteps, cutoff, unitcell, dt, ε, σ, mass, kB) = params
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

    # Initialize ParticleSystem
    sys = ParticleSystem(
        positions=x,
        unitcell=unitcell,
        cutoff=cutoff,
        output=EnergyAndForces(zero(cutoff), f), 
        output_name=:energy_and_forces,
    )   

    # Compute energy ans forces at initial point
    pairwise!((pair,uf) -> compute_energy_and_forces(pair,ε,σ,uf), sys)

    # Print data at initial point
    kinetic = compute_kinetic(v,mass)
    print_data(0,x,params,sys,kinetic,trajfile)

    # Simulate
    for istep in 1:nsteps
        f = sys.energy_and_forces.f

        # Update positions (velocity-verlet)
        @. x = x + v*dt + 0.5*(f/mass)*dt^2

        # Save forces in previous step
        flast .= f

        # Update forces
        pairwise!((pair,uf) -> compute_energy_and_forces(pair,ε,σ,uf), sys)
         
        # Update velocities
        @. v = v + 0.5*((flast + f)/mass)*dt 

        # Print data and output file
        kinetic = compute_kinetic(v,mass)
        print_data(istep,x,params,sys,kinetic,trajfile)

        # Isokinetic bath
        if istep%params.ibath == 0
            remove_drift!(v)
            temp = compute_temp(kinetic,kB,length(v))
            @. v = v * sqrt(temperature/temp)
        end

        # Update article system
        update!(sys; positions=x)

   end
   close(trajfile)

end

function (@main)(ARGS)
    params = Params()
    @time simulate(params)
end
