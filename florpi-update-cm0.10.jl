import Pkg
Pkg.activate(".")
Pkg.pkg"rm CellListMap"
Pkg.pkg"add CellListMap@0.10.3"
Pkg.pkg"pin CellListMap"
Pkg.update()
using CellListMap
using Plots, Plots.Measures
using DelimitedFiles
using LaTeXStrings
using StaticArrays
using Random
using Chairmarks
using ThreadPinning
import CellListMap: copy_output, reset_output!, reducer!

include("./print_line.jl")

pinthreads(:cores)

#
# florpi-update
#
# Companion benchmark to florpi-dev.jl. florpi-dev.jl measures cold cost:
# build a fresh ParticleSystem (full cell list construction) *and* run
# pairwise! once, every sample. That's dominated by construction, and
# construction for a freshly-built system is allocation-heavy (no
# `reset!`-retained scratch buffers to reuse yet) — for large parallel runs,
# most of that wall time can be GC, not compute, which can mask or distort
# genuine scheduling/algorithmic improvements to the parallel map/build
# machinery (see PERFORMANCE_NOTES_pairwise_scaling.md).
#
# This script instead measures *steady-state* cost: build the ParticleSystem
# once per (N, cd, parallel) point, then repeatedly perturb positions by a
# small, MD-step-sized displacement and time `update!` + `pairwise!` only.
# This is the cost that actually recurs on every step of a real simulation
# loop (e.g. Packmol, MD), and it is close to allocation-free since
# `update!`'s cell list rebuild reuses already-sized scratch buffers.
#
# Run as: julia --project -t N florpi-update-dev.jl
# Produces separate data files (cd_update_v*.dat / cv_update_v*.dat) so this
# never collides with florpi-dev.jl's own (cd_v*.dat / cv_v*.dat) baselines —
# the two measure different things and are not comparable to each other.
#
# A single update!+pairwise! step here is fast (single-digit ms up to a few
# hundred ms across the N range), and Chairmarks' `@b` defaults to a 0.1s
# time budget — at these per-call costs that's only 1-10 samples, which is
# not enough to be robust against ordinary system noise (OS scheduling
# jitter, thermal/frequency-scaling transients, etc: repeated runs were
# observed to report spurious 30-100%+ "regressions" at specific N that did
# not reproduce under more robust sampling). BENCHMARK_SECONDS below raises
# the per-point time budget so `@b` collects enough samples for its reported
# "fastest" to be a stable estimate rather than whatever one noisy window
# happened to catch.
#
const BENCHMARK_SECONDS = 2.0

copy_output(x::Tuple{Vector{Int},Vector{Float64}}) = (copy(x[1]), copy(x[2]))
function reset_output!(x::Tuple{Vector{Int},Vector{Float64}})
    x[1] .= 0
    x[2] .= 0.0
    return (x[1], x[2])
end
function reducer!(
    x::Tuple{Vector{Int},Vector{Float64}},
    y::Tuple{Vector{Int},Vector{Float64}},
)
    x[1] .+= y[1]
    x[2] .+= y[2]
    return (x[1], x[2])
end

@inline dotv(x::SVector{3,Float64}, y::SVector{3,Float64}) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

function compute_pairwise_mean_cell_lists!(pair, hist, velocities, rbins)
    (; x, y, i, j, d) = pair
    dx = x - y
    ibin = searchsortedfirst(rbins, d) - 1
    hist[1][ibin] += 1
    hist[2][ibin] += dotv(velocities[i]-velocities[j], dx)/d
    return hist
end

# Needs this to stabilize the type of velocities and hist, probably
function barrier(f::F, sys, velocities, rbins) where {F}
    hist = pairwise!((pair, hist) -> f(pair, hist, velocities, rbins), sys)
    return hist
end

# Builds the ParticleSystem once. Not timed: this is setup, not what this
# benchmark measures (that's florpi-dev.jl's job).
function florpi_update_setup(; N=100_000, cd=true, parallel=true, nbatches=(0, 0))

    n_halos = N

    if cd
        density = 10^5/250^3  # density of the original problem
        boxsize = (n_halos / density)^(1/3)
    else
        boxsize = 250.
    end

    Random.seed!(321)
    Lbox = [boxsize, boxsize, boxsize]
    positions = boxsize .* rand(Float64, 3, n_halos)
    velocities = rand(Float64, 3, n_halos)
    rbins = [0., 2., 4., 6., 8., 10.]
    r_max = maximum(rbins)

    n = size(positions)[2]
    positions = reshape(reinterpret(SVector{3,Float64}, positions), n)
    velocities = reshape(reinterpret(SVector{3,Float64}, velocities), n)

    hist = (zeros(Int, length(rbins)-1), zeros(Float64, length(rbins)-1))
    sys = ParticleSystem(
        positions=positions,
        unitcell=Lbox,
        cutoff=r_max,
        output=hist,
        parallel=parallel,
        nbatches=nbatches,
    )

    return sys, positions, velocities, rbins
end

# One simulation step: perturb positions by a small, MD-step-sized
# displacement (independent of box size — this is not a full reshuffle),
# then update! + pairwise!. Perturbations accumulate across repeated calls
# (a real, if crude, evolving trajectory) rather than resetting each time,
# which is intentional: it exercises `update!` the way a real simulation
# loop would, not a single isolated call.
function florpi_update_step!(sys, positions, velocities, rbins)
    displacement = 0.3
    for i in eachindex(positions)
        positions[i] = positions[i] .+ displacement .* (rand(SVector{3,Float64}) .- 0.5)
    end
    update!(sys; xpositions=positions)
    hist = barrier(compute_pairwise_mean_cell_lists!, sys, velocities, rbins)
    n_pairs = hist[1]
    mean_v_r = hist[2]
    mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0] ./ n_pairs[n_pairs .> 0]
    return mean_v_r
end

function plot_florpi_update(version, output=false)

    default(
        fontfamily="Computer Modern",
        label="", linewidth=2, framestyle=:box, legend=:topleft,
        margin=5mm
    )

    data = readdlm("./data/cd_update_v$version.dat", comments=true, comment_char=('#'))
    plot(data[:, 1], data[:, 2], label="Serial/CellListMap.jl $version")
    plot!(data[:, 1], data[:, 3], label="8 cores/CellListMap.jl $version")
    plot!(xlabel="Number of particles", ylabel="time / s")
    plot!(title=L"\textrm{update! + pairwise! - Constant\ density - \rho=(10^5/250^3) N/V;\ cutoff = 10}")
    if output
        savefig("./data/cd_update_v$version.png")
        println("created ./data/cd_update_v$version.png")
    end

    data = readdlm("./data/cv_update_v$version.dat", comments=true, comment_char=('#'))
    plot(data[:, 1], data[:, 2], label="Serial/CellListMap.jl $version")
    plot!(data[:, 1], data[:, 3], label="8 cores/CellListMap.jl $version")
    plot!(xlabel="Number of particles", ylabel="time / s")
    plot!(title=L"\textrm{update! + pairwise! - Constant\ volume - V=250^3;\ cutoff = 10}")
    if output
        savefig("./data/cv_update_v$version.png")
        println("created ./data/cv_update_v$version.png")
    end

end

function extract_data_version(filename::AbstractString, prefix::AbstractString)
    m = match(Regex("^" * prefix * raw"_v([0-9]+\.[0-9]+\.[0-9]+)\.dat$"), filename)
    return m === nothing ? nothing : VersionNumber(m.captures[1])
end

function pick_previous_data_file(prefix::AbstractString, current_version::VersionNumber; dir="./data")
    files = filter(x -> startswith(x, prefix * "_v") && endswith(x, ".dat"), readdir(dir))
    entries = Tuple{VersionNumber,String}[]
    for file in files
        v = extract_data_version(file, prefix)
        v === nothing && continue
        v == current_version && continue
        push!(entries, (v, file))
    end
    isempty(entries) && return nothing
    sort!(entries, by=first)
    selected = entries[end]
    return joinpath(dir, selected[2]), selected[1]
end

function run_update_benchmark(;
    output=true,
    last_cd=10_000_000,
    last_cv=3_000_000,
    types=[true, true, true, true],
    nbatches=(0, 0),
)

    if output && (Threads.nthreads() != 32 || !all(types))
        error("To save results, use julia -t 32 and set all types to true.")
    end

    ns = [10000
        50000
        100000
        200000
        300000
        400000
        500000
        600000
        700000
        800000
        900000
        1000000
        1500000
        2000000
        3000000
        4000000
        5000000
        6000000
        7000000
        8000000
        9000000
        10000000]

    ilast_cd = findfirst(isequal(last_cd), ns)
    ilast_cv = findfirst(isequal(last_cv), ns)

    version=filter(x -> x.second.name == "CellListMap", Pkg.dependencies()) |> x -> first(x)[2].version

    println(" Version: v$version ")
    println(" Measuring: update! + pairwise! (system built once, positions perturbed each sample) ")

    #
    # Reading previous data, if any (gracefully absent on the first run of
    # this benchmark — there's no historical baseline for this measurement
    # yet, unlike florpi-dev.jl's cd_v*.dat/cv_v*.dat).
    #
    prev_cd = pick_previous_data_file("cd_update", version)
    data_cd, previous_version_cd = if prev_cd === nothing
        println(" No previous cd_update data file found — reporting current times only.")
        nothing, nothing
    else
        path, v = prev_cd
        println(" Previous: $path ")
        readdlm(path, comments=true, comment_char=('#')), v
    end

    prev_cv = pick_previous_data_file("cv_update", version)
    data_cv, previous_version_cv = if prev_cv === nothing
        println(" No previous cv_update data file found — reporting current times only.")
        nothing, nothing
    else
        path, v = prev_cv
        println(" Previous: $path ")
        readdlm(path, comments=true, comment_char=('#')), v
    end

    new_cv = zeros(ilast_cv, 3)
    new_cv[:, 1] .= ns[1:ilast_cv]
    if data_cv !== nothing
        new_cv[1:min(size(data_cv, 1), ilast_cv), 1:2] .= data_cv[1:min(size(data_cv, 1), ilast_cv), 1:2]
    end

    new_cd = zeros(ilast_cd, 3)
    new_cd[:, 1] .= ns[1:ilast_cd]
    if data_cd !== nothing
        new_cd[1:min(size(data_cd, 1), ilast_cd), 1:2] .= data_cd[1:min(size(data_cd, 1), ilast_cd), 1:2]
    end

    #
    # Parallel
    #
    if types[1]
        println("Parallel (update), constant volume:")
        for i in 1:ilast_cv
            n = ns[i]
            sys, positions, velocities, rbins = florpi_update_setup(N=n, cd=false, parallel=true, nbatches=nbatches)
            florpi_update_step!(sys, positions, velocities, rbins) # warmup
            GC.gc()
            prev = data_cv === nothing ? 0 : (try data_cv[i, 3] catch; 0 end)
            t = @b florpi_update_step!($sys, $positions, $velocities, $rbins) seconds=BENCHMARK_SECONDS
            new_cv[i, 3] = t.time
            print_line(n, t.time, prev, string(previous_version_cv))
        end
    end

    if types[2]
        println("Parallel (update), constant density:")
        for i in 1:ilast_cd
            n = ns[i]
            sys, positions, velocities, rbins = florpi_update_setup(N=n, cd=true, parallel=true, nbatches=nbatches)
            florpi_update_step!(sys, positions, velocities, rbins) # warmup
            GC.gc()
            prev = data_cd === nothing ? 0 : (try data_cd[i, 3] catch; 0 end)
            t = @b florpi_update_step!($sys, $positions, $velocities, $rbins) seconds=BENCHMARK_SECONDS
            new_cd[i, 3] = t.time
            print_line(n, t.time, prev, string(previous_version_cd))
        end
    end

    #
    # Serial
    #
    if types[3]
        println("Serial (update), constant density:")
        for i in 1:ilast_cd
            n = ns[i]
            sys, positions, velocities, rbins = florpi_update_setup(N=n, cd=true, parallel=false, nbatches=nbatches)
            florpi_update_step!(sys, positions, velocities, rbins) # warmup
            GC.gc()
            prev = data_cd === nothing ? 0 : (try data_cd[i, 2] catch; 0 end)
            t = @b florpi_update_step!($sys, $positions, $velocities, $rbins) seconds=BENCHMARK_SECONDS
            new_cd[i, 2] = t.time
            print_line(n, t.time, prev, string(previous_version_cd))
        end
    end

    if types[4]
        println("Serial (update), constant volume:")
        for i in 1:ilast_cv
            n = ns[i]
            sys, positions, velocities, rbins = florpi_update_setup(N=n, cd=false, parallel=false, nbatches=nbatches)
            florpi_update_step!(sys, positions, velocities, rbins) # warmup
            GC.gc()
            prev = data_cv === nothing ? 0 : (try data_cv[i, 2] catch; 0 end)
            t = @b florpi_update_step!($sys, $positions, $velocities, $rbins) seconds=BENCHMARK_SECONDS
            new_cv[i, 2] = t.time
            print_line(n, t.time, prev, string(previous_version_cv))
        end
    end

    if output
        writedlm("./data/cd_update_v$version.dat", new_cd)
        writedlm("./data/cv_update_v$version.dat", new_cv)
        plot_florpi_update("$version", output)
    end

end

@main(args) = run_update_benchmark()
