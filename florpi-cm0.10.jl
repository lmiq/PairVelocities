import Pkg
Pkg.activate(".")
Pkg.pkg"dev CellListMap"
using CellListMap
using Plots, Plots.Measures
using DelimitedFiles
using LaTeXStrings
using BenchmarkTools
using StaticArrays
using Random
using Chairmarks
using ThreadPinning
import CellListMap: copy_output, reset_output!, reducer!

pinthreads(:cores)

#
# florpi
#

copy_output(x::Tuple{Vector{Int},Vector{Float64}}) = (copy(x[1]), copy(x[2]))
function reset_output!(x::Tuple{Vector{Int},Vector{Float64}}) 
    x[1] .= 0
    x[2] .= 0.0
    return (x[1],x[2])
end
function reducer!(
    x::Tuple{Vector{Int},Vector{Float64}},
    y::Tuple{Vector{Int},Vector{Float64}},
)
    x[1] .+= y[1]
    x[2] .+= y[2]
    return (x[1], x[2])
end

function florpi(;N=100_000,cd=true,parallel=true,nbatches=(0,0))

  @inline dot(x::SVector{3,Float64},y::SVector{3,Float64}) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
  
  function compute_pairwise_mean_cell_lists!(pair,hist,velocities,rbins)
    (; x, y, i, j, d) = pair
    dx = x - y
    ibin = searchsortedfirst(rbins, d) - 1
    hist[1][ibin] += 1
    hist[2][ibin] += dot(velocities[i]-velocities[j],dx)/d
    return hist
  end

  n_halos = N

  if cd
    density = 10^5/250^3  # density of the original problem
    boxsize = (n_halos / density)^(1/3)
  else
    boxsize = 250.
  end
  
  Random.seed!(321)
  Lbox = [boxsize,boxsize,boxsize]
  positions = boxsize .* rand(Float64, 3, n_halos)
  velocities = rand(Float64, 3, n_halos)
  rbins = [0.,2.,4.,6.,8.,10.]
  r_max = maximum(rbins)

  n = size(positions)[2]
  positions = reshape(reinterpret(SVector{3,Float64},positions),n)
  velocities = reshape(reinterpret(SVector{3,Float64},velocities),n)
  
  # Needs this to stabilize the type of velocities and hist, probably
  function barrier(f::F,sys,velocities,rbins) where {F}
    hist = pairwise!(
      (pair, hist) -> f(pair,hist,velocities,rbins),
      sys;
      update_lists=false,
    )
    return hist
  end

  hist = (zeros(Int,length(rbins)-1), zeros(Float64,length(rbins)-1))
  sys = ParticleSystem(
    positions=positions,
    unitcell=Lbox,
    cutoff=r_max,
    output=hist,
    parallel=parallel,
    nbatches=nbatches,
  )

  hist = barrier(compute_pairwise_mean_cell_lists!, sys,velocities, rbins)

  n_pairs = hist[1]
  mean_v_r = hist[2]
  mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0]./n_pairs[n_pairs .> 0]
  return mean_v_r

end

function plot_florpi(version,output=false)

  #
  # Plots
  #
  default(
    fontfamily="Computer Modern",
    label="",linewidth=2,framestyle=:box,legend=:topleft,
    margin=5mm
  )

  data = readdlm("./data/cd_v$version.dat", comments=true, comment_char='#')
  plot(data[:,1],data[:,2],label="Serial/halotools v0.8.2")
  plot!(data[:,1],data[:,3],label="8 cores/halotools v0.8.2")
  plot!(data[:,1],data[:,4],label="Serial/CellListMap.jl $version")
  plot!(data[:,1],data[:,5],label="8 cores/CellListMap.jl $version")
  plot!(xlabel="Number or particles",ylabel="time / s")
  plot!(title=L"\textrm{Constant\ density - \rho=(10^5/250^3) N/V;\ cutoff = 10}")
  if output
    savefig("./data/cd_v$version.png")          
    println("created ./data/cd_v$version.png")
  end

  data = readdlm("./data/cv_v$version.dat", comments=true, comment_char='#')
  plot(data[:,1],data[:,2],label="Serial/halotools v0.8.2")
  plot!(data[:,1],data[:,3],label="8 cores/halotools v0.8.2")
  plot!(data[:,1],data[:,4],label="Serial/CellListMap.jl $version")
  plot!(data[:,1],data[:,5],label="8 cores/CellListMap.jl $version")
  plot!(xlabel="Number or particles",ylabel="time / s")
  plot!(title=L"\textrm{Constant\ volume - V=250^3;\ cutoff = 10}")
  if output 
    savefig("./data/cv_v$version.png")
    println("created ./data/cv_v$version.png") 
  end

end

function run_benchmark(;
    output=false,
    last_cd=10_000_000,
    last_cv=3_000_000,
    types=[true, true, true, true],
    nbatches=(0,0),
)

  if output && (Threads.nthreads() != 8 || !all(types))
    error("To save results, use julia -t 8 and set all types to true.")
  end
  
  ns = [  10000   
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
          10000000 ]

  ilast_cd = findfirst(isequal(last_cd),ns)
  ilast_cv = findfirst(isequal(last_cv),ns)
  
  version=filter(x-> x.second.name == "CellListMap", Pkg.dependencies()) |> x -> first(x)[2].version
  
  println(" Version: v$version ")
  
  #
  # Reading data
  #
  
  previous_cd = sort!(filter(x -> occursin("cd",x) && occursin(".dat",x),readdir("./data")))[end]
  previous_version = previous_cd[5:9] 
  if previous_version == "$version"
    previous_cd = sort!(filter(x -> occursin("cd",x) && occursin(".dat",x),readdir("./data")))[end-1]
  end
  previous_cd = "./data/"*previous_cd
  
  data_cd = readdlm(previous_cd,comments=true,comment_char='#')
  println(" Previous: $previous_cd ")
  
  previous_cv = sort!(filter(x -> occursin("cv",x) && occursin(".dat",x),readdir("./data")))[end]
  previous_version = previous_cv[5:9] 
  if previous_version == "$version"
    previous_cv = sort!(filter(x -> occursin("cv",x) && occursin(".dat",x),readdir("./data")))[end-1]
  end
  previous_cv = "./data/"*previous_cv
  data_cv = readdlm(previous_cv,comments=true,comment_char='#')
  println(" Previous: $previous_cv ")
  
  new_cv = zeros(ilast_cv,5)
  new_cv[1:min(size(data_cv,1),ilast_cv),1:3] .= data_cv[1:min(size(data_cv,1),ilast_cv),1:3]

  new_cd = zeros(ilast_cd,5)
  new_cd[1:min(size(data_cd,1),ilast_cd),1:3] .= data_cd[1:min(size(data_cd,1),ilast_cd),1:3]

  #
  # Parallel 
  #
  if types[1]
      println("Parallel, constant volume:")
      for i in 1:ilast_cv
        GC.gc()
        n = ns[i]
        prev = try data_cv[i,5] catch; 0 end
        t = @b florpi(N=$n,cd=false,parallel=true,nbatches=nbatches);
        new_cv[i,5] = t.time
        println(n," ",t.time," prev ", prev)
      end 
  end
  
  if types[2]
      println("Parallel, constant density:")
      for i in 1:ilast_cd
        GC.gc()
        n = ns[i]
        prev = try data_cd[i,5] catch; 0 end
        t = @b florpi(N=$n,cd=true,parallel=true, nbatches=nbatches);
        new_cd[i,5] = t.time
        println(n," ",t.time," prev ", prev)
      end
  end
      
  #
  # Serial
  #
  if types[3]
      println("Serial, constant density:")
      for i in 1:ilast_cd
        GC.gc()
        n = ns[i]
        prev = try data_cd[i,4] catch; 0 end
        t = @b florpi(N=$n,cd=true,parallel=false,nbatches=nbatches);
        new_cd[i,4] = t.time
        println(n," ",t.time," prev ", prev)
      end
  end
  
  if types[4]
      println("Serial, constant volume:")
      for i in 1:ilast_cv
        GC.gc()
        n = ns[i]
        prev = try data_cv[i,4] catch; 0 end
        t = @b florpi(N=$n,cd=false,parallel=false,nbatches=nbatches);
        new_cv[i,4] = t.time
        println(n," ",t.time," prev ", prev)
      end
  end
  
  if output
    writedlm("./data/cd_v$version.dat",new_cd)
    writedlm("./data/cv_v$version.dat",new_cv)
    plot_florpi("$version",output)
  end

end

@main(args) = run_benchmark()

