import Pkg
using CellListMap
using Plots, Plots.Measures
using DelimitedFiles
using LaTeXStrings

ENV["GKSwstype"]="nul"

vers(s,c) = parse.(Int,[ split.(split(s,"."),c,keepempty=false)[i][1] for i in 1:3 ])

function prev(c,l::Int=0)
  file = sort!(
    filter(x -> occursin(c,x) && occursin(".dat",x),readdir("./data")),
    by= s->vers(s,"$(c)_v"))[end-l]
  return file[5:end-4] 
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

  data = readdlm("./data/cd_v$version.dat")
  plot(data[:,1],data[:,2],label="Serial/halotools v0.7")
  plot!(data[:,1],data[:,3],label="4 cores/halotools v0.7")
  plot!(data[:,1],data[:,4],label="Serial/CellListMap.jl $version")
  plot!(data[:,1],data[:,5],label="4 cores/CellListMap.jl $version")
  plot!(xlabel="Number or particles",ylabel="time / s")
  plot!(title=L"\textrm{Constant\ density - \rho=(10^5/250^3) N/V;\ cutoff = 10}")
  if output
    savefig("./data/cd_v$version.png")          
    println("created ./data/cd_v$version.png")
  end

  data = readdlm("./data/cv_v$version.dat")
  plot(data[:,1],data[:,2],label="Serial/halotools v0.7")
  plot!(data[:,1],data[:,3],label="4 cores/halotools v0.7")
  plot!(data[:,1],data[:,4],label="Serial/CellListMap.jl $version")
  plot!(data[:,1],data[:,5],label="4 cores/CellListMap.jl $version")
  plot!(xlabel="Number or particles",ylabel="time / s")
  plot!(title=L"\textrm{Constant\ volume - V=250^3;\ cutoff = 10}")
  if output 
    savefig("./data/cv_v$version.png")
    println("created ./data/cv_v$version.png") 
  end

end

function plot_latest() 
  
  latest_version = prev("cd")
  latest_cd = "./data/cd_v"*latest_version*".dat"
  
  data_cd = readdlm(latest_cd,comments=true,comment_char='#')
  println(" Latest: $latest_cd ")
  
  latest_version = prev("cv")
  latest_cv = "./data/cv_v"*latest_version*".dat"
  data_cv = readdlm(latest_cv,comments=true,comment_char='#')
  println(" Latest: $latest_cv ")
  
  plot_florpi("$latest_version",true)

end



