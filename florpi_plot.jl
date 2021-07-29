using Plots, Plots.Measures
using DelimitedFiles
using LaTeXStrings

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


