import Pkg
using CellListMap
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

function run_benchmark(output=false,last_cd=10_000_000,last_cv=3_000_000)

  if Threads.nthreads() != 8 
    error(" Run with julia -t auto ")
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
  println("Parallel, constant volume:")
  CellListMap.florpi(N=1000,cd=false,parallel=true);
  for i in 1:ilast_cv
    n = ns[i]
    prev = try data_cv[i,5] catch; 0 end
    t = @elapsed CellListMap.florpi(N=n,cd=false,parallel=true);
    new_cv[i,5] = t
    println(n," ",t," prev ", prev)
  end
  
  println("Parallel, constant density:")
  CellListMap.florpi(N=1000,cd=true,parallel=true);
  for i in 1:ilast_cd
    n = ns[i]
    prev = try data_cd[i,5] catch; 0 end
    t = @elapsed CellListMap.florpi(N=n,cd=true,parallel=true);
    new_cd[i,5] = t
    println(n," ",t," prev ", prev)
  end
  
  #
  # Serial
  #
  println("Serial, constant density:")
  CellListMap.florpi(N=1000,cd=true,parallel=false);
  for i in 1:ilast_cd
    n = ns[i]
    prev = try data_cd[i,4] catch; 0 end
    t = @elapsed CellListMap.florpi(N=n,cd=true,parallel=false);
    new_cd[i,4] = t
    println(n," ",t," prev ", prev)
  end
  
  println("Serial, constant volume:")
  CellListMap.florpi(N=1000,cd=false,parallel=false);
  for i in 1:ilast_cv
    n = ns[i]
    prev = try data_cv[i,4] catch; 0 end
    t = @elapsed CellListMap.florpi(N=n,cd=false,parallel=false);
    new_cv[i,4] = t
    println(n," ",t," prev ", prev)
  end
  
  if output
    writedlm("./data/cd_v$version.dat",new_cd)
    writedlm("./data/cv_v$version.dat",new_cv)
    plot_florpi("$version",output)
  end

end



