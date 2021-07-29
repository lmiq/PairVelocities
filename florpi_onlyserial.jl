import Pkg
using CellListMap
using Plots, Plots.Measures
using DelimitedFiles
using LaTeXStrings

function run_benchmark()

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
  
  #
  # Serial
  #
  println("Serial, constant density:")
  CellListMap.florpi(N=1000,cd=true,parallel=false);
  times1 = Float64[]
  for i in 1:length(ns)
    n = ns[i]
    prev = data_cd[i,4]
    t = @elapsed CellListMap.florpi(N=n,cd=true,parallel=false);
    push!(times1,t)
    println(n," ",t," prev ", prev)
  end
  
  println("Serial, constant volume:")
  CellListMap.florpi(N=1000,cd=false,parallel=false);
  times3 = Float64[]
  for i in 1:15
    n = ns[i]
    prev = data_cv[i,4]
    t = @elapsed CellListMap.florpi(N=n,cd=false,parallel=false);
    push!(times3,t)
    println(n," ",t," prev ", prev)
  end
  
  return nothing
  
end



