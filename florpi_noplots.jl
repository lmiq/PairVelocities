import Pkg
using CellListMap
import CellListMap.Examples: florpi
using DelimitedFiles
using LaTeXStrings
using BenchmarkTools
using Base.Threads

ENV["GKSwstype"]="nul"

vers(s,c) = parse.(Int,[ split.(split(s,"."),c,keepempty=false)[i][1] for i in 1:3 ])

function prev(c,l::Int=0)
  file = sort!(
    filter(x -> occursin(c,x) && occursin(".dat",x),readdir("./data")),
    by= s->vers(s,"$(c)_v"))[end-l]
  return file[5:end-4] 
end

function get_elapsed(n,cd,parallel)
    ntrials = 1
    t = +Inf
    for _ in 1:ntrials
        GC.gc()
        t = min(t,@belapsed florpi(N=$n,cd=$cd,parallel=$parallel,nbatches=(0,0)))
    end
    return t
end

function run_benchmark(output=false,last_cd=10_000_000,last_cv=3_000_000;nthreads=nothing)

  if nthreads == nothing 
      if Threads.nthreads() != 8 
         error(" Run with julia -t 8")
      end
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
  
  previous_version = prev("cd")
  if previous_version == "$version"
    previous_version = prev("cd",1)
  end
  previous_cd = "./data/cd_v"*previous_version*".dat"
  
  data_cd = readdlm(previous_cd,comments=true,comment_char='#')
  println(" Previous: $previous_cd ")
  
  previous_previous_version = prev("cv")
  if previous_version == "$version"
    previous_version = prev("cv",1)
  end
  previous_cv = "./data/cv_v"*previous_version*".dat"
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
  florpi(N=1000,cd=false,parallel=true);
  for i in 1:ilast_cv
    n = ns[i]
    prev = try data_cv[i,5] catch; 0 end
    t = get_elapsed(n,false,true)
    new_cv[i,5] = t
    println(n," ",t," prev ", prev)
  end
  
  println("Parallel, constant density:")
  florpi(N=1000,cd=true,parallel=true);
  for i in 1:ilast_cd
    n = ns[i]
    prev = try data_cd[i,5] catch; 0 end
    t = get_elapsed(n,true,true)
    new_cd[i,5] = t
    println(n," ",t," prev ", prev)
  end
  
  #
  # Serial
  #
  println("Serial, constant density:")
  florpi(N=1000,cd=true,parallel=false);
  for i in 1:ilast_cd
    n = ns[i]
    prev = try data_cd[i,4] catch; 0 end
    t = get_elapsed(n,true,false)
    new_cd[i,4] = t
    println(n," ",t," prev ", prev)
  end
  
  println("Serial, constant volume:")
  florpi(N=1000,cd=false,parallel=false);
  for i in 1:ilast_cv
    n = ns[i]
    prev = try data_cv[i,4] catch; 0 end
    t = get_elapsed(n,false,false)
    new_cv[i,4] = t
    println(n," ",t," prev ", prev)
  end
  
  if output
    writedlm("./data/cd_v$version.dat",new_cd)
    writedlm("./data/cv_v$version.dat",new_cv)
  end

end
