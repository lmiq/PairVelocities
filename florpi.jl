import Pkg
using CellListMap
using Plots, Plots.Measures
using DelimitedFiles
using LaTeXStrings

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

version=filter(x-> x.second.name == "CellListMap", Pkg.dependencies()) |> x -> first(x)[2].version

println(" Version: v$version ")

#
# CD
#
previous_cd = sort!(filter(x -> occursin("cd",x) && occursin(".dat",x),readdir("./data")))[end]
previous_version = previous_cd[5:9] 
if previous_version == "$version"
  previous_cd = sort!(filter(x -> occursin("cd",x) && occursin(".dat",x),readdir("./data")))[end-1]
end
previous_cd = "./data/"*previous_cd

data_cd = readdlm(previous_cd,comments=true,comment_char='#')
println(" Previous: $previous_cd ")

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

println("Parallel, constant density:")
CellListMap.florpi(N=1000,cd=true,parallel=true);
times2 = Float64[]
for i in 1:length(ns)
  n = ns[i]
  prev = data_cd[i,5]
  t = @elapsed CellListMap.florpi(N=n,cd=true,parallel=true);
  push!(times2,t)
  println(n," ",t," prev ", prev)
end

data_cd[:,4] .= times1
data_cd[:,5] .= times2

default(fontfamily="Computer Modern",label="",linewidth=2,framestyle=:box,legend=:topleft,margin=5mm)
data = data_cd
plot(data[:,1],data[:,2],label="Serial/halotools v0.7")
plot!(data[:,1],data[:,3],label="4 cores/halotools v0.7")
plot!(data[:,1],data[:,4],label="Serial/CellListMap.jl $version")
plot!(data[:,1],data[:,5],label="4 cores/CellListMap.jl $version")
plot!(xlabel="Number or particles",ylabel="time / s")
plot!(title=L"\textrm{Constant\ density - \rho=(10^5/250^3) N/V;\ cutoff = 10}")
savefig("./data/cd_v$version.png")
writedlm("./data/cd_v$version.dat",data)

println("created ./data/cd_v$version.png")
println("wrote ./data/cd_v$version.dat")

#
# CV
#
previous_cv = sort!(filter(x -> occursin("cv",x) && occursin(".dat",x),readdir("./data")))[end]
previous_version = previous_cv[5:9] 
if previous_version == "$version"
  previous_cv = sort!(filter(x -> occursin("cv",x) && occursin(".dat",x),readdir("./data")))[end-1]
end
previous_cv = "./data/"*previous_cv
data_cv = readdlm(previous_cv,comments=true,comment_char='#')
println(" Previous: $previous_cv ")

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

println("Parallel, constant volume:")
CellListMap.florpi(N=1000,cd=false,parallel=true);
times4 = Float64[]
for i in 1:15
  n = ns[i]
  prev = data_cv[i,5]
  t = @elapsed CellListMap.florpi(N=n,cd=false,parallel=true);
  push!(times4,t)
  println(n," ",t," prev ", prev)
end


data_cv[:,4] .= times3
data_cv[:,5] .= times4

data = data_cv
plot(data[:,1],data[:,2],label="Serial/halotools v0.7")
plot!(data[:,1],data[:,3],label="4 cores/halotools v0.7")
plot!(data[:,1],data[:,4],label="Serial/CellListMap.jl $version")
plot!(data[:,1],data[:,5],label="4 cores/CellListMap.jl $version")
plot!(xlabel="Number or particles",ylabel="time / s")
plot!(title=L"\textrm{Constant\ volume - V=250^3;\ cutoff = 10}")
savefig("./data/cv_v$version.png")
writedlm("./data/cv_v$version.dat",data)

println("created ./data/cv_v$version.png")
println("wrote ./data/cv_v$version.dat")



