using DelimitedFiles, Plots, Plots.Measures, LaTeXStrings

data = readdlm("./cv.dat",comments=true,comment_char='#')
default(fontfamily="Computer Modern",label="",linewidth=2,framestyle=:box,legend=:topleft,margin=5mm)
plot(data[:,1],data[:,2],label="Serial Python/halotools")
plot!(data[:,1],data[:,3],label="4 cores Python/halotools")
plot!(data[:,1],data[:,4],label="Serial Julia/CellListMap.jl")
plot!(data[:,1],data[:,5],label="4 cores Julia/CellListMap.jl")
plot!(xlabel="Number or particles",ylabel="time / s")
plot!(title=L"\textrm{Constant\ volume - V=250^3;\ cutoff = 10}")
savefig("./cv.png")

data = readdlm("./cd.dat",comments=true,comment_char='#')
default(fontfamily="Computer Modern",label="",linewidth=2,framestyle=:box,legend=:topleft,margin=5mm)
plot(data[:,1],data[:,2],label="Serial Python/halotools")
plot!(data[:,1],data[:,3],label="4 cores Python/halotools")
plot!(data[:,1],data[:,4],label="Serial Julia/CellListMap.jl")
plot!(data[:,1],data[:,5],label="4 cores Julia/CellListMap.jl")
plot!(xlabel="Number or particles",ylabel="time / s")
plot!(title=L"\textrm{Constant\ density - \rho=(10^5/250^3) N/V;\ cutoff = 10}")
savefig("./cd.png")


