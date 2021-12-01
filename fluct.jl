using CellListMap, UnicodePlots
CellListMap.Examples.florpi()
t = Float64[]
for i in 1:30
    e = @elapsed CellListMap.Examples.florpi(N=2_000_000,cd=false,nbatches=(8,2^13))
    @show i, e
    push!(t,e)
end
histogram(t)
