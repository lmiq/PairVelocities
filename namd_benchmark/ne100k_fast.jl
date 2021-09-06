include("./simulate_fast.jl")

side = 99.8795922298781 
p = Params(
    dt = 1.0,
    x0 = getcoor("./ne100k_initial.pdb"),
    box = Box([ side, side, side ], 12.),
    trajfile = "ne100k_traj.xyz"
)
@time simulate_fast(p)

