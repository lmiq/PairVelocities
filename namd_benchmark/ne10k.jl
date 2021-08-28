include("./simulate.jl")

side =  46.35999999999999
p = Params(
    dt = 1.0,
    x0 = getcoor("./ne10k_initial.pdb"),
    box = Box([ side, side, side ], 12.0),
    trajfile = "ne10k_traj.xyz"
)
@time simulate(p)
