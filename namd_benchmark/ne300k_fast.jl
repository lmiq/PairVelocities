include("./simulate_fast.jl")

side = 144.05129897602086

p = Params(
    dt = 1.0,
    x0 = getcoor("./ne300k_initial.pdb"),
    box = Box([ side, side, side ], 12.),
    trajfile = "ne300k_traj.xyz"
)
@time simulate_fast(p)
