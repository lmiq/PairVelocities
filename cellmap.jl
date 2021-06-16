using CellListMap
using DelimitedFiles
using StaticArrays
using BenchmarkTools

#
# florpi
#

@inline dot(x::SVector{3,Float64},y::SVector{3,Float64}) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

function compute_pairwise_mean_cell_lists!(x,y,i,j,d2,hist,velocities, rbins,sides)
    d = x - y
    r = sqrt(d2)
    ibin = searchsortedfirst(rbins, r) - 1
    hist[1][ibin] += 1
    hist[2][ibin] += dot(velocities[i]-velocities[j],d)/r
    return hist
end

function get_pairwise_velocity_radial_mean_cell_lists(
        positions, velocities,
        rbins, r_max,
        boxsize, n
        )
    positions = reshape(reinterpret(SVector{3,Float64},positions),n)
    velocities = reshape(reinterpret(SVector{3,Float64},velocities),n)
    lc = LinkedLists(n)
    box = Box(boxsize, r_max)
    initlists!(positions,box,lc)
    hist = (zeros(Int,length(rbins)-1), zeros(Float64,length(rbins)-1))
    hist = map_pairwise!(
      (x,y,i,j,d2,hist) -> compute_pairwise_mean_cell_lists!(x,y,i,j,d2,hist,velocities, rbins, boxsize),
      hist, positions, box, lc,
    )
    n_pairs = hist[1]
    mean_v_r = hist[2]
    mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0]./n_pairs[n_pairs .> 0]
    return mean_v_r
end

function florpi(;N=100_000,cd=true)
  
  n_halos = N

  println(" constant density = ", cd)
  if cd
    density = 10^4/250^3  # density of the original problem
    boxsize = (n_halos / density)^(1/3)
  else
    boxsize = 250.
  end

  if N == 100_000
      positions = readdlm("./florpi100k_positions.dat")
      velocities = readdlm("./florpi100k_velocities.dat")
      positions = positions'
      velocities = velocities'
  else
      Random.seed!(321)
      positions = boxsize .* rand(Float64, 3, n_halos)
      velocities = rand(Float64, 3, n_halos)
  end
  Lbox = [boxsize,boxsize,boxsize]
  rbins = [0.,2.,4.,6.,8.,10.]

  n = size(positions)[2]
  println(n)
  r_max = maximum(rbins)
  
  mean_v_r = get_pairwise_velocity_radial_mean_cell_lists(
    positions,
    velocities,
    rbins, r_max,
    Lbox, n
  )
  println(mean_v_r)

  println(@btime get_pairwise_velocity_radial_mean_cell_lists(
    $positions,
    $velocities,
    $rbins, $r_max,
    $Lbox, $n
  ))

end

florpi()


