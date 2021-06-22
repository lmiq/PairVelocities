using CellListMap

function readdata()
  println("Reading positions...")
  file = open("~/Downloads/data_enrique/pos.dat","r") 
  positions = SVector{3,Float64}[]
  for line in eachline(file)
    push!(positions,SVector{3,Float64}(parse.(Float64,split(line))))
  end
  close(file)
  
  println("Reading velocities...")
  file = open("~/Downloads/data_enrique/vel.dat","r") 
  velocities= SVector{3,Float64}[]
  for line in eachline(file)
    push!(velocities,SVector{3,Float64}(parse.(Float64,split(line))))
  end
  close(file)
  return positions,velocities
end

function compute_mean(
  positions,velocities;
  parallel=true,
  show_progress=true,
  lcell=2
)

  @inline dot(x::SVector{3,Float64},y::SVector{3,Float64}) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

  function compute_pairwise_mean_cell_lists!(x,y,i,j,d2,hist,velocities,rbins,sides)
    d = x - y
    r = sqrt(d2)
    ibin = searchsortedfirst(rbins, r) - 1
    hist[1][ibin] += 1
    hist[2][ibin] += dot(velocities[i]-velocities[j],d)/r
    return hist
  end

  function reduce_hist(hist,hist_threaded)
    hist = hist_threaded[1]
    for i in 2:Threads.nthreads()
      hist[1] .+= hist_threaded[i][1]
      hist[2] .+= hist_threaded[i][2]
    end
    return hist
  end

  N = 8_000_000
  pos = @view(positions[1:N])
  vel = @view(velocities[1:N])
  
  println("Building lists...")
  lbox = [2000,2000,2000]
  cutoff = 200
  box = Box(lbox,cutoff,lcell=lcell)
  cl = CellList(pos,box)

  rbins = 0.:20.:200.
  hist = (zeros(Int,length(rbins)-1), zeros(Float64,length(rbins)-1))      

  println("Computing histogram...")
  hist = map_pairwise!(
    (x,y,i,j,d2,hist) -> compute_pairwise_mean_cell_lists!(
      x,y,i,j,d2,hist,vel,rbins,lbox
    ),
    hist, box, cl,
    reduce=reduce_hist,
    parallel=parallel,
    show_progress=show_progress
  )
  n_pairs = hist[1]
  mean_v_r = hist[2]
  mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0]./n_pairs[n_pairs .> 0]

  return mean_v_r

end








