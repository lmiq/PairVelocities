using BenchmarkTools, Test
using NearestNeighbors
using CellListMap, StaticArrays
using LinearAlgebra

function nl_NN(x,y,r)
    balltree = BallTree(x)
    return inrange(balltree,y,r,true)
end

function nl_CL(x,y,r;parallel=true,autoswap=false)
    box = Box(limits(x,y),r)
    cl = CellList(x,y,box,parallel=parallel,autoswap=autoswap)
    return CellListMap.neighbourlist(box,cl,parallel=parallel)
end

function compare_result(list_CL,list_NN)
    if length(list_CL) <= 100
        for (i,list) in pairs(list_NN)
            cl = filter(tup -> tup[2] == i, list_CL)
            length(cl) != length(list) &&  return false
            for j in list
               length(findall(tup -> tup[1] == j, cl)) != 1 && return false
            end
        end
    else
        npairs = 0
        for i in list_NN
            npairs += length(list_NN[i])
        end
        npairs != length(list_CL) && return false
    end
    return true
end

function naive(x,y,cutoff)
  pair_list = Int[]
  for vx in x
    for (i,vy) in pairs(y) 
      if norm(vx - vy) <= cutoff  
        push!(pair_list,i)
      end
    end
  end
  pair_list
end

function nl_run()

  r = 0.05

  for N1 in [1, 10, 100, 1_000, 10_000, 100_000]
    for N2 in [10^6]

       x = [ rand(SVector{3,Float64}) for i in 1:N1 ]
       y = [ rand(SVector{3,Float64}) for i in 1:N2 ]
       
       list_CL = nl_CL(x,y,r,parallel=true)
       list_NN = nl_NN(x,y,r)
       
       println("----------------------------------")
       print("N1 = $N1 ; N2 = $N2 ; PASS TEST = ")
       println(compare_result(list_CL,list_NN))
       print("nl (x,y): "); @btime nl_NN($x,$y,$r) samples=1 
       print("nl (y,x): "); @btime nl_NN($y,$x,$r) samples=1 
       print("cl serial (x,y): "); @btime nl_CL($x,$y,$r,parallel=false, autoswap=false) samples=1
       print("cl parallel (x,y): "); @btime nl_CL($x,$y,$r,parallel=true, autoswap=false) samples=1
       print("cl serial (y,x): "); @btime nl_CL($y,$x,$r,parallel=false, autoswap=false) samples=1
       print("cl parallel (y,x): "); @btime nl_CL($y,$x,$r,parallel=true, autoswap=false) samples=1

    end
  end

  return 
end







