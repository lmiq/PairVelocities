import Chemfiles
using CellListMap
using FastPow
using StaticArrays
using Printf

const ε = 0.0441795
const σ = 2*1.64009

function lj_NE(d2,u) 
    d = sqrt(d2)
    @fastpow u += ε*( (σ/d)^12 - 2*(σ/d)^6 )
end

function getcoor(file)
    traj = redirect_stdout(()->Chemfiles.Trajectory(file),devnull)
    frame = Chemfiles.read_step(traj,0)
    Chemfiles.close(traj)
    return reinterpret(reshape,SVector{3,Float64},Chemfiles.positions(frame))
end

function test(file,unit_cell_matrix,correct)
    coordinates = getcoor(file)
    box = Box(unit_cell_matrix, 10.)
    cl = CellList(coordinates,box)
    u = map_pairwise!((x,y,i,j,d2,u) -> lj_NE(d2,u),0.0,box,cl)
    println("$u $correct DIFF = $(@sprintf("%8.5f",100*(u-correct)/u))%")
    nothing
end

unit_cell_matrix = [ 50.      0.      0.
                      0.     50.      0. 
                      0.      0.     50. ]
test("./o1.dcd", unit_cell_matrix, 32230.0218)

unit_cell_matrix = [ 80.      0.      0.
                      0.     70.      0. 
                      0.      0.     50. ]
test("./o2.dcd", unit_cell_matrix, 1093.7230)

unit_cell_matrix = [ 50.      0.     50.
                     50.     50.      0. 
                      0.     50.     50. ]
test("./o3.dcd", unit_cell_matrix, 1724.3208)

unit_cell_matrix = transpose([ 70.7107   0.0      0.0
                               35.3553  61.2372   0.0
                               35.3553  20.4124  57.735 ])
test("./o4.dcd", unit_cell_matrix, 1754.0812)

unit_cell_matrix = transpose([ 70.7107   0.0      0.0
                               35.3553  61.2372   0.0
                               35.3553  20.4124  57.735 ])
test("./o5.dcd", unit_cell_matrix, 1765.1400)

unit_cell_matrix = [ 80.      0.      0.
                      0.     80.      0. 
                      0.      0.     80. ]
test("./o6.dcd", unit_cell_matrix, -158.0470)

unit_cell_matrix = [ 45.      0.      0.
                      0.     45.      0. 
                      0.      0.     45. ]
test("./o7.dcd", unit_cell_matrix, 132542.1195)

unit_cell_matrix = [ 80.      0.     30.
                     30.     80.      0. 
                      0.     40.     80. ]
test("./t1.dcd", unit_cell_matrix, -116.5313)

unit_cell_matrix = [ 50.      0.      0.
                     50.     50.      0. 
                      0.     50.     50. ]
test("./t2.dcd", unit_cell_matrix, 32096.4923)




