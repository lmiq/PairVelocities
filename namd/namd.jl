using CellListMap
using PDBTools

#9999 so we do not require two segments
x = [ 100*rand(3) for i in 1:9_999 ]

box = Box([ 50. 10.  0.
             0. 50. 20.
            10.  0. 50.],10)

#cl = CellList(x,box)

p = CellListMap.wrap_to_first.(x,Ref(box.unit_cell.matrix))

atoms = [ Atom(
            index=i,
            chain="A",
            resnum=i,
            residue=i,
            name="NE",
            resname="NE",
            x=p[i][1],
            y=p[i][2],
            z=p[i][3],
            segname="NEON"
          ) for i in 1:length(p) ]

writePDB(atoms,"triclinic1.pdb")








