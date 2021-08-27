using PDBTools

for n in [ 10000, 100000, 300000 ]

    x = [ 300*rand(3) for i in 1:n ]

    d = 10_000/(46.36^3)
    v = n / d
    l = v^(1/3)
    println(l)
    
    atoms = [ Atom(
                index=i,
                chain="A",
                resnum = i%10000,
                residue = i%10000,
                name="NE",
                resname="NE",
                x=x[i][1],
                y=x[i][2],
                z=x[i][3],
                segname="N$(div(i,10000)+1)"
              ) for i in 1:n ]
    
    writePDB(atoms,"$n.pdb")

end







