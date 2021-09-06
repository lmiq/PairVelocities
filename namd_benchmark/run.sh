charmrun=/home/leandro/Downloads/NAMD_2.14_Linux-x86_64-multicore/charmrun
namd2=/home/leandro/Downloads/NAMD_2.14_Linux-x86_64-multicore/namd2
julia=/home/leandro/programs/julia/julia-1.6.2/bin/julia

for sys in 10k 100k 300k; do
    echo "$sys 1 ..."
    time $namd2 ./ne$sys.namd >& ne$sys\_namd1.log 
    time $julia --check-bounds=no ./ne$sys.jl >& ne$sys\_julia1.log
    echo "--------"
done
for np in 2 4 8 16 32; do
    for sys in 10k 100k 300k; do
        echo "$sys $np ..."
        time $charmrun +p$np $namd2 ./ne$sys.namd >& ne$sys\_namd$np.log 
        time $julia --check-bounds=no -t $np ./ne$sys.jl >& ne$sys\_julia$np.log
        echo "--------"
    done
done

