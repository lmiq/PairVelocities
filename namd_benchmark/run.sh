charmrun=/home/leandro/Downloads/NAMD_2.14_Linux-x86_64-multicore/charmrun
namd2=/home/leandro/Downloads/NAMD_2.14_Linux-x86_64-multicore/namd2
julia=/home/leandro/programs/julia/julia-1.6.2/bin/julia

for sys in 10k 100k 300k; do
    for np in 32 16 8 4 2; do
        echo "$sys $np ..."
        time $charmrun +p$np $namd2 ./ne$sys.namd >& ne$sys\_namd$np.log 
        time $julia -t $np ./ne$sys.jl >& ne$sys\_julia$np.log
        echo "--------"
    done
    echo "$sys 1 ..."
    time $namd2 ./ne$sys.namd >& ne$sys\_namd1.log 
    time $julia ./ne$sys.jl >& ne$sys\_julia1.log
    echo "--------"
done

