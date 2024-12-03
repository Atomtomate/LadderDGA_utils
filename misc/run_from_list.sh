while read p; do
    vars=( $p )
    cd "$p/lDGA_julia"
    ls
    sed -i '9s/.*/Nk = [10,20]/' config.toml || true
    sbatch lDGA_j.sh || true
    cd "../.."
done < todo.txt
