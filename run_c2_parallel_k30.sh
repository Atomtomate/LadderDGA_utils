#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 41
#SBATCH -p standard96
#SBATCH -A hhp00048

julia --check-bounds=no -p 40 run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_annealing lDGA_ntc_k30_ c2_k30_annealing.jld2 30 > run_k30.out 2> run_k30.err
