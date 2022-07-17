#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --ntasks 2
#SBATCH -p standard96:test
#SBATCH -A hhp00048

julia --check-bounds=no -p 2 run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_annealing/b12.3_U3.5 lDGA_ntc_k10_ c2_k10_testing_single.jld2 10 > run_single.out 2> run_single.err
