#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 96
#SBATCH -p standard96
#SBATCH -A hhp00048

julia --check-bounds=no -p 96 run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_annealing lDGA_ntc_k10_ c2_k10_annealing.jld2 10 > run_p96_k10.out 2> run_p96_k10.err
