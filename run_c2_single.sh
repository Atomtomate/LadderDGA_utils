#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 1
#SBATCH -p large96:shared
#SBATCH -A hhp00048

julia run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_annealing/b14.0_U2.75 lDGA_ntc_k10_ c2_k10_annealing_single_2.jld2 10 > run_1_2.out 2> run_1_2.err
