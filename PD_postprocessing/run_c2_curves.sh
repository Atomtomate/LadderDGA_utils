#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 96
#SBATCH -p large96
###SBATCH -A hhp00048
julia -p 95 run_c2_curves.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_sb_01 lDGA_rtc_k10 c2_ntc_10.jld2 10 > run_n_10.out 2> run_n_10.err
