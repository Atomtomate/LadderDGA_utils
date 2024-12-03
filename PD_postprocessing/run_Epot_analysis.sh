#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 96
#SBATCH -p large96
#SBATCH -A hhp00048
julia -p 95 /scratch/projects/hhp00048/codes/LadderDGA_utils/E_pot_analysis.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_sb_01/U1.0 lDGA_rtc_k10 Epot_analysis_rtc_k10_testU1.jld2 1 rtc > run_rtc_10.out 2> run_rtc_10.err
