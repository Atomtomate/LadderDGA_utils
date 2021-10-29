#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --ntasks 1
#SBATCH -p large96:shared
#SBATCH -A hhp00048

julia PD_gather_2_ext_l.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_sb_01/U1.0 lDGA_rtc_k10_ c2_k10_U1.0.jld2 10 > run.out 2> run.err
