#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 1
#SBATCH -p large96:shared
###SBATCH -A hhp00048
julia PD_gather_ext_l.jl /scratch/projects/hhp00048/lDGA/PD/data_PD_hf_sb_01/U1.0 lDGA_rtc_k10 res_rtc_U1.0_k10 > run.out 2> run.err 
