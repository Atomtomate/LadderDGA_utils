#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --ntasks 1
#SBATCH -p large96:shared
#SBATCH -A hhp00048

julia collect_c2_curve_results.jl /scratch/projects/hhp00048/codes/scripts/LadderDGA_utils/c2_data_all c2_all.jld2
