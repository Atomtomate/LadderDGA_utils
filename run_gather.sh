#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --ntasks 1
#SBATCH -p standard96:test
#SBATCH -A hhp00048

julia PD_gather.jl /scratch/projects/hhp00048/lDGA/PD/data_3Dsc_n1.0_c2Grid lDGA_ntc_k20 c2Grid
