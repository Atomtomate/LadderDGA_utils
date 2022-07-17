#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --ntasks 1
#SBATCH -p large96:shared
#SBATCH -A hhp00048

julia PD_gather_kConv.jl /scratch/projects/hhp00048/lDGA/PD/data_3Dsc_n1.0_c2Grid lDGA_ntc_kConv c2Grid_kConv
