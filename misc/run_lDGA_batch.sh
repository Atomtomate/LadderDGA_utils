#!/bin/bash


python run_lDGA.py /scratch/projects/hhp00048/lDGA/cuprates/ED_New/beta030 /scratch/projects/hhp00048/codes/LadderDGA.jl/scripts/script_02b_ldm_sc.jl 200 res_TailComp_tcEoM_Nk10.jld2 1 
python run_lDGA.py /scratch/projects/hhp00048/lDGA/cuprates/ED_New/beta030 /scratch/projects/hhp00048/codes/LadderDGA.jl/scripts/script_02b_ldm_sc.jl 200 res_TailComp_tcExpStep08_Nk10.jld2 2 0.8 
python run_lDGA.py /scratch/projects/hhp00048/lDGA/cuprates/ED_New/beta030 /scratch/projects/hhp00048/codes/LadderDGA.jl/scripts/script_02b_ldm_sc.jl 200 res_TailComp_tcExpStep02_Nk10.jld2 2 0.2 
python run_lDGA.py /scratch/projects/hhp00048/lDGA/cuprates/ED_New/beta030 /scratch/projects/hhp00048/codes/LadderDGA.jl/scripts/script_02b_ldm_sc.jl 200 res_TailComp_tcExpStep001_Nk10.jld2 2 0.01 
python run_lDGA.py /scratch/projects/hhp00048/lDGA/cuprates/ED_New/beta030 /scratch/projects/hhp00048/codes/LadderDGA.jl/scripts/script_02b_ldm_sc.jl 200 res_TailComp_tcFull_Nk10.jld2 3 
