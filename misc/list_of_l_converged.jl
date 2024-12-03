"""
    Description:
Read through directory tree given by ARG1 and search for chi_lDGA output starting with ARG2. 
Write of convergence and value of new lambda to ARG3.

    Example:
julia list_of_lconverged.jl /scratch/projects/hhp00048/lDGA/PD/data lDGA_ res
"""

using Pkg
Pkg.activate("/scratch/projects/hhp00048/codes/LadderDGA.jl",io=devnull)
using LadderDGA
using JLD2
using CSV

include("helpers/run_lDGA_dir.jl")

dir = ARGS[1]
fname_pre = ARGS[2]
out_fname = ARGS[3]
println(ARGS)
flush(stdout)

function gather_lambda(fname)
    res = []
    jldopen(fname, "r") do f
        #res = (f["λsp"], f["λspch"]..., )
        res = (f["λsp"], f["λspch"].zero..., f["λspch"].f_converged)
    end
    return res
end

paths = valid_paths(dir, fname_pre)
println(paths, "dbg")
res_paths = String[]
for i in 1:length(paths)
    path = paths[i]
    #TODO: config from input
    in_files = filter(x->startswith(x,fname_pre), readdir(path))
    in_files = filter(x->endswith(x,".jld2"), in_files)
    for in_file in in_files
        println("combining $(joinpath(path,in_file))")
        flush(stdout)
        push!(res_paths, joinpath(path,in_file))
    end
end
println("constructed list of files: ")
println(res_paths)
flush(stdout)
res = map(gather_lambda, res_paths)
println("------")
#println(res)
open(out_fname, "w") do f
    for row in res
        println(row)
        write(f, "$(row[1]) $(row[2]) $(row[3]) $(row[4])\n")
      end
end
