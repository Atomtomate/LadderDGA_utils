"""
    Description:
Combines multiple similar jld2 files into one. ARG1 is the target file, ARG2 ... ARGN are the source files.

    Example:
julia res_combined res_U1.0 res_U1.1
"""

using JLD2

function combine_data!(destination::JLD2.JLDFile{JLD2.MmapIO}, source::JLD2.JLDFile{JLD2.MmapIO}; prefix::String = "")
    ks = prefix == "" ? keys(source) : keys(source[prefix])
    for k in  ks
        pk = prefix == "" ? k : prefix*"/"*k
        if typeof(source[pk]) <: JLD2.Group
            combine_data!(destination, source, prefix=pk)
        else
            destination[pk] = source[pk]
        end
    end
end

length(ARGS) < 2 && error("Not enough arguments. Usage: julia combine_jld2.jl DESTINATION SOURCE1 SOURCE2 ...")

dest = ARGS[1]
sources = ARGS[2:end]
isfile(dest*".jld2") && error("Output file already exists!")

written_anything = false
jldopen(dest*".jld2", "w") do fnew
for source in sources
    !isfile(source*".jld2") && (println("File $(source).jld2 not found. Skipping") ; continue)
    global written_anything = true
    println("Copying $(source).jld2 to $(dest).jld2")
    jldopen(source*".jld2","r") do fold
        combine_data!(fnew, fold)
    end
end
end

!written_anything && rm(dest*".jld2")
