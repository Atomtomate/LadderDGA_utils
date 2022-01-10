"""
    Description:
Helper functions to read through directory tree containing lDGA results.
Data directory are expected to contain a jld2 and config.toml file.
The result jld2 file is expected to contain _k<kn>_ in its name, where <kn>
is the number of k points of the lDGA grid.
#TODO: in the future this may be read from the contents of the result file, probably when more k grid types are added

"""
function walk_lDGA_dir(dir::String, fname_pre::String, transform!::Function, res; IOh=stdout)
    dpw = displaysize(stdout)
    println(IOh, "Walking through $dir")
    for (root, dirs, files) in walkdir(dir)
        if "config.toml" in files
            flist = filter(x->startswith(x,fname_pre), files)
            print(IOh, "\r$(repeat(" ", dpw[2]))")
            for file in flist
                m = match(r"_k(?<kn>\d+)_",file)
                if m !== nothing
                    kn = parse(Int,m[:kn])     
                    f = joinpath(root, file)
                    transform!(res, f, root)
                end
            end
        end
    end
end


"""
    valid_paths(dir::String, fname_pre::String)

This is an example how to use `walk_lDGA_dir`. It will return all
subdirectories of `dir` containing files starting with `fname_pre`
which also fulfill the requirements of `walk_lDGA_dir`.
"""
function valid_paths(dir::String, fname_pre::String)
    paths = String[]
    tf(res, file, root) = push!(res, root)
    walk_lDGA_dir(dir, fname_pre, tf, paths)
    return unique(paths)
end
